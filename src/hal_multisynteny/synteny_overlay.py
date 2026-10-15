"""Experimental all-edge halSynteny chains clipped to existing child blocks."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import shutil
import subprocess
from bisect import bisect_left
from collections import defaultdict
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from .extract import (
    EdgeExtractionResult,
    ExtractionError,
    HalSyntenyLeafInitializer,
    HalSyntenyLeafResult,
    _leaf_blocks_from_rows,
    run_halsynteny,
)
from .hal import _parse_parent_map
from .models import EdgeBlockState, EdgeMappingRun, UnmappedEdgeEvidence
from .psl_extract import STATUS_POLICY_VERSION, edge_block_state
from .source_diagnostic import PSL_PARSER_VERSION, PSLComponent, parse_psl_record
from .utils import sha256_file

BACKEND = "hal-synteny-overlay"
OVERLAY_SCHEMA_VERSION = "1"
OVERLAY_STATUS_POLICY_VERSION = "1"
WARNING = (
    "Experimental halSynteny overlay: chain alternatives are unresolved local mapping paths, "
    "not biological copies. Uncovered sequence is not a confirmed deletion. "
    "halSynteny transitivity across tree levels has not been established."
)


@dataclass(frozen=True, slots=True)
class HalSyntenyChain:
    synteny_chain_id: str
    components: tuple[PSLComponent, ...]
    raw_line_numbers: tuple[int, ...]


def parse_synteny_chains(raw: str) -> tuple[HalSyntenyChain, ...]:
    """Use the audited DNA PSL normalizer; identical mappings share one chain.

    Preserve all raw rows separately. Line numbers never participate in IDs.
    The existing leaf-family constructor remains compatible and unchanged.
    """
    chains = {}
    for line_number, line in enumerate(raw.splitlines(), 1):
        if not line.strip() or line.lstrip().startswith(("#", "psLayout", "match", "---")):
            continue
        _, record = parse_psl_record("synteny\t" + line, line_number)
        key = tuple(
            (
                c.q_name,
                c.q_size,
                c.q_start,
                c.q_end,
                c.t_name,
                c.t_start,
                c.t_end,
                c.relative_orientation,
            )
            for c in record.components
        )
        chain_id = (
            "HSC_" + hashlib.sha256(json.dumps(key, separators=(",", ":")).encode()).hexdigest()
        )
        if chain_id in chains:
            previous = chains[chain_id]
            representative = min((previous.components, record.components), key=lambda cs: cs[0].raw)
            chains[chain_id] = replace(
                previous,
                components=representative,
                raw_line_numbers=previous.raw_line_numbers + (line_number,),
            )
        else:
            chains[chain_id] = HalSyntenyChain(chain_id, record.components, (line_number,))
    return tuple(chains[key] for key in sorted(chains))


def tsv(rows, fields):
    handle = io.StringIO()
    writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue()


class HalSyntenyOverlayLeafInitializer(HalSyntenyLeafInitializer):
    """Keep legacy leaf families/outputs, using the audited component normalizer."""

    def probe_version(self):
        if not self.executable_version:
            self.executable_version = _tool_version(self.executable, "halSynteny")
        return self.executable_version

    def initialize(self, child_node, parent_node, query_genome, target_genome):
        version = self.probe_version()
        raw, command, stdout, stderr = run_halsynteny(self, query_genome, target_genome)
        rows = []
        for chain in parse_synteny_chains(raw):
            cs = chain.components
            query = tuple((c.q_start, c.q_end) for c in cs)
            rows.append(
                {
                    "fields": tuple(cs[0].raw.split("\t")[1:]),
                    "qkey": (cs[0].q_name, cs[0].q_size, query),
                    "qname": cs[0].q_name,
                    "tname": cs[0].t_name,
                    "q_intervals": query,
                    "t_intervals": tuple((c.t_start, c.t_end) for c in cs),
                    "strand": cs[0].relative_orientation,
                }
            )
        blocks, runs = _leaf_blocks_from_rows(rows, child_node, parent_node, command)
        return HalSyntenyLeafResult(
            blocks,
            runs,
            raw,
            command,
            stdout,
            stderr,
            {"halSynteny": version},
            self.executable,
            version,
        )


def overlay_node_blocks_on_chains(blocks, chains, *, child_node, parent_node, command=()):
    """Project each exact intersection; input block identity is never replaced."""
    index = defaultdict(list)
    for chain in chains:
        for component in chain.components:
            index[component.q_name].append((component.q_start, component.q_end, chain, component))
    for entries in index.values():
        entries.sort(key=lambda item: (item[0], item[1], item[2].synteny_chain_id))
    starts = {chrom: [item[0] for item in entries] for chrom, entries in index.items()}
    max_lengths = {
        chrom: max(item[1] - item[0] for item in entries) for chrom, entries in index.items()
    }
    runs, unmapped, uncovered, states, observations = [], [], [], [], []
    for block in sorted(blocks, key=lambda b: (b.chrom, b.start, b.end, b.block_id)):
        entries = index.get(block.chrom, [])
        lo = bisect_left(starts.get(block.chrom, []), block.start - max_lengths.get(block.chrom, 0))
        hi = bisect_left(starts.get(block.chrom, []), block.end)
        hits = [
            (chain, c, max(block.start, c.q_start), min(block.end, c.q_end))
            for _, end, chain, c in entries[lo:hi]
            if end > block.start
        ]
        alternatives = any(
            first[0].synteny_chain_id != second[0].synteny_chain_id
            and first[2] < second[3]
            and second[2] < first[3]
            for i, first in enumerate(hits)
            for second in hits[i + 1 :]
        )
        inherited = edge_block_state(block, 1)
        raw = "ambiguous" if alternatives else "unique" if hits else "unaligned"
        effective = (
            (inherited.inherited_barrier_status or ("duplicated" if alternatives else "unique"))
            if hits
            else "unaligned"
        )
        state = EdgeBlockState(
            block.block_id,
            raw,
            effective,
            inherited.inherited_barrier_status,
            inherited.inherited_barrier_reason,
            block.mapping_status,
            block.classification,
        )
        states.append(state)
        gap_index = len(uncovered)
        cursor = block.start
        for start, end in sorted((a, b) for _, _, a, b in hits):
            if start > cursor:
                uncovered.append(
                    {
                        "child_block_id": block.block_id,
                        "child_chrom": block.chrom,
                        "child_start": cursor,
                        "child_end": start,
                        "reason": "no_chain_component",
                    }
                )
            cursor = max(cursor, end)
        if cursor < block.end:
            uncovered.append(
                {
                    "child_block_id": block.block_id,
                    "child_chrom": block.chrom,
                    "child_start": cursor,
                    "child_end": block.end,
                    "reason": "no_chain_component",
                }
            )
        if not hits:
            unmapped.append(
                UnmappedEdgeEvidence(
                    child_node,
                    parent_node,
                    block.block_id,
                    block.chrom,
                    block.start,
                    block.end,
                    "unaligned",
                    "no halSynteny component overlap; not a deletion call",
                    f"{child_node}:{block.block_id}:unaligned",
                    "halSynteny",
                    json.dumps(command),
                    raw,
                    state.inherited_barrier_status,
                    state.inherited_barrier_reason,
                )
            )
        for chain, c, a, b in hits:
            if c.relative_orientation == "+":
                pstart, pend = c.t_start + a - c.q_start, c.t_start + b - c.q_start
            else:
                pstart, pend = c.t_end - (b - c.q_start), c.t_end - (a - c.q_start)
            if b - a != pend - pstart:
                raise ExtractionError("overlay projection changed component length")
            component_id = f"{chain.synteny_chain_id}:component{c.component_order:06d}"
            occurrence = f"{child_node}:{block.block_id}:overlay-input"
            runs.append(
                EdgeMappingRun(
                    child_node,
                    block.block_id,
                    occurrence,
                    block.chrom,
                    a,
                    b,
                    parent_node,
                    c.t_name,
                    pstart,
                    pend,
                    c.relative_orientation,
                    "",
                    effective,
                    f"{occurrence}:{component_id}",
                    "halSynteny",
                    json.dumps(command),
                    chain.synteny_chain_id,
                    c.query_strand,
                    c.target_strand_explicit or "",
                    c.raw_strand_field,
                    c.psl_line_number,
                    c.component_order,
                    raw,
                    state.inherited_barrier_status,
                    state.inherited_barrier_reason,
                    chain.synteny_chain_id,
                    component_id,
                )
            )
        uncovered_bp = sum(row["child_end"] - row["child_start"] for row in uncovered[gap_index:])
        observations.append(
            {
                "child_block_id": block.block_id,
                "chain_count": len({ch.synteny_chain_id for ch, _, _, _ in hits}),
                "component_run_count": len(hits),
                "overlapping_alternatives": alternatives,
                "target_chromosome_count": len({c.t_name for _, c, _, _ in hits}),
                "relative_orientation_count": len({c.relative_orientation for _, c, _, _ in hits}),
                "covered_source_bases": block.end - block.start - uncovered_bp,
                "uncovered_source_bases": uncovered_bp,
            }
        )
    runs.sort(
        key=lambda r: (
            r.parent_chrom,
            r.parent_start,
            r.parent_end,
            r.child_block_id,
            r.synteny_chain_id,
            r.component_id,
        )
    )
    return tuple(runs), tuple(unmapped), tuple(states), uncovered, observations


def _tool_version(tool, name):
    help_run = subprocess.run((tool, "--help"), check=False, capture_output=True, text=True)
    lines = (help_run.stdout + help_run.stderr).splitlines()
    marker = f"{name} v"
    return next(
        (line[line.index(marker) :] for line in lines if marker in line),
        next(iter(lines), "unknown"),
    )


class HalSyntenyEdgeExtractor:
    backend = BACKEND

    def __init__(
        self,
        hal_path,
        *,
        min_block_size=50,
        max_anchor_distance=1000,
        executable=None,
        hal_stats=None,
        timeout_seconds=None,
    ):
        if min_block_size < 1 or max_anchor_distance < 0:
            raise ExtractionError("minBlockSize must be positive and maxAnchorDistance nonnegative")
        self.runner = HalSyntenyLeafInitializer(
            hal_path,
            min_block_size=min_block_size,
            max_anchor_distance=max_anchor_distance,
            executable=executable,
        )
        self.runner.hal_path = self.runner.hal_path.resolve()
        self.runner.execution_timeout = timeout_seconds
        self.hal_checksum = self.runner.hal_checksum
        stats = hal_stats or shutil.which("halStats")
        if not stats:
            raise ExtractionError("missing HAL tool: halStats")
        self.hal_stats = str(Path(stats).resolve())
        tree = subprocess.run(
            (self.hal_stats, str(self.runner.hal_path), "--tree"),
            check=False,
            capture_output=True,
            text=True,
        )
        if tree.returncode:
            raise ExtractionError(f"halStats --tree failed: {tree.stderr}")
        self.parent_map = _parse_parent_map(tree.stdout)
        self.runner.executable_version = _tool_version(self.runner.executable, "halSynteny")
        self.stats_version = _tool_version(self.hal_stats, "halStats")

    def identity(self):
        return {
            "backend": BACKEND,
            "hal_sha256": self.hal_checksum,
            "parent_map": self.parent_map,
            "parser_version": PSL_PARSER_VERSION,
            "schema_version": OVERLAY_SCHEMA_VERSION,
            "inherited_status_policy_version": STATUS_POLICY_VERSION,
            "status_policy_version": OVERLAY_STATUS_POLICY_VERSION,
            "minBlockSize": self.runner.min_block_size,
            "maxAnchorDistance": self.runner.max_anchor_distance,
            "executable": self.runner.executable,
            "executable_sha256": self.runner.executable_sha256,
            "version": self.runner.executable_version,
            "hal_stats": self.hal_stats,
            "hal_stats_sha256": sha256_file(self.hal_stats),
            "hal_stats_version": self.stats_version,
        }

    def extract(
        self, child_node_id, parent_node_id, child_hal_genome, parent_hal_genome, child_blocks
    ):
        self._direct(child_hal_genome, parent_hal_genome)
        blocks = tuple(child_blocks)
        # Chains describe the whole edge, independently of current block coverage.
        raw, command, stdout, stderr = run_halsynteny(
            self.runner, child_hal_genome, parent_hal_genome
        )
        return self.from_raw(
            child_node_id,
            parent_node_id,
            child_hal_genome,
            parent_hal_genome,
            blocks,
            raw,
            command,
            stdout,
            stderr,
        )

    def _direct(self, child, parent):
        if self.parent_map.get(child) != parent:
            raise ExtractionError(f"not a direct HAL child-to-parent edge: {child}->{parent}")

    def from_raw(
        self,
        child_node,
        parent_node,
        child_genome,
        parent_genome,
        blocks,
        raw,
        command=(),
        stdout="",
        stderr="",
    ):
        self._direct(child_genome, parent_genome)
        blocks = tuple(blocks)
        if len({b.block_id for b in blocks}) != len(blocks) or any(
            b.node != child_node for b in blocks
        ):
            raise ExtractionError("duplicate block IDs or wrong child node")
        chains = parse_synteny_chains(raw)
        runs, unmapped, states, uncovered, observations = overlay_node_blocks_on_chains(
            blocks, chains, child_node=child_node, parent_node=parent_node, command=command
        )
        metrics = {
            "input_child_block_count": len(blocks),
            "block_outcome_denominator": len(blocks),
            "hal_synteny_chain_count": len(chains),
            "component_count": sum(len(ch.components) for ch in chains),
            "blocks_no_chain_overlap": sum(o["chain_count"] == 0 for o in observations),
            "blocks_one_chain": sum(o["chain_count"] == 1 for o in observations),
            "blocks_source_disjoint_chains": sum(
                o["chain_count"] > 1 and not o["overlapping_alternatives"] for o in observations
            ),
            "blocks_overlapping_alternative_chains": sum(
                o["overlapping_alternatives"] for o in observations
            ),
            "blocks_multiple_target_chromosomes": sum(
                o["target_chromosome_count"] > 1 for o in observations
            ),
            "blocks_multiple_relative_orientations": sum(
                o["relative_orientation_count"] > 1 for o in observations
            ),
            "covered_source_bases": sum(o["covered_source_bases"] for o in observations),
            "uncovered_source_bases": sum(o["uncovered_source_bases"] for o in observations),
            "inherited_ambiguity_blocks": sum(
                s.inherited_barrier_status == "ambiguous" for s in states
            ),
            "inherited_duplication_blocks": sum(
                s.inherited_barrier_status == "duplicated" for s in states
            ),
            "raw_mapping_outcome_counts": {
                status: sum(s.raw_mapping_outcome == status for s in states)
                for status in ("unique", "ambiguous", "unaligned")
            },
            "effective_status_counts": {
                status: sum(s.effective_status == status for s in states)
                for status in ("unique", "duplicated", "ambiguous", "unaligned")
            },
        }
        if (
            sum(metrics["raw_mapping_outcome_counts"].values()) != len(blocks)
            or sum(metrics["effective_status_counts"].values()) != len(blocks)
            or metrics["covered_source_bases"] + metrics["uncovered_source_bases"]
            != sum(b.end - b.start for b in blocks)
        ):
            raise ExtractionError("overlay block outcomes or source coverage do not reconcile")
        chain_rows = [
            {
                "synteny_chain_id": chain.synteny_chain_id,
                "component_id": f"{chain.synteny_chain_id}:component{c.component_order:06d}",
                "query_chrom": c.q_name,
                "query_start": c.q_start,
                "query_end": c.q_end,
                "target_chrom": c.t_name,
                "target_start": c.t_start,
                "target_end": c.t_end,
                "query_strand": c.query_strand,
                "target_strand_explicit": c.target_strand_explicit or "",
                "relative_orientation": c.relative_orientation,
                "raw_strand_field": c.raw_strand_field,
                "raw_line_numbers": ",".join(map(str, chain.raw_line_numbers)),
            }
            for chain in chains
            for c in chain.components
        ]
        fields = (
            "synteny_chain_id",
            "component_id",
            "query_chrom",
            "query_start",
            "query_end",
            "target_chrom",
            "target_start",
            "target_end",
            "query_strand",
            "target_strand_explicit",
            "relative_orientation",
            "raw_strand_field",
            "raw_line_numbers",
        )
        provenance = {
            "identity": self.identity(),
            "command": command,
            "child_node": child_node,
            "parent_node": parent_node,
            "query_genome": child_genome,
            "target_genome": parent_genome,
            "raw_psl_sha256": hashlib.sha256(raw.encode()).hexdigest(),
            "warnings": [WARNING],
        }
        artifacts = (
            ("raw_halsynteny.psl", raw),
            ("halsynteny_stdout.txt", stdout),
            ("halsynteny_stderr.txt", stderr),
            ("synteny_chains.tsv", tsv(chain_rows, fields)),
            (
                "uncovered_edge_evidence.tsv",
                tsv(
                    uncovered,
                    ("child_block_id", "child_chrom", "child_start", "child_end", "reason"),
                ),
            ),
            (
                "overlay_block_states.tsv",
                tsv([asdict(s) for s in states], tuple(EdgeBlockState.__dataclass_fields__)),
            ),
            (
                "overlay_block_observations.tsv",
                tsv(
                    observations,
                    (
                        "child_block_id",
                        "chain_count",
                        "component_run_count",
                        "overlapping_alternatives",
                        "target_chromosome_count",
                        "relative_orientation_count",
                        "covered_source_bases",
                        "uncovered_source_bases",
                    ),
                ),
            ),
            ("overlay_provenance.json", json.dumps(provenance, indent=2, sort_keys=True) + "\n"),
            ("overlay_statistics.json", json.dumps(metrics, indent=2, sort_keys=True) + "\n"),
        )
        return EdgeExtractionResult(
            runs,
            unmapped,
            BACKEND,
            {"halSynteny": self.runner.executable_version, "halStats": self.stats_version},
            (command,) if command else (),
            (WARNING,),
            artifacts,
            states,
            metrics,
        )
