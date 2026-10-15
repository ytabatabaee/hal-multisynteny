"""Conservative, source-aware direct-edge HAL extractor; no multi-record resolver."""
from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
from collections import Counter
from collections.abc import Iterable
from dataclasses import asdict
from itertools import pairwise
from pathlib import Path

from .audit import AuditBlock
from .extract import EdgeExtractionResult, ExtractionError, HalEdgeExtractor
from .models import EdgeBlockState, EdgeMappingRun, NodeBlock, UnmappedEdgeEvidence
from .source_diagnostic import PSL_PARSER_VERSION, classify_source_mapping, read_source_psl
from .utils import sha256_file

BACKEND = "hal-psl-conservative"
PSL_EDGE_SCHEMA_VERSION = "2"
STATUS_POLICY_VERSION = "1"
WARNING = (
    "Multiple PSL records are unresolved hypotheses, not biological copies; "
    "their components are ambiguity barriers and are not propagated as biological occurrences. "
    "halSynteny transitivity and general mapping-path resolution remain unestablished."
)

# Orientation/order conflicts alone are not inherited barriers. A compound
# conflict is already classified as complex and conservatively remains ambiguous.
CLASSIFICATION_BARRIERS = {"ambiguous": "ambiguous", "complex": "ambiguous",
                           "duplication_conflict": "duplicated"}
MAPPING_STATUS_BARRIERS = {"ambiguous": "ambiguous", "duplicated": "duplicated"}


def edge_block_state(
    block: NodeBlock, record_count: int, *, raw_outcome: str | None = None
) -> EdgeBlockState:
    barriers = [
        (MAPPING_STATUS_BARRIERS.get(block.mapping_status), f"inherited mapping_status={block.mapping_status}"),
        (CLASSIFICATION_BARRIERS.get(block.classification), f"inherited classification={block.classification}"),
    ]
    inherited = ("ambiguous" if any(status == "ambiguous" for status, _ in barriers) else
                 "duplicated" if any(status == "duplicated" for status, _ in barriers) else "")
    raw = raw_outcome or ("unaligned" if record_count == 0 else "unique" if record_count == 1 else "ambiguous")
    effective = "unaligned" if record_count == 0 else "ambiguous" if raw == "ambiguous" else inherited or "unique"
    return EdgeBlockState(block.block_id, raw, effective, inherited,
                          "; ".join(reason for status, reason in barriers if status),
                          block.mapping_status, block.classification)


class HalPslConservativeExtractor(HalEdgeExtractor):
    """One logical input occurrence per block; one hypothesis per output record.

    Spatial ambiguous runs enter the sweep as barriers only. All records are
    retained (including identical rows); no PSL-record count becomes a copy count.
    """

    def __init__(self, hal_path: str | Path, **kwargs) -> None:
        super().__init__(hal_path, **kwargs)
        from .hal import _parse_parent_map

        self.hal_path = self.hal_path.resolve()
        self.hal_liftover = str(Path(self.hal_liftover).resolve())
        self.hal_stats = str(Path(self.hal_stats).resolve())
        tree = subprocess.run((self.hal_stats, str(self.hal_path), "--tree"),
                              check=False, capture_output=True, text=True)
        if tree.returncode:
            raise ExtractionError(f"halStats --tree failed: {tree.stderr.strip()}")
        self.parent_map = _parse_parent_map(tree.stdout)
        if not self.parent_map:
            raise ExtractionError("HAL tree has no validated direct edges")
        self._tools = {}
        for name, tool in (("halStats", self.hal_stats), ("halLiftover", self.hal_liftover)):
            help_run = subprocess.run((tool, "--help"), check=False, capture_output=True, text=True)
            lines = [line.strip() for line in (help_run.stdout + help_run.stderr).splitlines() if line.strip()]
            version = next((s for s in lines if " v" in s), lines[0] if lines else "unknown")
            if f"{name} v" in version:
                version = version[version.index(f"{name} v"):]
            self._tools[name] = {"path": tool, "version": version,
                                 "sha256": sha256_file(tool)}

    def identity(self) -> dict[str, object]:
        return {"backend": BACKEND, "parser_version": PSL_PARSER_VERSION,
                "edge_schema_version": PSL_EDGE_SCHEMA_VERSION, "hal_sha256": self.hal_checksum,
                "status_policy_version": STATUS_POLICY_VERSION,
                "parent_map": self.parent_map, "tools": self._tools}

    def tool_versions(self) -> dict[str, str]:
        return {name: value["version"] for name, value in self._tools.items()}

    def extract(self, child_node_id: str, parent_node_id: str, child_hal_genome: str,
                parent_hal_genome: str, child_blocks: Iterable[NodeBlock]) -> EdgeExtractionResult:
        if self.parent_map.get(child_hal_genome) != parent_hal_genome:
            raise ExtractionError(f"not a direct HAL child-to-parent edge: {child_hal_genome}->{parent_hal_genome}")
        blocks = tuple(sorted(child_blocks, key=lambda b: (b.chrom, b.start, b.end, b.block_id)))
        if len({b.block_id for b in blocks}) != len(blocks) or any(b.node != child_node_id for b in blocks):
            raise ExtractionError("duplicate block IDs or blocks from a different child node")
        if not blocks:
            return self._result((), (), (), "", "", "", ())
        command = (self.hal_liftover, "--outPSLWithName", str(self.hal_path), child_hal_genome,
                   "child.bed", parent_hal_genome, "parent.psl")
        # Relative file names and an isolated cwd preserve the exact invoked
        # command without leaking random temporary directory names into identity.
        with tempfile.TemporaryDirectory(prefix="hms-psl-") as temp_name:
            temp = Path(temp_name)
            (temp / "child.bed").write_text("".join(
                f"{b.chrom}\t{b.start}\t{b.end}\t{b.block_id}\t0\t+\n" for b in blocks), encoding="utf-8")
            completed = subprocess.run(command, cwd=temp, check=False, capture_output=True, text=True)
            if completed.returncode:
                raise ExtractionError(f"halLiftover PSL failed: {completed.stderr.strip() or completed.stdout.strip()}")
            path = temp / "parent.psl"
            if not path.is_file():
                raise ExtractionError("halLiftover did not create parent.psl")
            raw = path.read_text(encoding="utf-8")
            grouped = read_source_psl(path)
        if set(grouped) - {b.block_id for b in blocks}:
            raise ExtractionError("halLiftover PSL returned unknown block IDs")
        runs, unmapped, states = [], [], []
        for block in blocks:
            records = grouped.get(block.block_id, [])
            diagnostic = classify_source_mapping(AuditBlock(block.block_id, block.chrom, block.start, block.end), records)
            if diagnostic["classification"] == "invalid_or_ambiguous":
                raise ExtractionError(f"invalid PSL source identity/bounds for {block.block_id}")
            state = edge_block_state(block, len(records))
            states.append(state)
            if not records:
                unmapped.append(UnmappedEdgeEvidence(child_node_id, parent_node_id, block.block_id,
                    block.chrom, block.start, block.end, "unaligned", "halLiftover PSL returned zero records",
                    f"{child_node_id}:{block.block_id}:psl:unaligned", "halLiftover", json.dumps(command),
                    state.raw_mapping_outcome, state.inherited_barrier_status, state.inherited_barrier_reason))
                continue
            # Prior unresolved evidence must remain a barrier even when this edge
            # has a single alignment record. It cannot be rehabilitated by liftover.
            status = state.effective_status
            occurrence_id = f"{child_node_id}:{block.block_id}:psl-input"
            keyed = []
            for record in records:
                key = tuple((c.q_name, c.q_size, c.q_start, c.q_end, c.t_name, c.t_start, c.t_end,
                             c.query_strand, c.target_strand_explicit or "", c.relative_orientation)
                            for c in record.components)
                keyed.append((json.dumps(key, separators=(",", ":")), record))
            duplicates = Counter()
            for key, record in sorted(keyed, key=lambda item: (item[0], item[1].raw)):
                digest = hashlib.sha256(key.encode()).hexdigest()[:24]
                duplicates[digest] += 1
                hypothesis = f"PSLH_{digest}_{duplicates[digest]:04d}"
                for c in record.components:
                    runs.append(EdgeMappingRun(child_node_id, block.block_id, occurrence_id, c.q_name, c.q_start, c.q_end,
                        parent_node_id, c.t_name, c.t_start, c.t_end, c.relative_orientation,
                        "1" if status == "unique" else "", status,
                        f"{occurrence_id}:{hypothesis}:component{c.component_order:06d}", "halLiftover", json.dumps(command),
                        hypothesis, c.query_strand, c.target_strand_explicit or "", c.raw_strand_field,
                        c.psl_line_number, c.component_order, state.raw_mapping_outcome,
                        state.inherited_barrier_status, state.inherited_barrier_reason))
        runs.sort(key=lambda r: (r.parent_chrom, r.parent_start, r.parent_end, r.child_block_id,
                                 r.hypothesis_id, r.component_order))
        return self._result(tuple(runs), tuple(unmapped), (command,), raw, completed.stdout, completed.stderr, tuple(states))

    def _result(self, runs, unmapped, commands, raw, stdout, stderr, states):
        provenance = {**self.identity(), "commands": commands, "warnings": [WARNING],
                      "raw_psl_sha256": hashlib.sha256(raw.encode()).hexdigest(),
                      "outcomes": {state.child_block_id: state.raw_mapping_outcome for state in states},
                      "block_states": [asdict(state) for state in states]}
        return EdgeExtractionResult(runs, unmapped, BACKEND, self.tool_versions(), commands, (WARNING,),
            (("raw.psl", raw), ("stdout.txt", stdout), ("stderr.txt", stderr),
             ("psl_provenance.json", json.dumps(provenance, sort_keys=True, indent=2) + "\n")), states)


PIECEWISE_BACKEND = "hal-psl-piecewise"
PIECEWISE_SCHEMA_VERSION = "1"
PIECEWISE_STATUS_POLICY_VERSION = "1"


def _component_key(component) -> tuple[object, ...]:
    return (
        component.q_name, component.q_start, component.q_end,
        component.t_name, component.t_start, component.t_end,
        component.relative_orientation, component.query_strand,
        component.target_strand_explicit or "",
    )


def _piece_projection(component, start: int, end: int) -> tuple[int, int]:
    if component.relative_orientation == "+":
        parent_start = component.t_start + (start - component.q_start)
        parent_end = component.t_start + (end - component.q_start)
    else:
        parent_start = component.t_end - (end - component.q_start)
        parent_end = component.t_end - (start - component.q_start)
    if end - start != parent_end - parent_start:
        raise ExtractionError("piecewise PSL projection is not length preserving")
    return parent_start, parent_end


def _interval_overlap(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return left[0] < right[1] and right[0] < left[1]


class HalPslPiecewiseExtractor(HalPslConservativeExtractor):
    """Resolve source-disjoint liftover records into exact mapping fragments.

    PSL records are output records, not biological copies.  This backend
    partitions each input block by normalized component boundaries and calls a
    block unique only when every covered source atom has one distinct target
    projection.  Overlapping projections remain ambiguous; source-disjoint
    pieces remain spatial evidence with one logical child occurrence.
    """

    def identity(self) -> dict[str, object]:
        identity = super().identity()
        identity.update({
            "backend": PIECEWISE_BACKEND,
            "edge_schema_version": PIECEWISE_SCHEMA_VERSION,
            "status_policy_version": PIECEWISE_STATUS_POLICY_VERSION,
        })
        return identity

    def extract(self, child_node_id: str, parent_node_id: str, child_hal_genome: str,
                parent_hal_genome: str, child_blocks: Iterable[NodeBlock]) -> EdgeExtractionResult:
        if self.parent_map.get(child_hal_genome) != parent_hal_genome:
            raise ExtractionError(f"not a direct HAL child-to-parent edge: {child_hal_genome}->{parent_hal_genome}")
        blocks = tuple(sorted(child_blocks, key=lambda b: (b.chrom, b.start, b.end, b.block_id)))
        if len({b.block_id for b in blocks}) != len(blocks) or any(b.node != child_node_id for b in blocks):
            raise ExtractionError("duplicate block IDs or blocks from a different child node")
        if not blocks:
            return self._piecewise_result((), (), (), "", "", "", ())

        command = (self.hal_liftover, "--outPSLWithName", str(self.hal_path), child_hal_genome,
                   "child.bed", parent_hal_genome, "parent.psl")
        with tempfile.TemporaryDirectory(prefix="hms-psl-piecewise-") as temp_name:
            temp = Path(temp_name)
            (temp / "child.bed").write_text("".join(
                f"{b.chrom}\t{b.start}\t{b.end}\t{b.block_id}\t0\t+\n" for b in blocks), encoding="utf-8")
            completed = subprocess.run(command, cwd=temp, check=False, capture_output=True, text=True)
            if completed.returncode:
                raise ExtractionError(f"halLiftover PSL failed: {completed.stderr.strip() or completed.stdout.strip()}")
            path = temp / "parent.psl"
            if not path.is_file():
                raise ExtractionError("halLiftover did not create parent.psl")
            raw = path.read_text(encoding="utf-8")
            grouped = read_source_psl(path)
        if set(grouped) - {b.block_id for b in blocks}:
            raise ExtractionError("halLiftover PSL returned unknown block IDs")

        runs: list[EdgeMappingRun] = []
        unmapped: list[UnmappedEdgeEvidence] = []
        states: list[EdgeBlockState] = []
        block_statistics: list[dict[str, object]] = []
        for block in blocks:
            records = grouped.get(block.block_id, [])
            diagnostic = classify_source_mapping(AuditBlock(block.block_id, block.chrom, block.start, block.end), records)
            if diagnostic["classification"] == "invalid_or_ambiguous":
                raise ExtractionError(f"invalid PSL source identity/bounds for {block.block_id}")
            component_by_key = {}
            for component in (component for record in records for component in record.components):
                component_by_key.setdefault(_component_key(component), component)
            components = [component_by_key[key] for key in sorted(component_by_key)]
            components = [component for component in components
                          if component.q_start < block.end and block.start < component.q_end]
            boundaries = sorted({block.start, block.end} | {
                point for component in components
                for point in (max(block.start, component.q_start), min(block.end, component.q_end))
            })
            atoms: list[tuple[int, int, tuple[object, ...]]] = []
            unique_atoms = 0
            overlapping_atoms = 0
            for start, end in pairwise(boundaries):
                projections = []
                for component in components:
                    if component.q_start <= start and end <= component.q_end:
                        parent_start, parent_end = _piece_projection(component, start, end)
                        projections.append((component.t_name, parent_start, parent_end, component.relative_orientation))
                distinct = tuple(sorted(set(projections)))
                if distinct:
                    if len(distinct) == 1:
                        unique_atoms += 1
                    else:
                        overlapping_atoms += 1
                atoms.append((start, end, distinct))
            mapped_bases = sum(end - start for start, end, projections in atoms if projections)
            uncovered_bases = (block.end - block.start) - mapped_bases
            source_disjoint = all(not _interval_overlap((a.q_start, a.q_end), (b.q_start, b.q_end))
                                  for index, a in enumerate(components) for b in components[index + 1:])
            target_overlap_conflict = any(
                _interval_overlap((a.t_start, a.t_end), (b.t_start, b.t_end))
                and not _interval_overlap((a.q_start, a.q_end), (b.q_start, b.q_end))
                for index, a in enumerate(components) for b in components[index + 1:]
            )
            ordered_components = sorted(components, key=lambda component: (component.q_start, component.q_end, _component_key(component)))
            chromosome_boundary_count = sum(
                left.t_name != right.t_name
                for left, right in pairwise(ordered_components)
            )
            orientation_boundary_count = sum(
                left.relative_orientation != right.relative_orientation
                for left, right in pairwise(ordered_components)
            )
            order_boundary_count = sum(
                left.t_name == right.t_name and
                left.relative_orientation == right.relative_orientation and
                ((left.t_start > right.t_start) if left.relative_orientation == "+" else (left.t_start < right.t_start))
                for left, right in pairwise(ordered_components)
            )
            raw_outcome = "unaligned" if not components else "ambiguous" if overlapping_atoms else "unique"
            resolution = (
                "unaligned" if raw_outcome == "unaligned" else
                "overlapping_alternatives" if overlapping_atoms else
                "target_overlap_conflict" if target_overlap_conflict else
                "source_disjoint_fragments" if len(records) > 1 or len(components) > 1 else
                "single_mapping"
            )
            state = edge_block_state(block, len(records), raw_outcome=raw_outcome)
            if target_overlap_conflict and state.effective_status == "unique":
                state = EdgeBlockState(state.child_block_id, state.raw_mapping_outcome, "duplicated",
                    "duplicated", "source-disjoint components overlap in target coordinates",
                    state.inherited_mapping_status, state.inherited_classification)
            states.append(state)
            line_numbers = sorted({record.psl_line_number for record in records})
            block_statistics.append({
                "block_id": block.block_id, "input_block_bases": block.end - block.start,
                "raw_psl_record_count": len(records), "normalized_distinct_mapping_count": len(components),
                "component_count": len(components), "covered_source_bases": mapped_bases,
                "uncovered_source_bases": uncovered_bases, "source_atom_count": max(0, len(boundaries) - 1),
                "uniquely_mapped_atom_count": unique_atoms, "overlapping_alternative_atom_count": overlapping_atoms,
                "records_source_disjoint": source_disjoint, "target_chromosome_count": len({c.t_name for c in components}),
                "relative_orientation_count": len({c.relative_orientation for c in components}),
                "target_overlap_conflict": target_overlap_conflict,
                "chromosome_boundary_count": chromosome_boundary_count,
                "orientation_boundary_count": orientation_boundary_count,
                "order_boundary_count": order_boundary_count,
                "raw_mapping_outcome": state.raw_mapping_outcome, "effective_status": state.effective_status,
                "inherited_barrier_status": state.inherited_barrier_status,
                "inherited_barrier_reason": state.inherited_barrier_reason,
                "resolution_category": "inherited_barrier" if state.inherited_barrier_status else resolution,
                "psl_line_numbers": line_numbers,
            })
            if not components:
                unmapped.append(UnmappedEdgeEvidence(child_node_id, parent_node_id, block.block_id,
                    block.chrom, block.start, block.end, "unaligned", "halLiftover PSL returned no mapped component",
                    f"{child_node_id}:{block.block_id}:piecewise:unaligned", "halLiftover", json.dumps(command),
                    state.raw_mapping_outcome, state.inherited_barrier_status, state.inherited_barrier_reason))
                continue
            occurrence_id = f"{child_node_id}:{block.block_id}:piecewise-input"
            effective_status = state.effective_status
            for component_index, component in enumerate(components, start=1):
                digest = hashlib.sha256(json.dumps(_component_key(component), separators=(",", ":")).encode()).hexdigest()[:24]
                source_anchor = f"{occurrence_id}:fragment:{digest}:{component_index:06d}"
                hypothesis = ""
                if effective_status == "ambiguous":
                    hypothesis = f"PSLP_{digest}"
                runs.append(EdgeMappingRun(
                    child_node_id, block.block_id, occurrence_id, component.q_name, component.q_start,
                    component.q_end, parent_node_id, component.t_name, component.t_start, component.t_end,
                    component.relative_orientation, "1" if effective_status == "unique" else "", effective_status,
                    source_anchor, "halLiftover", json.dumps(command), hypothesis, component.query_strand,
                    component.target_strand_explicit or "", component.raw_strand_field, component.psl_line_number,
                    component.component_order, state.raw_mapping_outcome, state.inherited_barrier_status,
                    state.inherited_barrier_reason, "", f"piece_{digest}_{component_index:06d}"))
        raw_counts = {status: sum(state.raw_mapping_outcome == status for state in states)
                      for status in ("unique", "ambiguous", "unaligned")}
        effective_counts = {status: sum(state.effective_status == status for state in states)
                            for status in ("unique", "duplicated", "ambiguous", "unaligned")}
        if sum(raw_counts.values()) != len(blocks) or sum(effective_counts.values()) != len(blocks):
            raise ExtractionError("piecewise mapping status counts do not reconcile with input blocks")
        if sum(int(row["covered_source_bases"]) + int(row["uncovered_source_bases"]) for row in block_statistics) != sum(b.end - b.start for b in blocks):
            raise ExtractionError("piecewise covered and uncovered source bases do not reconcile")
        statistics = {
            "denominator_definition": "one row per input NodeBlock; source bases are interval-union totals",
            "input_block_count": len(blocks), "raw_mapping_outcome_counts": raw_counts,
            "effective_status_counts": effective_counts, "blocks": block_statistics,
            "component_count": sum(int(row["component_count"]) for row in block_statistics),
            "normalized_distinct_mapping_count": sum(int(row["normalized_distinct_mapping_count"]) for row in block_statistics),
            "covered_source_bases": sum(int(row["covered_source_bases"]) for row in block_statistics),
            "uncovered_source_bases": sum(int(row["uncovered_source_bases"]) for row in block_statistics),
            "blocks_with_no_chain_overlap": sum(row["raw_mapping_outcome"] == "unaligned" for row in block_statistics),
            "blocks_with_overlapping_alternative_atoms": sum(int(row["overlapping_alternative_atom_count"]) > 0 for row in block_statistics),
            "target_overlap_conflict_count": sum(bool(row["target_overlap_conflict"]) for row in block_statistics),
            "chromosome_boundary_count": sum(int(row["chromosome_boundary_count"]) for row in block_statistics),
            "orientation_boundary_count": sum(int(row["orientation_boundary_count"]) for row in block_statistics),
            "order_boundary_count": sum(int(row["order_boundary_count"]) for row in block_statistics),
        }
        return self._piecewise_result(tuple(runs), tuple(unmapped), (command,), raw,
                                      completed.stdout, completed.stderr, tuple(states), statistics)

    def _piecewise_result(self, runs, unmapped, commands, raw, stdout, stderr, states, statistics=None):
        provenance = {**self.identity(), "commands": commands, "warnings": [
            "Source-disjoint PSL records are retained as mapping fragments; this does not establish biological synteny or copy correspondence."
        ], "raw_psl_sha256": hashlib.sha256(raw.encode()).hexdigest(),
            "block_states": [asdict(state) for state in states], "statistics": statistics or {}}
        return EdgeExtractionResult(tuple(sorted(runs, key=lambda r: (r.parent_chrom, r.parent_start,
            r.parent_end, r.child_block_id, r.source_anchor_id))), tuple(unmapped), PIECEWISE_BACKEND,
            self.tool_versions(), commands, ("PSL record count is not treated as mapping multiplicity",),
            (("raw.psl", raw), ("stdout.txt", stdout), ("stderr.txt", stderr),
             ("piecewise_statistics.json", json.dumps(statistics or {}, sort_keys=True, indent=2) + "\n"),
             ("psl_provenance.json", json.dumps(provenance, sort_keys=True, indent=2) + "\n")),
            states, statistics or {})
