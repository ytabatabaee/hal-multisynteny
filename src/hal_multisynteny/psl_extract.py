"""Conservative, source-aware direct-edge HAL extractor; no multi-record resolver."""
from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
from collections import Counter
from collections.abc import Iterable
from pathlib import Path

from .audit import AuditBlock
from .extract import EdgeExtractionResult, ExtractionError, HalEdgeExtractor
from .models import EdgeMappingRun, NodeBlock, UnmappedEdgeEvidence
from .source_diagnostic import PSL_PARSER_VERSION, classify_source_mapping, read_source_psl
from .utils import sha256_file

BACKEND = "hal-psl-conservative"
PSL_EDGE_SCHEMA_VERSION = "1"
WARNING = (
    "Multiple PSL records are unresolved hypotheses, not biological copies; "
    "their components are ambiguity barriers and are not propagated as biological occurrences. "
    "halSynteny transitivity and general mapping-path resolution remain unestablished."
)


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
        runs, unmapped = [], []
        for block in blocks:
            records = grouped.get(block.block_id, [])
            diagnostic = classify_source_mapping(AuditBlock(block.block_id, block.chrom, block.start, block.end), records)
            if diagnostic["classification"] == "invalid_or_ambiguous":
                raise ExtractionError(f"invalid PSL source identity/bounds for {block.block_id}")
            if not records:
                unmapped.append(UnmappedEdgeEvidence(child_node_id, parent_node_id, block.block_id,
                    block.chrom, block.start, block.end, "unaligned", "halLiftover PSL returned zero records",
                    f"{child_node_id}:{block.block_id}:psl:unaligned", "halLiftover", json.dumps(command)))
                continue
            # Prior unresolved evidence must remain a barrier even when this edge
            # has a single alignment record. It cannot be rehabilitated by liftover.
            ambiguous = len(records) > 1 or block.mapping_status == "ambiguous" or block.classification in {"ambiguous", "complex"}
            status = "ambiguous" if ambiguous else "unique"
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
                        "" if ambiguous else "1", status,
                        f"{occurrence_id}:{hypothesis}:component{c.component_order:06d}", "halLiftover", json.dumps(command),
                        hypothesis, c.query_strand, c.target_strand_explicit or "", c.raw_strand_field,
                        c.psl_line_number, c.component_order))
        runs.sort(key=lambda r: (r.parent_chrom, r.parent_start, r.parent_end, r.child_block_id,
                                 r.hypothesis_id, r.component_order))
        provenance = {**self.identity(), "commands": [command], "warnings": [WARNING],
                      "raw_psl_sha256": hashlib.sha256(raw.encode()).hexdigest(),
                      "outcomes": {b.block_id: ("unaligned" if not grouped.get(b.block_id) else
                                   "unique" if len(grouped[b.block_id]) == 1 else "ambiguous") for b in blocks}}
        return EdgeExtractionResult(tuple(runs), tuple(unmapped), BACKEND, self.tool_versions(), (command,), (WARNING,),
            (("raw.psl", raw), ("stdout.txt", completed.stdout), ("stderr.txt", completed.stderr),
             ("psl_provenance.json", json.dumps(provenance, sort_keys=True, indent=2) + "\n")))
