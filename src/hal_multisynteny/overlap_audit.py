"""Reproducible multiplicity audit of preserved direct-edge PSL evidence."""

from __future__ import annotations

import gzip
import hashlib
import json
import tempfile
from collections import Counter
from itertools import pairwise
from pathlib import Path

from .audit import _read_blocks, _write_tsv
from .source_diagnostic import PSL_PARSER_VERSION, classify_source_mapping, read_source_psl
from .utils import sha256_file

OVERLAP_AUDIT_POLICY = "1"


def geometry(c):
    return (c.q_name, c.q_start, c.q_end, c.t_name, c.t_start, c.t_end, c.relative_orientation)


def analyze_block(block, records, inherited_status=None):
    """Count distinct projections per source atom, independently of PSL rows."""
    if classify_source_mapping(block, records)["classification"] == "invalid_or_ambiguous":
        raise ValueError(f"invalid PSL identity/bounds for {block.block_id}")
    raw = [c for record in records for c in record.components]
    components = {geometry(c): c for c in raw}
    ordered = [components[key] for key in sorted(components)]
    boundaries = sorted(
        {block.start, block.end} | {p for c in ordered for p in (c.q_start, c.q_end)}
    )
    atoms = Counter()
    bases = Counter()
    maximum = 0
    for start, end in pairwise(boundaries):
        projections = set()
        for c in ordered:
            if c.q_start <= start and end <= c.q_end:
                target_start = (
                    c.t_start + start - c.q_start
                    if c.relative_orientation == "+"
                    else c.t_end - (end - c.q_start)
                )
                projections.add(
                    (c.t_name, target_start, target_start + end - start, c.relative_orientation)
                )
        multiplicity = len(projections)
        category = "zero" if multiplicity == 0 else "one" if multiplicity == 1 else "multiple"
        atoms[category] += 1
        bases[category] += end - start
        maximum = max(maximum, multiplicity)
    target_overlaps = sum(
        a.t_name == b.t_name
        and a.t_start < b.t_end
        and b.t_start < a.t_end
        and not (a.q_start < b.q_end and b.q_start < a.q_end)
        for i, a in enumerate(ordered)
        for b in ordered[i + 1 :]
    )
    record_geometries = [
        tuple(sorted(geometry(c) for c in record.components)) for record in records
    ]
    new_ambiguity = maximum > 1 or target_overlaps > 0
    row = {
        "block_id": block.block_id,
        "source_chrom": block.chrom,
        "source_start": block.start,
        "source_end": block.end,
        "input_source_bases": block.end - block.start,
        "raw_psl_record_count": len(records),
        "raw_component_count": len(raw),
        "normalized_component_count": len(ordered),
        "covered_source_bases": bases["one"] + bases["multiple"],
        "uncovered_source_bases": bases["zero"],
        **{f"source_atoms_{name}": atoms[name] for name in ("zero", "one", "multiple")},
        **{f"source_bases_{name}": bases[name] for name in ("zero", "one", "multiple")},
        "maximum_source_projection_multiplicity": maximum,
        "target_chromosome_count": len({c.t_name for c in ordered}),
        "relative_orientation_count": len({c.relative_orientation for c in ordered}),
        "relative_orientations": ",".join(sorted({c.relative_orientation for c in ordered})),
        "chromosome_boundary_count": sum(a.t_name != b.t_name for a, b in pairwise(ordered)),
        "orientation_boundary_count": sum(
            a.relative_orientation != b.relative_orientation for a, b in pairwise(ordered)
        ),
        "non_monotone_order_boundary_count": sum(
            a.t_name == b.t_name
            and a.relative_orientation == b.relative_orientation
            and (a.t_start > b.t_start if a.relative_orientation == "+" else a.t_start < b.t_start)
            for a, b in pairwise(ordered)
        ),
        "source_disjoint_target_overlap_pair_count": target_overlaps,
        "exact_duplicate_psl_row_count": len(records) - len({r.raw for r in records}),
        "duplicate_record_geometry_count": len(records) - len(set(record_geometries)),
        "duplicate_component_geometry_count": len(raw) - len(ordered),
        "inherited_ambiguity": "unavailable"
        if inherited_status is None
        else inherited_status == "ambiguous",
        "inherited_status": inherited_status if inherited_status is not None else "unavailable",
        "newly_observed_ambiguity": new_ambiguity,
        "new_source_multiplicity": maximum > 1,
        "atomwise_unique": maximum == 1,
    }
    if row["covered_source_bases"] + row["uncovered_source_bases"] != row["input_source_bases"]:
        raise ValueError("overlap audit source-base conservation failed")
    return row


def _totals(rows, child, parent):
    numeric = [
        key
        for key, value in rows[0].items()
        if type(value) is int
        and key not in {"source_start", "source_end", "maximum_source_projection_multiplicity"}
    ]
    problematic = sum(row["newly_observed_ambiguity"] for row in rows)
    return {
        "child_genome": child,
        "parent_genome": parent,
        "input_block_count": len(rows),
        **{key: sum(row[key] for row in rows) for key in numeric},
        "maximum_source_projection_multiplicity": max(
            row["maximum_source_projection_multiplicity"] for row in rows
        ),
        "blocks_with_source_multiplicity": sum(row["new_source_multiplicity"] for row in rows),
        "blocks_with_new_ambiguity": problematic,
        "new_ambiguity_fraction": problematic / len(rows),
        "inherited_status_unavailable_blocks": sum(
            row["inherited_status"] == "unavailable" for row in rows
        ),
    }


def audit_preserved_psl(audit_dir, output_prefix):
    """Discover every preserved source-aware capture; fail on changed raw bytes.

    No subprocesses, HAL files, or biological reconstruction are used. Sampling
    inputs have no upward block status: inherited ambiguity is unavailable, never
    silently inferred absent from an ancestral genome label.
    """
    audit_dir, prefix = Path(audit_dir), Path(output_prefix)
    evidence_dir = audit_dir.parent / "evidence"
    raw_paths = sorted(evidence_dir.glob("*.raw.psl.gz"))
    portable = bool(raw_paths and (evidence_dir / "manifest.json").exists())
    if not raw_paths:
        raw_paths = sorted(audit_dir.glob("results*/source-aware/*.raw.psl"))
    if not raw_paths:
        raise ValueError("no preserved source-aware PSL captures found")
    all_rows, edge_rows, evidence = [], [], []
    seen = set()
    for raw in raw_paths:
        if portable:
            manifest = json.loads((evidence_dir / "manifest.json").read_text())
            entry = next((row for row in manifest["captures"] if row["compressed_file"] == raw.name), None)
            if entry is None:
                raise ValueError(f"evidence manifest has no entry for {raw.name}")
            compressed_checksum = sha256_file(raw)
            with gzip.open(raw, "rb") as handle:
                raw_bytes = handle.read()
            checksum = hashlib.sha256(raw_bytes).hexdigest()
            if compressed_checksum != entry["compressed_sha256"] or checksum != entry["raw_sha256"]:
                raise ValueError(f"compressed or raw PSL checksum changed: {raw}")
            child, parent = entry["child_genome"], entry["parent_genome"]
            sampled = evidence_dir
            paths = [sampled / entry["sampled_file"]]
            provenance = {"hal_sha256": entry.get("hal_sha256", manifest["hal_sha256"])}
            with tempfile.NamedTemporaryFile(suffix=".psl") as handle:
                handle.write(raw_bytes)
                handle.flush()
                grouped = read_source_psl(handle.name)
        else:
            provenance = json.loads(raw.with_name(raw.name.replace(".raw.psl", ".provenance.json")).read_text())
            checksum = sha256_file(raw)
            if checksum != provenance["raw_psl_sha256"]:
                raise ValueError(f"raw PSL checksum changed: {raw}")
            child, parent = provenance["child_genome"], provenance["parent_genome"]
            sampled = raw.parent.parent / "sampled"
            paths = sorted(sampled.glob(f"*__{child}__to__{parent}.blocks.tsv"))
            grouped = read_source_psl(raw)
        if len(paths) != 1:
            raise ValueError(f"expected one sampled input for {child}->{parent}")
        blocks_path = paths[0]
        if portable:
            if sha256_file(blocks_path) != entry["sampled_sha256"]:
                raise ValueError(f"sampled block checksum changed: {blocks_path}")
        else:
            sample_manifest = json.loads((sampled / "sampling-manifest.json").read_text())
            if sample_manifest["identity"]["hal_sha256"] != provenance["hal_sha256"]:
                raise ValueError("sampling and PSL HAL identities disagree")
            if not any(Path(value).name == blocks_path.name for value in sample_manifest["outputs"].values()):
                raise ValueError("sampled file absent from original manifest")
        capture = (child, parent, checksum, sha256_file(blocks_path))
        if capture in seen:
            continue
        seen.add(capture)
        blocks = _read_blocks(blocks_path)
        known = {b.block_id for b in blocks}
        if len(known) != len(blocks) or set(grouped) - known:
            raise ValueError("duplicate sampled block IDs or unknown PSL names")
        rows = [
            {
                "child_genome": child,
                "parent_genome": parent,
                **analyze_block(block, grouped.get(block.block_id, [])),
            }
            for block in blocks
        ]
        rows.sort(key=lambda row: row["block_id"])
        all_rows.extend(rows)
        edge_rows.append(_totals(rows, child, parent))
        evidence.append(
            {
                "child_genome": child,
                "parent_genome": parent,
                "raw_psl": str(raw.relative_to(evidence_dir if portable else audit_dir)),
                "raw_psl_sha256": checksum,
                "sampled_blocks": blocks_path.name,
                "sampled_blocks_sha256": sha256_file(blocks_path),
                "hal_sha256": provenance["hal_sha256"],
            }
        )
    planned_edges = set()
    for plan_path in sorted(audit_dir.glob("results*/plan.json")):
        for edge in json.loads(plan_path.read_text()).get("edges", []):
            planned_edges.add((edge["child_hal_genome"], edge["parent_hal_genome"]))
    covered_edges = {(edge["child_genome"], edge["parent_genome"]) for edge in edge_rows}
    missing_edges = [
        {"child_genome": child, "parent_genome": parent}
        for child, parent in sorted(planned_edges - covered_edges)
    ]
    multiple = [row for row in all_rows if row["new_source_multiplicity"]]
    problematic = []
    for original in all_rows:
        if not (original["new_source_multiplicity"] or original["source_disjoint_target_overlap_pair_count"]):
            continue
        row = dict(original)
        row["problem_type"] = (
            "both" if row["new_source_multiplicity"] and row["source_disjoint_target_overlap_pair_count"]
            else "source_multiplicity" if row["new_source_multiplicity"]
            else "same_target_overlap"
        )
        problematic.append(row)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix, rows, fields in (
        ("blocks", all_rows, list(all_rows[0])),
        ("edges", edge_rows, list(edge_rows[0])),
        ("multiplicity", multiple, list(all_rows[0])),
        ("problematic", problematic, list(problematic[0]) if problematic else list(all_rows[0])),
    ):
        _write_tsv(Path(f"{prefix}.{suffix}.tsv"), fields, rows)
    payload = {
        "policy_version": OVERLAP_AUDIT_POLICY,
        "parser_version": PSL_PARSER_VERSION,
        "method": "reparse existing direct-edge PSL; no HAL invocation",
        "inherited_ambiguity": "unavailable: sampled intervals lack upward block/occurrence states",
        "definitions": {
            "normalized_component": "distinct source/target geometry including chromosome and relative orientation",
            "boundaries": "adjacent components sorted by source start/end then target geometry; alternatives included",
            "atomwise_unique": "one distinct projection per covered source atom; target overlap may remain",
            "duplicate_counts": "excess rows/components beyond first identical member",
            "source_atoms": "disjoint after atomization; target projections may overlap or be one-to-many",
        },
        "input_block_count": len(all_rows),
        "edges": edge_rows,
        "evidence": evidence,
        "planned_edges_without_psl": missing_edges,
        "blocks_with_source_multiplicity": multiple,
        "problematic_block_count": len(problematic),
        "problematic_component_overlap_pair_count": sum(row["source_disjoint_target_overlap_pair_count"] for row in problematic),
    }
    Path(f"{prefix}.json").write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n")
    lines = [
        "# Preserved Palaeognathae PSL multiplicity audit",
        "",
        f"Audited {len(all_rows):,} input intervals across {len(edge_rows)} preserved direct-edge captures.",
        "No HAL commands were executed. Coverage is limited to edges with existing PSL evidence; BED6-only edges are excluded.",
        "",
        "Source atoms are disjoint after atomization; target projections may overlap or be one-to-many. Atomwise unique does not establish biological uniqueness or synteny.",
        "",
        "| Child → parent | Inputs | Problematic blocks | Source multiplicity >1 | Same-target overlap pairs | Max multiplicity |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for edge in edge_rows:
        lines.append(
            f"| {edge['child_genome']} → {edge['parent_genome']} | {edge['input_block_count']} | {edge['blocks_with_new_ambiguity']} | {edge['blocks_with_source_multiplicity']} | {edge['source_disjoint_target_overlap_pair_count']} | {edge['maximum_source_projection_multiplicity']} |"
        )
    total_problematic = sum(edge["blocks_with_new_ambiguity"] for edge in edge_rows)
    total_overlap_pairs = sum(edge["source_disjoint_target_overlap_pair_count"] for edge in edge_rows)
    lines += [
        "",
        f"Current-edge conflict affects {total_problematic:,} blocks and contains {total_overlap_pairs:,} overlapping component pairs. Per-edge concentration (share of all problematic blocks; within-edge rate):",
    ]
    for edge in sorted(
        edge_rows, key=lambda e: (-e["blocks_with_new_ambiguity"], e["child_genome"])
    ):
        share = edge["blocks_with_new_ambiguity"] / total_problematic if total_problematic else 0
        lines.append(
            f"- {edge['child_genome']} → {edge['parent_genome']}: {share:.2%} of problematic inputs; {edge['new_ambiguity_fraction']:.2%} within edge."
        )
    if missing_edges:
        lines += ["", "Planned direct edges without preserved PSL (excluded from this audit):"]
        lines += [f"- {edge['child_genome']} → {edge['parent_genome']}" for edge in missing_edges]
    lines += [
        "",
        f"Every problematic input ({len(problematic):,} blocks; {sum(row['source_disjoint_target_overlap_pair_count'] for row in problematic):,} overlapping component pairs) is listed in `{prefix.name}.problematic.tsv`; source-multiplicity inputs are listed separately in `{prefix.name}.multiplicity.tsv`.",
        "Inherited ambiguity is unavailable for these randomly sampled intervals because no upward state was captured. Newly observed source multiplicity and same-target overlap are reported separately; ancestor names do not imply inherited ambiguity.",
        "Raw PSL checksums are verified against capture provenance. Sampling-file checksums are recorded now; the historical sampling manifest did not checksum its outputs.",
        "Statistics describe sampled mappings, not full-genome rates, biological copy identities, or liftover transitivity.",
        "",
    ]
    Path(f"{prefix}.md").write_text("\n".join(lines))
    return payload
