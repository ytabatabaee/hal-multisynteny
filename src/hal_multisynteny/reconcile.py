"""Deterministic reconciliation of two child block systems at one parent node."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from hashlib import sha256
from itertools import pairwise

from .models import ParentMappedRun

CLASSIFICATIONS = (
    "shared_consistent",
    "left_only",
    "right_only",
    "orientation_conflict",
    "order_conflict",
    "duplication_conflict",
    "ambiguous",
    "unaligned",
    "complex",
)

COPY_ID_SCOPES = ("local", "global")


def _ids(records: Iterable[ParentMappedRun], attribute: str) -> str:
    return ",".join(sorted({getattr(record, attribute) for record in records}))


@dataclass(frozen=True, slots=True)
class ReconcileConfig:
    min_block_length: int = 50
    boundary_tolerance: int = 1
    max_merge_gap: int = 0
    guide_tree_id: str | None = None
    copy_id_scope: str = "local"

    def __post_init__(self) -> None:
        if self.min_block_length < 1:
            raise ValueError("min_block_length must be at least 1")
        if self.boundary_tolerance < 0:
            raise ValueError("boundary_tolerance must be non-negative")
        if self.max_merge_gap < 0:
            raise ValueError("max_merge_gap must be non-negative")
        if self.copy_id_scope not in COPY_ID_SCOPES:
            raise ValueError(f"copy_id_scope must be one of {list(COPY_ID_SCOPES)}")


@dataclass(frozen=True, slots=True)
class AtomicInterval:
    atomic_interval_id: str
    parent_node: str
    parent_chrom: str
    parent_start: int
    parent_end: int
    left_record_ids: str
    right_record_ids: str
    left_block_ids: str
    right_block_ids: str
    classification: str
    left_orientations: str
    right_orientations: str
    left_copies: str
    right_copies: str
    disposition: str
    merge_reason: str
    final_parent_block_id: str


@dataclass(frozen=True, slots=True)
class ParentBlock:
    block_id: str
    parent_node: str
    parent_chrom: str
    parent_start: int
    parent_end: int
    classification: str
    left_block_ids: str
    right_block_ids: str
    species_count: int
    occurrence_count: int
    copy_count: int
    mapping_quality: str
    atomic_interval_count: int


@dataclass(frozen=True, slots=True)
class ReconciledOccurrence:
    block_id: str
    occurrence_id: str
    species: str
    chrom: str
    start: int
    end: int
    strand: str
    copy_id: str
    status: str
    ancestor: str
    anc_chrom: str
    anc_start: int
    anc_end: int
    coverage: str
    source_anchor_id: str
    child_side: str
    child_node: str
    child_block_id: str
    child_occurrence_id: str


@dataclass(frozen=True, slots=True)
class ProvenanceRecord:
    parent_block_id: str
    child_side: str
    child_node: str
    child_block_id: str
    child_occurrence_id: str
    source_anchor_id: str
    relationship: str


@dataclass(frozen=True, slots=True)
class ConflictRecord:
    parent_node: str
    parent_chrom: str
    parent_start: int
    parent_end: int
    conflict_type: str
    left_block_ids: str
    right_block_ids: str
    left_copies: str
    right_copies: str
    explanation: str
    resolution: str


@dataclass(frozen=True, slots=True)
class MergeDecision:
    parent_node: str
    parent_chrom: str
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    can_merge: bool
    accepted_reason: str
    rejection_reasons: tuple[str, ...]
    parent_gap: int
    left_child_gaps: str
    right_child_gaps: str


@dataclass(frozen=True, slots=True)
class ReconcileResult:
    blocks: tuple[ParentBlock, ...]
    occurrences: tuple[ReconciledOccurrence, ...]
    atomic_intervals: tuple[AtomicInterval, ...]
    provenance: tuple[ProvenanceRecord, ...]
    conflicts: tuple[ConflictRecord, ...]
    input_records: int
    splits: int
    merges: int
    warnings: tuple[str, ...]
    merge_decisions: tuple[MergeDecision, ...] = ()
    statistics: dict[str, object] | None = None


@dataclass(slots=True)
class _Atom:
    chrom: str
    start: int
    end: int
    left: tuple[ParentMappedRun, ...]
    right: tuple[ParentMappedRun, ...]
    classification: str
    retained: bool
    block_id: str = ""


@dataclass(frozen=True, slots=True)
class _ProjectedPiece:
    record: ParentMappedRun
    parent_start: int
    parent_end: int
    start: int
    end: int


@dataclass(frozen=True, slots=True)
class _TrackSummary:
    occurrence_keys: frozenset[tuple[str, str, str, str, str]]
    copy_ids_by_species: dict[str, frozenset[str]]
    duplicated_species: frozenset[str]
    repeated_copy_species: frozenset[str]


def _record_sort_key(record: ParentMappedRun) -> tuple[object, ...]:
    return (
        record.parent_start,
        record.parent_end,
        record.child_node,
        record.child_block_id,
        record.child_occurrence_id,
        record.species,
        record.chrom,
        record.start,
        record.end,
        record.copy_id,
        record.source_anchor_id,
    )


def _project_run_to_parent_block(
    record: ParentMappedRun, block_start: int, block_end: int
) -> _ProjectedPiece:
    overlap_start = max(record.parent_start, block_start)
    overlap_end = min(record.parent_end, block_end)
    if overlap_end <= overlap_start:
        raise ValueError("cannot project an empty parent intersection")

    if record.strand == "+":
        child_start = record.start + (overlap_start - record.parent_start)
        child_end = record.start + (overlap_end - record.parent_start)
    else:
        child_start = record.end - (overlap_end - record.parent_start)
        child_end = record.end - (overlap_start - record.parent_start)

    if child_end - child_start != overlap_end - overlap_start:
        raise ValueError("projected child interval length does not match parent overlap")
    return _ProjectedPiece(record, overlap_start, overlap_end, child_start, child_end)


def _order_conflicts(runs: Iterable[ParentMappedRun]) -> set[tuple[str, str, str, str]]:
    grouped: dict[tuple[str, str, str, str], list[ParentMappedRun]] = defaultdict(list)
    for run in runs:
        grouped[(run.child_node, run.child_block_id, run.child_occurrence_id, run.copy_id)].append(run)
    conflicts: set[tuple[str, str, str, str]] = set()
    for key, group in grouped.items():
        ordered = sorted(group, key=lambda r: (r.parent_chrom, r.parent_start, r.parent_end))
        if len({run.parent_chrom for run in ordered}) > 1 or len({run.chrom for run in ordered}) > 1:
            conflicts.add(key)
            continue
        for previous, current in pairwise(ordered):
            same_interval = (
                previous.parent_start == current.parent_start
                and previous.parent_end == current.parent_end
                and previous.start == current.start
                and previous.end == current.end
            )
            if same_interval:
                continue
            monotonic = (
                previous.end <= current.start
                if previous.strand == current.strand == "+"
                else current.end <= previous.start
                if previous.strand == current.strand == "-"
                else False
            )
            if not monotonic:
                conflicts.add(key)
                break
    return conflicts


def _summarize_tracks(records: tuple[ParentMappedRun, ...]) -> _TrackSummary:
    """Summarize records by biological occurrence track.

    Rows sharing child node, species, child block, child occurrence, and copy ID
    are fragments or source anchors for one occurrence. Rows with different
    occurrence identities within one species are separate candidate copies.
    """
    occurrence_keys = frozenset(
        (
            record.child_node,
            record.species,
            record.child_block_id,
            record.child_occurrence_id,
            record.copy_id,
        )
        for record in records
    )
    copy_ids_by_species: dict[str, set[str]] = defaultdict(set)
    species_copy_to_occurrences: dict[
        tuple[str, str], set[tuple[str, str, str, str, str]]
    ] = defaultdict(set)
    duplicated_species: set[str] = set()
    for record in records:
        key = (
            record.child_node,
            record.species,
            record.child_block_id,
            record.child_occurrence_id,
            record.copy_id,
        )
        copy_ids_by_species[record.species].add(record.copy_id)
        species_copy_to_occurrences[(record.species, record.copy_id)].add(key)
        if record.status == "duplicated":
            duplicated_species.add(record.species)
    repeated_copy_species = {
        species
        for (species, _copy_id), keys in species_copy_to_occurrences.items()
        if len(keys) > 1
    }
    return _TrackSummary(
        occurrence_keys,
        {species: frozenset(copy_ids) for species, copy_ids in copy_ids_by_species.items()},
        frozenset(duplicated_species),
        frozenset(repeated_copy_species),
    )


def _has_species_level_duplication(summary: _TrackSummary) -> bool:
    occurrence_counts: dict[str, int] = defaultdict(int)
    for _child_node, species, _block_id, _occurrence_id, _copy_id in summary.occurrence_keys:
        occurrence_counts[species] += 1
    return any(count > 1 for count in occurrence_counts.values()) or bool(
        summary.duplicated_species or summary.repeated_copy_species
    )


def _global_copy_ids_resolve(
    left_summary: _TrackSummary,
    right_summary: _TrackSummary,
) -> bool:
    if left_summary.repeated_copy_species or right_summary.repeated_copy_species:
        return False
    for summary in (left_summary, right_summary):
        if any(len(summary.copy_ids_by_species[species]) < 2 for species in summary.duplicated_species):
            return False
    left_copies = set().union(*left_summary.copy_ids_by_species.values())
    right_copies = set().union(*right_summary.copy_ids_by_species.values())
    if not left_copies or not right_copies:
        return False
    return left_copies == right_copies


def _biological_track_count(occurrences: Iterable[ReconciledOccurrence]) -> int:
    return len(
        {
            (
                occurrence.child_side,
                occurrence.child_node,
                occurrence.species,
                occurrence.child_block_id,
                occurrence.child_occurrence_id,
                occurrence.copy_id,
            )
            for occurrence in occurrences
        }
    )


def _classify(
    left: tuple[ParentMappedRun, ...],
    right: tuple[ParentMappedRun, ...],
    order_conflicts: set[tuple[str, str, str, str]],
    config: ReconcileConfig,
) -> str:
    all_records = left + right
    usable_left = tuple(record for record in left if record.status not in {"ambiguous", "unaligned"})
    usable_right = tuple(record for record in right if record.status not in {"ambiguous", "unaligned"})

    ambiguous = any(record.status == "ambiguous" for record in all_records)
    unaligned = any(record.status == "unaligned" for record in all_records)
    has_usable = bool(usable_left or usable_right)
    if not usable_left and not usable_right:
        if ambiguous and unaligned:
            return "complex"
        if ambiguous:
            return "ambiguous"
        return "unaligned"
    if has_usable and (ambiguous or unaligned):
        return "complex"
    if usable_left and not usable_right:
        return "left_only"
    if usable_right and not usable_left:
        return "right_only"

    conflict_types: list[str] = []
    left_orientations = {record.strand for record in usable_left}
    right_orientations = {record.strand for record in usable_right}
    if left_orientations != right_orientations or len(left_orientations) != 1:
        conflict_types.append("orientation_conflict")
    if any(
        (record.child_node, record.child_block_id, record.child_occurrence_id, record.copy_id)
        in order_conflicts
        for record in usable_left + usable_right
    ):
        conflict_types.append("order_conflict")

    left_summary = _summarize_tracks(usable_left)
    right_summary = _summarize_tracks(usable_right)
    has_candidates = _has_species_level_duplication(left_summary) or _has_species_level_duplication(
        right_summary
    )
    if has_candidates:
        copies_resolve = config.copy_id_scope == "global" and _global_copy_ids_resolve(
            left_summary, right_summary
        )
        if not copies_resolve:
            conflict_types.append("duplication_conflict")

    if len(conflict_types) > 1:
        return "complex"
    return conflict_types[0] if conflict_types else "shared_consistent"


def _sweep_chromosome(
    chrom: str,
    left: list[ParentMappedRun],
    right: list[ParentMappedRun],
    config: ReconcileConfig,
    order_conflicts: set[tuple[str, str, str, str]],
) -> list[_Atom]:
    records = [("left", record) for record in left] + [("right", record) for record in right]
    boundaries = sorted({point for _, record in records for point in (record.parent_start, record.parent_end)})
    starts: dict[int, list[tuple[str, ParentMappedRun]]] = defaultdict(list)
    ends: dict[int, list[tuple[str, ParentMappedRun]]] = defaultdict(list)
    for side, record in records:
        starts[record.parent_start].append((side, record))
        ends[record.parent_end].append((side, record))
    active: dict[str, dict[tuple[str, str, str, str], ParentMappedRun]] = {
        "left": {},
        "right": {},
    }
    atoms: list[_Atom] = []
    for start, end in pairwise(boundaries):
        for side, record in ends[start]:
            active[side].pop(record.record_id, None)
        for side, record in starts[start]:
            active[side][record.record_id] = record
        left_active = tuple(sorted(active["left"].values(), key=_record_sort_key))
        right_active = tuple(sorted(active["right"].values(), key=_record_sort_key))
        if not left_active and not right_active:
            continue
        classification = _classify(left_active, right_active, order_conflicts, config)
        atoms.append(
            _Atom(
                chrom,
                start,
                end,
                left_active,
                right_active,
                classification,
                end - start >= config.min_block_length,
            )
        )
    return atoms


def _membership(atom: _Atom) -> tuple[object, ...]:
    def signature(records: tuple[ParentMappedRun, ...]) -> tuple[tuple[str, ...], ...]:
        return tuple(
            sorted(
                (
                    record.child_block_id,
                    record.child_occurrence_id,
                    record.species,
                    record.chrom,
                    record.strand,
                    record.copy_id,
                    record.status,
                    record.source_anchor_id,
                    str(record.start),
                    str(record.end),
                )
                for record in records
            )
        )

    return atom.classification, signature(atom.left), signature(atom.right)


def merge_decision(previous: _Atom, current: _Atom, config: ReconcileConfig, parent_node: str = "") -> MergeDecision:
    reasons: list[str] = []
    gap = current.start - previous.end
    if previous.chrom != current.chrom:
        reasons.append("different_parent_chromosome")
    if gap > config.max_merge_gap:
        reasons.append("parent_gap_exceeds_threshold")
    if reasons:
        return MergeDecision(parent_node, previous.chrom, previous.start, previous.end, current.start, current.end,
                             False, "", tuple(reasons), gap, _child_gap_text(previous, current, "left"), _child_gap_text(previous, current, "right"))
    if _membership(previous) == _membership(current):
        # Exact record membership makes descendant projections collinear and monotonic.
        ok = all(
            record.parent_start <= previous.start and record.parent_end >= current.end
            for record in previous.left + previous.right
        )
        return MergeDecision(parent_node, previous.chrom, previous.start, previous.end, current.start, current.end,
                             ok, "merge_compatible" if ok else "intervening_filtered_atom", () if ok else ("intervening_filtered_atom",),
                             gap, _child_gap_text(previous, current, "left"), _child_gap_text(previous, current, "right"))
    if (
        previous.classification != current.classification
        or current.start != previous.end
        or min(previous.end - previous.start, current.end - current.start)
        > config.boundary_tolerance
    ):
        if previous.classification != current.classification:
            reasons.append("different_classification")
        if current.start != previous.end:
            reasons.append("intervening_filtered_atom")
        if min(previous.end - previous.start, current.end - current.start) > config.boundary_tolerance:
            reasons.append("boundary_atom_exceeds_tolerance")
        return MergeDecision(parent_node, previous.chrom, previous.start, previous.end, current.start, current.end,
                             False, "", tuple(reasons), gap, _child_gap_text(previous, current, "left"), _child_gap_text(previous, current, "right"))

    def side_is_collinear(
        before: tuple[ParentMappedRun, ...], after: tuple[ParentMappedRun, ...]
    ) -> bool:
        def key(record: ParentMappedRun) -> tuple[str, ...]:
            return (
                record.species,
                record.chrom,
                record.strand,
                record.copy_id,
                record.status,
            )

        before_by_key = {key(record): record for record in before}
        after_by_key = {key(record): record for record in after}
        if len(before_by_key) != len(before) or len(after_by_key) != len(after):
            return False
        if before_by_key.keys() != after_by_key.keys():
            return False
        for record_key, first in before_by_key.items():
            second = after_by_key[record_key]
            if first == second:
                continue
            if first.strand == "+" and first.end != second.start:
                return False
            if first.strand == "-" and second.end != first.start:
                return False
        return True

    left_ok = side_is_collinear(previous.left, current.left)
    right_ok = side_is_collinear(previous.right, current.right)
    if not left_ok:
        reasons.append("non_collinear_forward" if any(r.strand == "+" for r in previous.left + current.left) else "non_collinear_reverse")
    if not right_ok:
        reasons.append("non_collinear_forward" if any(r.strand == "+" for r in previous.right + current.right) else "non_collinear_reverse")
    ok = left_ok and right_ok
    return MergeDecision(parent_node, previous.chrom, previous.start, previous.end, current.start, current.end,
                         ok, "merge_compatible" if ok else "", tuple(dict.fromkeys(reasons)), gap,
                         _child_gap_text(previous, current, "left"), _child_gap_text(previous, current, "right"))


def _child_gap_text(previous: _Atom, current: _Atom, side: str) -> str:
    before = previous.left if side == "left" else previous.right
    after = current.left if side == "left" else current.right
    values = []
    for first in before:
        for second in after:
            if first.chrom == second.chrom and first.strand == second.strand:
                values.append(str(abs(second.start - first.end) if first.strand == "+" else abs(first.start - second.end)))
    return ",".join(sorted(set(values)))


def _can_merge(previous: _Atom, current: _Atom, config: ReconcileConfig) -> bool:
    return merge_decision(previous, current, config).can_merge


def _mapping_quality(records: Iterable[ParentMappedRun]) -> str:
    statuses = {record.status for record in records}
    if "unaligned" in statuses:
        return "unaligned"
    if "ambiguous" in statuses:
        return "ambiguous"
    if "duplicated" in statuses:
        return "duplicated"
    return "unique"


def _relationship(record: ParentMappedRun, group: list[_Atom], atoms: list[_Atom]) -> str:
    if record.status in {"duplicated", "ambiguous", "unaligned"}:
        return record.status

    record_length = record.parent_end - record.parent_start
    filtered_coverage = sum(
        atom.end - atom.start
        for atom in atoms
        if not atom.retained and (record in atom.left or record in atom.right)
    )
    retained_atoms = [
        atom for atom in atoms if atom.retained and (record in atom.left or record in atom.right)
    ]
    retained_coverage = sum(atom.end - atom.start for atom in retained_atoms)
    retained_blocks = {atom.block_id for atom in retained_atoms if atom.block_id}
    group_coverage = sum(
        atom.end - atom.start
        for atom in group
        if record in atom.left or record in atom.right
    )

    if filtered_coverage:
        return "filtered_partial"
    if len(retained_blocks) > 1:
        return "split_partial" if retained_coverage < record_length else "split"
    if len(group) > 1 and sum(record in atom.left or record in atom.right for atom in group) > 1:
        return "merged"
    if group_coverage < record_length:
        return "split_partial"
    return "preserved"


def _project_occurrences(
    block_id: str,
    block_start: int,
    block_end: int,
    records: tuple[ParentMappedRun, ...],
    left_runs: list[ParentMappedRun],
    parent_node: str,
    parent_chrom: str,
) -> list[ReconciledOccurrence]:
    grouped: dict[tuple[object, ...], list[_ProjectedPiece]] = defaultdict(list)
    for record in records:
        side = "left" if record in left_runs else "right"
        piece = _project_run_to_parent_block(record, block_start, block_end)
        key = (
            side,
            record.child_node,
            record.child_block_id,
            record.child_occurrence_id,
            record.species,
            record.chrom,
            record.strand,
            record.copy_id,
            record.status,
        )
        grouped[key].append(piece)

    coalesced: list[tuple[tuple[object, ...], list[_ProjectedPiece]]] = []
    for key, pieces in sorted(grouped.items()):
        ordered = sorted(
            pieces,
            key=lambda piece: (
                piece.parent_start,
                piece.parent_end,
                piece.start,
                piece.end,
                piece.record.source_anchor_id,
            ),
        )
        for piece in ordered:
            if not coalesced or coalesced[-1][0] != key:
                coalesced.append((key, [piece]))
                continue
            previous = coalesced[-1][1][-1]
            same_interval = (
                previous.parent_start == piece.parent_start
                and previous.parent_end == piece.parent_end
                and previous.start == piece.start
                and previous.end == piece.end
            )
            adjacent_in_parent = previous.parent_end == piece.parent_start
            if piece.record.strand == "+":
                adjacent_in_child = previous.end == piece.start
            else:
                adjacent_in_child = piece.end == previous.start
            if same_interval or (adjacent_in_parent and adjacent_in_child):
                coalesced[-1][1].append(piece)
            else:
                coalesced.append((key, [piece]))

    occurrences: list[ReconciledOccurrence] = []
    for occurrence_index, (key, pieces) in enumerate(coalesced, start=1):
        (
            side,
            child_node,
            child_block_id,
            child_occurrence_id,
            species,
            chrom,
            strand,
            copy_id,
            status,
        ) = key
        anchors = ",".join(sorted({piece.record.source_anchor_id for piece in pieces}))
        anc_start = min(piece.parent_start for piece in pieces)
        anc_end = max(piece.parent_end for piece in pieces)
        start = min(piece.start for piece in pieces)
        end = max(piece.end for piece in pieces)
        if end - start != anc_end - anc_start:
            raise ValueError("coalesced occurrence lengths disagree between child and parent")
        coverage = "full" if anc_start == block_start and anc_end == block_end else "partial"
        occurrences.append(
            ReconciledOccurrence(
                block_id,
                f"{block_id}.{occurrence_index}",
                species,
                chrom,
                start,
                end,
                strand,
                copy_id,
                status,
                parent_node,
                parent_chrom,
                anc_start,
                anc_end,
                coverage,
                anchors,
                side,
                child_node,
                child_block_id,
                child_occurrence_id,
            )
        )
    return occurrences


def reconcile_node(
    left_runs: list[ParentMappedRun],
    right_runs: list[ParentMappedRun],
    *,
    parent_node: str,
    left_node: str,
    right_node: str,
    config: ReconcileConfig | None = None,
) -> ReconcileResult:
    """Reconcile two child systems already mapped into one common parent."""
    config = config or ReconcileConfig()
    all_runs = left_runs + right_runs
    if any(run.parent_node != parent_node for run in all_runs):
        found = sorted({run.parent_node for run in all_runs})
        raise ValueError(f"all runs must use requested parent node {parent_node!r}; found {found}")
    if any(run.child_node != left_node for run in left_runs):
        raise ValueError(f"left input contains records not from left node {left_node!r}")
    if any(run.child_node != right_node for run in right_runs):
        raise ValueError(f"right input contains records not from right node {right_node!r}")
    record_ids = [run.record_id for run in all_runs]
    if len(record_ids) != len(set(record_ids)):
        raise ValueError("duplicate record identity (child node/block/occurrence/source anchor)")

    order_conflicts = _order_conflicts(all_runs)
    by_chrom: dict[str, dict[str, list[ParentMappedRun]]] = defaultdict(
        lambda: {"left": [], "right": []}
    )
    for run in left_runs:
        by_chrom[run.parent_chrom]["left"].append(run)
    for run in right_runs:
        by_chrom[run.parent_chrom]["right"].append(run)
    atoms: list[_Atom] = []
    for chrom in sorted(by_chrom):
        sides = by_chrom[chrom]
        atoms.extend(
            _sweep_chromosome(chrom, sides["left"], sides["right"], config, order_conflicts)
        )

    groups: list[list[_Atom]] = []
    decisions: list[MergeDecision] = []
    hit_barrier = False
    for atom in atoms:
        if not atom.retained:
            hit_barrier = True
            continue
        if groups and hit_barrier:
            previous = groups[-1][-1]
            base = merge_decision(previous, atom, config, parent_node)
            decisions.append(MergeDecision(base.parent_node, base.parent_chrom, base.left_start, base.left_end,
                base.right_start, base.right_end, False, "", tuple(dict.fromkeys(base.rejection_reasons + ("intervening_filtered_atom",))),
                base.parent_gap, base.left_child_gaps, base.right_child_gaps))
            groups.append([atom])
        elif groups:
            decision = merge_decision(groups[-1][-1], atom, config, parent_node)
            decisions.append(decision)
            if decision.can_merge:
                groups[-1].append(atom)
            else:
                groups.append([atom])
        else:
            groups.append([atom])
        hit_barrier = False

    blocks: list[ParentBlock] = []
    occurrences: list[ReconciledOccurrence] = []
    provenance: list[ProvenanceRecord] = []
    for index, group in enumerate(groups, start=1):
        block_id = f"HMSP{index:08d}"
        for atom in group:
            atom.block_id = block_id

    for index, group in enumerate(groups, start=1):
        block_id = f"HMSP{index:08d}"
        records = tuple(
            sorted({record for atom in group for record in atom.left + atom.right}, key=_record_sort_key)
        )
        occurrences.extend(
            _project_occurrences(
                block_id,
                group[0].start,
                group[-1].end,
                records,
                left_runs,
                parent_node,
                group[0].chrom,
            )
        )
        block_occurrences = [occurrence for occurrence in occurrences if occurrence.block_id == block_id]
        blocks.append(
            ParentBlock(
                block_id,
                parent_node,
                group[0].chrom,
                group[0].start,
                group[-1].end,
                group[0].classification,
                _ids((record for atom in group for record in atom.left), "child_block_id"),
                _ids((record for atom in group for record in atom.right), "child_block_id"),
                len({occurrence.species for occurrence in block_occurrences}),
                len(block_occurrences),
                _biological_track_count(block_occurrences),
                _mapping_quality(records),
                len(group),
            )
        )
        for record in records:
            side = "left" if record in left_runs else "right"
            provenance.append(
                ProvenanceRecord(
                    block_id,
                    side,
                    record.child_node,
                    record.child_block_id,
                    record.child_occurrence_id,
                    record.source_anchor_id,
                    _relationship(record, group, atoms),
                )
            )

    atomic_rows: list[AtomicInterval] = []
    for index, atom in enumerate(atoms, start=1):
        merged = bool(atom.block_id and next(b for b in blocks if b.block_id == atom.block_id).atomic_interval_count > 1)
        atomic_rows.append(
            AtomicInterval(
                f"ATOM{index:08d}",
                parent_node,
                atom.chrom,
                atom.start,
                atom.end,
                _ids(atom.left, "source_anchor_id"),
                _ids(atom.right, "source_anchor_id"),
                _ids(atom.left, "child_block_id"),
                _ids(atom.right, "child_block_id"),
                atom.classification,
                _ids(atom.left, "strand"),
                _ids(atom.right, "strand"),
                _ids(atom.left, "copy_id"),
                _ids(atom.right, "copy_id"),
                "filtered" if not atom.retained else "merged" if merged else "retained",
                "below_min_block_length"
                if not atom.retained
                else "identical evidence; collinear monotonic projections"
                if merged
                else next((d.rejection_reasons[0] if d.rejection_reasons else d.accepted_reason
                           for d in decisions if d.right_start == atom.start and d.left_end == atom.start), "not_merged"),
                atom.block_id,
            )
        )

    explanations = {
        "orientation_conflict": "child support has incompatible parent-relative orientations",
        "order_conflict": "a child occurrence is non-collinear, non-monotonic, or changes chromosome",
        "duplication_conflict": "candidate copies cannot be paired uniquely by copy identity",
        "ambiguous": "input mapping is explicitly ambiguous",
        "unaligned": "one or both sides have explicit unaligned evidence or no usable mapping",
        "complex": "multiple conflict conditions occur in the same atomic interval",
    }
    conflicts = tuple(
        ConflictRecord(
            parent_node,
            atom.chrom,
            atom.start,
            atom.end,
            atom.classification,
            _ids(atom.left, "child_block_id"),
            _ids(atom.right, "child_block_id"),
            _ids(atom.left, "copy_id"),
            _ids(atom.right, "copy_id"),
            explanations[atom.classification],
            "unresolved",
        )
        for atom in atoms
        if atom.classification in explanations
    )
    split_records = sum(
        1
        for record in all_runs
        if sum(record in atom.left or record in atom.right for atom in atoms if atom.retained) > 1
    )
    warnings = (
        (
            "boundary tolerance is used only when evaluating near-boundary compatibility; "
            "the exact atomic partition is unchanged"
        ),
        "halSynteny transitivity has not been established; this run records provenance only",
    )
    classification_counts = {name: sum(1 for atom in atoms if atom.classification == name) for name in CLASSIFICATIONS}
    classification_bp = {name: sum(atom.end - atom.start for atom in atoms if atom.classification == name) for name in CLASSIFICATIONS}
    retained = sum(atom.retained for atom in atoms)
    mapped_blocks = len({(run.child_node, run.child_block_id) for run in all_runs})
    statistics = {
        "left_input_block_count": len({run.child_block_id for run in left_runs}),
        "right_input_block_count": len({run.child_block_id for run in right_runs}),
        "left_mapped_child_block_count": len({run.child_block_id for run in left_runs}),
        "right_mapped_child_block_count": len({run.child_block_id for run in right_runs}),
        "left_mapping_run_count": len(left_runs), "right_mapping_run_count": len(right_runs),
        "left_unmapped_child_block_count": 0, "right_unmapped_child_block_count": 0,
        "distinct_input_boundaries": len({point for run in all_runs for point in (run.parent_start, run.parent_end)}),
        "covered_atomic_intervals": len(atoms), "below_minimum_length_atoms": sum(not a.retained for a in atoms),
        "retained_atoms_before_merging": retained, "final_parent_blocks": len(blocks),
        "accepted_pairwise_merges": sum(d.can_merge for d in decisions),
        "compression_ratio": len(blocks) / retained if retained else 0.0,
        "fragmentation_expansion": retained / mapped_blocks if mapped_blocks else 0.0,
        "classification_counts": classification_counts, "classification_bp": classification_bp,
        "unique_mappings": sum(run.status == "unique" for run in all_runs),
        "duplicated_mappings": sum(run.status == "duplicated" for run in all_runs),
        "ambiguous_mappings": sum(run.status == "ambiguous" for run in all_runs),
        "unaligned_mappings": sum(run.status == "unaligned" for run in all_runs),
        "leaf_occurrences_represented": len(occurrences),
        "merge_rejection_counts": {reason: sum(reason in d.rejection_reasons for d in decisions) for reason in sorted({r for d in decisions for r in d.rejection_reasons})},
    }
    return ReconcileResult(
        tuple(blocks),
        tuple(occurrences),
        tuple(atomic_rows),
        tuple(sorted(provenance, key=lambda row: (
            row.parent_block_id, row.child_side, row.child_node, row.child_block_id,
            row.child_occurrence_id, row.source_anchor_id,
        ))),
        conflicts,
        len(all_runs),
        split_records,
        sum(max(0, len(group) - 1) for group in groups),
        warnings,
        tuple(decisions),
        statistics,
    )


def checksum(path_bytes: bytes) -> str:
    """Return a stable SHA-256 checksum for summary metadata."""
    return sha256(path_bytes).hexdigest()
