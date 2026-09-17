from hal_multisynteny.builder import BuildConfig, build_blocks
from hal_multisynteny.models import AlignmentRun


def run(anchor, species, start, end, anc_start, anc_end, strand="+", **kwargs):
    return AlignmentRun(
        anchor_id=anchor,
        species=species,
        chrom="chr1",
        start=start,
        end=end,
        ancestor="avianAncestor",
        anc_chrom="anc1",
        anc_start=anc_start,
        anc_end=anc_end,
        strand=strand,
        **kwargs,
    )


def test_two_species_shared_block():
    result = build_blocks(
        [run("a", "chicken", 100, 200, 0, 100), run("b", "eagle", 500, 600, 0, 100)],
        BuildConfig(min_block_length=10, min_species=2),
    )
    assert len(result.blocks) == 1
    assert result.blocks[0].species_count == 2
    assert {(x.species, x.start, x.end) for x in result.occurrences} == {
        ("chicken", 100, 200),
        ("eagle", 500, 600),
    }


def test_boundaries_create_atomic_blocks():
    result = build_blocks(
        [run("a", "A", 0, 100, 0, 100), run("b", "B", 100, 150, 25, 75)],
        BuildConfig(min_block_length=1, min_species=1),
    )
    assert [(b.anc_start, b.anc_end, b.species_count) for b in result.blocks] == [
        (0, 25, 1),
        (25, 75, 2),
        (75, 100, 1),
    ]


def test_reverse_projection():
    result = build_blocks(
        [run("a", "A", 1000, 1100, 0, 100, strand="-")],
        BuildConfig(min_block_length=1, min_species=1),
    )
    occurrence = result.occurrences[0]
    assert (occurrence.start, occurrence.end, occurrence.strand) == (1000, 1100, "-")


def test_duplicate_copies_are_preserved():
    result = build_blocks(
        [
            run("a1", "A", 0, 100, 0, 100, status="duplicated", copy_id="1"),
            run("a2", "A", 200, 300, 0, 100, status="duplicated", copy_id="2"),
        ],
        BuildConfig(min_block_length=1, min_species=1),
    )
    assert len(result.occurrences) == 2
    assert {x.copy_id for x in result.occurrences} == {"1", "2"}


def test_ambiguous_runs_excluded_by_default():
    result = build_blocks(
        [run("a", "A", 0, 100, 0, 100, status="ambiguous")],
        BuildConfig(min_block_length=1, min_species=1),
    )
    assert result.blocks == ()


def test_empty_input():
    assert build_blocks([]).blocks == ()
