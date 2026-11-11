# Bounded Palaeognathae native-vs-halSynteny comparison

Status: **blocked before reconstruction**. No Palaeognathae traversal output is
claimed by this report.

## Candidate subtree

The smallest verified sibling subtree remains:

```text
(GCF_036417845.1, GCF_036370855.1)PalaeognathaeAnc3;
```

The HAL is `data/cache/palaeognathae/PalaeognathaeAnc0.hal` (7,470,552,826
bytes; SHA-256
`33d29439caed203523fcb2435328a6fb9f0718d72e52c78f7142c84b406df831`). The
direct parent for both leaves is `PalaeognathaeAnc3`; the guide tree and node
map are the committed files in `experiments/downward-pilot/`.

## Selection attempt and blocker

The required fixed regions must be selected from direct native child→Anc3
evidence, not from the earlier arbitrary pilot intervals. The exact commands
attempted were:

```bash
env PATH="$PWD/.hal-env/bin:$PATH" \
  halAlignedExtract data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  GCF_036417845.1 --viewParentCoords > /tmp/pal-native-select/left.tsv

env PATH="$PWD/.hal-env/bin:$PATH" \
  halAlignedExtract data/cache/palaeognathae/PalaeognathaeAnc0.hal \
  GCF_036370855.1 --viewParentCoords > /tmp/pal-native-select/right.tsv
```

`halAlignedExtract v2.2` exposes no chromosome or coordinate restriction. The
first complete direct edge produced approximately 595 MB of text (about 8.3
million rows). The second exceeded 2 GB while still scanning the complete HAL
edge and was interrupted. A bounded streaming summary attempt likewise could
not finish in the available execution window; no region was selected from
partial output. No native or halSynteny reconstruction was run, and no
arbitrary intervals were reused.

This is a precise tooling blocker: the current native backend is streaming and
disk-cached, but its source CLI still requires a complete child-edge scan before
the requested region can be isolated. Selecting a genuinely homologous,
reproducible region therefore requires a new bounded preprocessing/indexing
stage (or a HAL API exposing region-restricted aligned segments).

The new SQLite index builder was then exercised directly on the first edge. The
command ran for 204.35 seconds before being interrupted during the bulk insert;
the temporary database had reached approximately 3.5 GB and was removed. No
final index or segment count was recorded, so this remains an incomplete build,
not usable evidence. Secondary indexes are now created only after bulk loading
to reduce construction overhead on future attempts.

## Required next step

Implement a resumable disk-backed native segment index that scans each required
direct edge once, records chromosome/parent-coordinate offsets, and supports
region queries without loading the edge into RAM. Then select fixed overlapping
Anc3 intervals for both leaves, run Method 3 and Method 1 on identical regions,
and compare their compact checkpoints. Until that stage exists, any claimed
Palaeognathae native-vs-halSynteny comparison would violate the region-selection
and evidence requirements.

The runner now accepts `--native-index-dir`, `--region-chrom`,
`--region-start`, `--region-end`, and `--region-node`; these options are
validated and exercised on controlled indexed fixtures. The current bounded
runner safely supports a selected guide-tree root (the two-leaf Anc3 pilot
shape) and rejects deeper non-root targets rather than reusing root coordinates
incorrectly. Full multi-level region propagation remains a follow-up extension.
The options require a complete validated index and therefore do not bypass the
blocker above.
