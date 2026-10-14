# HAL PSL diagnostic semantics

This diagnostic is an engineering investigation, not the production liftover
backend. It does not resolve paralogous alternatives versus rearranged mapping
pieces. The eight-edge BED6 audit and reconciliation semantics are independent.

## Verified implementation

The source review uses HAL revision
`acdd8586885583d6c5e5e8e383a7f15c494a1afc`:

- [Liftover::assignBlocksToIntervals, compatible, and flipBlocks](https://github.com/ComparativeGenomicsToolkit/hal/blob/acdd8586885583d6c5e5e8e383a7f15c494a1afc/liftover/impl/halLiftover.cpp)
  sort source blocks, start a new output record for overlapping source blocks or
  incompatible strand, target chromosome, or target order, and order arrays for
  PSL output. Consequently, non-collinearity within one record is noncanonical
  input, not evidence for a rearrangement emitted by HAL.
- [BedLine::writePSL and validatePSL](https://github.com/ComparativeGenomicsToolkit/hal/blob/acdd8586885583d6c5e5e8e383a7f15c494a1afc/liftover/impl/halBedLine.cpp)
  write the query sign followed by the target sign, forward enclosing bounds,
  and strand-relative component starts. Their validation checks base totals,
  enclosing extrema, and inserted-base totals.
- [PSLInfo](https://github.com/ComparativeGenomicsToolkit/hal/blob/acdd8586885583d6c5e5e8e383a7f15c494a1afc/liftover/inc/halBedLine.h)
  stores query strand and absolute query starts before the writer converts them.
- The [UCSC PSL specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format2)
  defines ordinary DNA PSL's one-character query orientation and negative query
  coordinates. HAL additionally writes both strand signs for DNA. This parser
  is limited to nucleotide coordinates, not translated protein PSL.

`halLiftover --help` documents `--outPSLWithName`: the input BED name precedes
each PSL row. A row is an alignment output record, not a biological chain. The
line-based `psl_record_id` carries no alternative/family identity.

## Strand and coordinate rules

| Raw strand | Query sign | Explicit target sign | Relative orientation |
|---|---|---|---|
| `+` | `+` | absent | `+` |
| `-` | `-` | absent | `-` |
| `++` | `+` | `+` | `+` |
| `+-` | `+` | `-` | `-` |
| `-+` | `-` | `+` | `-` |
| `--` | `-` | `-` | `+` |

The target coordinate frame for ordinary one-character DNA PSL is forward;
its *explicit* target sign remains null (empty TSV cell). Mapping orientation
is the query sign in this case. With two signs it is positive when the signs
agree and negative when they differ. A query or target strand sign is never
substituted for relative orientation.

`qStart/qEnd` and `tStart/tEnd` are forward, zero-based, half-open bounds. For
negative query/target signs normalize the corresponding component start as
`sequence_size - strand_start - block_size`; a positive sign leaves it intact.
Both raw component-start arrays must ascend without overlaps in their own strand
coordinate frames. Their normalized extrema must match enclosing bounds. The
parser rejects noncanonical order, overlap, bad arrays/bounds, inconsistent
aligned-base counts, and inconsistent insertion counts. It does not reinterpret
an illegal single record as a rearrangement.

After validation components are sorted in normalized source-forward order;
`component_order` describes that diagnostic order. For negative query signs this
reverses the raw array order. `source_gap_before` is the difference between
adjacent normalized components in the same record; it is empty for each first
component. Raw filename plus actual PSL line number uniquely refers back to the
preserved evidence.

## Classification and observations

The primary categories are `unmapped`,
`single_record_collinear_full_coverage`,
`single_record_collinear_with_source_gaps`, `multiple_records_unresolved`, and
`invalid_or_ambiguous`. A partial single component uses the source-gaps category;
it is not described as split. A syntactically valid row whose source does not
match the sampled input is invalid/ambiguous; malformed records cause an
explicit parser error rather than silently disappearing.

Multiple records always remain unresolved. Source overlap, mutually disjoint
source coverage, different target chromosomes, and different *relative*
orientations are independent observations, not biological classifications.
Joint collinearity requires one target chromosome and relative orientation,
non-overlapping source components, and directionally non-overlapping target
components. It does not establish alternative/path grouping.

Leading, internal, and trailing uncovered bases are measured on the union of
normalized source coverage within the sampled interval. If no component maps,
all input bases are counted as leading uncovered by convention. Unmapped does
not prove biological deletion or non-homology.

## Compact summary schemas and denominators

`palaeognathae-source-aware-by-edge.tsv` and `...-by-bin.tsv` contain only
`edge/scope/length_bin/classification/count/fraction`. Each edge has 6,000 input
intervals; each edge/bin has 1,000. Counts sum exactly to that denominator and
fractions use it, including unmapped intervals.

`...-metrics-by-edge.tsv` and `...-metrics-by-bin.tsv` have exactly one row per
edge or edge/bin. `input_interval_count` is the denominator for every
`intervals_with_*` field. Zero/one/multiple PSL-record interval counts reconcile
to input intervals. `psl_record_count` sums emitted records; `component_count`
sums alignment components, not occurrences or loci. Interval flag counts do not
count record pairs. Disjoint-record coverage means at least two records with no
cross-record source overlap; overlapping means at least one overlapping pair.
These flag counts need not sum to all inputs.

`input_source_bases`, `mapped_source_bases`, `uncovered_source_bases`, and
`excess_source_coverage_bases` are base totals, summed independently for each
sampled interval. Mapped coverage is a source union; excess is summed component
length minus that union, without a biological duplication interpretation.
Mapped plus uncovered equals input bases. Uncovered bases include completely
unmapped inputs. Because sampled intervals may overlap, base totals are not
unique bases across the genome. Aggregate fields are never repeated on
classification rows. Evidence JSON records raw PSL and sampled-input hashes.

Regenerate these summaries from the existing ignored raw PSL and sampled files:

```bash
python -m scripts.regenerate_palaeognathae_source_diagnostic
```

This command does not execute HAL, resample, or rerun the BED6 audit. It checks
raw PSL against its original checksum and HAL identity against the sampling
manifest, then records sampled-input checksums and regenerated-output hashes.

## Controlled real HAL validation

`tests/fixtures/tiny_hal/psl-orientation.maf` is a separate two-genome MAF,
imported with the same `maf2hal` infrastructure as the existing fixture. It has
a forward block, three reversed collinear components, a source-gapped forward
mapping, and a target sequence change. The test executes:

```bash
maf2hal tests/fixtures/tiny_hal/psl-orientation.maf orientation.hal \
  --refGenome PSL_parent --targetGenomes PSL_child
halLiftover --outPSLWithName orientation.hal PSL_child \
  tests/fixtures/tiny_hal/psl-orientation.bed PSL_parent orientation.psl
```

Negative BED inputs exercise both query signs. Compact raw PSL output is
preserved in `psl-orientation.expected.psl`; real integration tests compare all
records and check raw arrays, normalized coordinates, relative orientation,
primary classification, and multi-record observations. Parser-only CI also
checks this captured output. No controlled paralogous duplication is claimed:
this MAF does not encode alternative evolutionary mapping paths, so the
unresolved grouping problem persists. **Readiness: NO-GO** for production
reconstruction based on this diagnostic.
