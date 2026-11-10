# Preserved Palaeognathae PSL multiplicity audit

Audited 24,000 input intervals across 4 preserved direct-edge captures.
No HAL commands were executed. Coverage is limited to edges with existing PSL evidence; BED6-only edges are excluded.

Source atoms are disjoint after atomization; target projections may overlap or be one-to-many. Atomwise unique does not establish biological uniqueness or synteny.

| Child → parent | Inputs | Source multiplicity >1 | Newly ambiguous | Same-target overlap pairs | Max multiplicity |
|---|---:|---:|---:|---:|---:|
| GCF_028389875.1 → PalaeognathaeAnc2 | 6000 | 0 | 76 | 5794 | 1 |
| PalaeognathaeAnc1 → PalaeognathaeAnc0 | 6000 | 0 | 139 | 6036 | 1 |
| PalaeognathaeAnc2 → PalaeognathaeAnc1 | 6000 | 0 | 54 | 3433 | 1 |
| PalaeognathaeAnc3 → PalaeognathaeAnc2 | 6000 | 0 | 30 | 999 | 1 |

New ambiguity affects 299 inputs. Per-edge concentration (share of all problematic inputs; within-edge rate):
- PalaeognathaeAnc1 → PalaeognathaeAnc0: 46.49% of problematic inputs; 2.32% within edge.
- GCF_028389875.1 → PalaeognathaeAnc2: 25.42% of problematic inputs; 1.27% within edge.
- PalaeognathaeAnc2 → PalaeognathaeAnc1: 18.06% of problematic inputs; 0.90% within edge.
- PalaeognathaeAnc3 → PalaeognathaeAnc2: 10.03% of problematic inputs; 0.50% within edge.

Planned direct edges without preserved PSL (excluded from this audit):
- GCA_047922985.1 → PalaeognathaeAnc1
- GCF_036370855.1 → PalaeognathaeAnc3
- GCF_036417845.1 → PalaeognathaeAnc3
- GCF_040807025.1 → PalaeognathaeAnc0

Every input with source multiplicity >1 (0 inputs) is listed in `palaeognathae.multiplicity.tsv` and the JSON summary.
Inherited ambiguity is unavailable for these randomly sampled intervals because no upward state was captured. Newly observed source multiplicity and same-target overlap are reported separately; ancestor names do not imply inherited ambiguity.
Raw PSL checksums are verified against capture provenance. Sampling-file checksums are recorded now; the historical sampling manifest did not checksum its outputs.
Statistics describe sampled mappings, not full-genome rates, biological copy identities, or liftover transitivity.
