# Piecewise Palaeognathae comparison

This compact comparison reclassifies a fresh direct-edge `halLiftover --outPSLWithName` capture. Its bytes differ from the prior capture only in PSL row order; normalized records, coordinates, strands, and components are identical. Runtime and peak memory were not measured for this capture. No halSynteny command was run.

Edge: GCF_028389875.1 -> PalaeognathaeAnc2

| backend | mapped | unmapped | covered bp | uncovered bp | raw unique | raw ambiguous | effective unique | effective duplicated | effective ambiguous | parent blocks | runs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| hal-psl-conservative | 8 | 4 | 2339 | 747 | 4 | 4 | 4 | 0 | 4 | 59 | 59 |
| hal-psl-piecewise | 8 | 4 | 2339 | 747 | 8 | 0 | 8 | 0 | 0 | 59 | 59 |

The four previously ambiguous blocks are in real-comparison.tsv. Source-disjoint records are geometric fragments and do not establish biological synteny or copy correspondence.
