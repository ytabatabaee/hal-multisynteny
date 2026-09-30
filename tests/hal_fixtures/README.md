# Tiny real-HAL fixture protocol

The default test suite does not require HAL tools or HAL files. Real-HAL
integration coverage is manifest-driven so small local fixtures can be tested
without committing genome data.

Set `HAL_MULTISYNTENY_TINY_HAL_MANIFEST` to a JSON file with this shape:

```json
{
  "hal": "/path/to/tiny.hal",
  "cases": [
    {
      "name": "identity",
      "child_node_id": "child",
      "parent_node_id": "parent",
      "child_hal_genome": "HAL_CHILD",
      "parent_hal_genome": "HAL_PARENT",
      "block_id": "B1",
      "chrom": "chr1",
      "start": 0,
      "end": 10,
      "expect": "spatial",
      "statuses": ["unique"],
      "strands": ["+"]
    }
  ]
}
```

Recommended tiny fixtures and assertions:

- identity: one exact whole-interval output, `statuses=["unique"]`,
  `strands=["+"]`;
- inversion: one exact whole-interval output, `statuses=["unique"]`,
  `strands=["-"]`;
- deletion/unmapped: no output for the source interval, `expect="unaligned"`;
- one-to-many duplication: multiple full-length distinct outputs,
  `statuses=["duplicated", "duplicated"]`;
- gapped mapping without source subinterval correspondence:
  `expect_error="shorter or gapped"`.

These fixtures should be generated locally from synthetic sequences only. Do not
commit large or copyrighted HAL files. The fake-executable tests remain separate
from this real-HAL path.
