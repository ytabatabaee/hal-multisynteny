# Fake four-leaf pilot

This example runs the 0.3.0 bottom-up tree runner without HAL tools. It is a
technical fixture for checkpointing, propagation, and deterministic traversal;
it is not a biological block definition.

```bash
hal-multisynteny run-tree \
  --tree examples/pilot_tree/pilot.nwk \
  --node-map examples/pilot_tree/node-map.tsv \
  --leaf-blocks examples/pilot_tree/leaf-blocks \
  --output-dir pilot-run \
  --min-block-length 1 \
  --backend fake \
  --fake-mappings examples/pilot_tree/fake-edge-mappings.tsv
```

The final root package is written to `pilot-run/nodes/root/`. Repeating the
same command with `--resume` reuses complete checkpoints only when manifests and
checksums still match.
