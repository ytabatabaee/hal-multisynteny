#!/usr/bin/env bash
set -euo pipefail

out_dir=${1:-tiny-hal-output}
repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixture_dir="${repo_dir}/tests/fixtures/tiny_hal"
hal="${fixture_dir}/build/tiny.hal"

if [[ ! -s "${hal}" ]]; then
  "${repo_dir}/scripts/build_tiny_hal_fixture.sh"
fi

rm -rf "${out_dir}"
mkdir -p "${out_dir}/audit"

hal-multisynteny hal-info \
  --hal "${hal}" \
  --metadata-level basic \
  --output "${out_dir}/hal-info.json"

python - "${fixture_dir}/manifest.json" <<'PY' | while IFS=$'\t' read -r child parent label blocks; do
import json
import sys
from pathlib import Path
manifest = json.loads(Path(sys.argv[1]).read_text())
for edge in manifest["audit_edges"]:
    print("\t".join([edge["child"], edge["parent"], edge["label"], edge["blocks"]]))
PY
  hal-multisynteny audit-liftover \
    --hal "${hal}" \
    --child-genome "${child}" \
    --parent-genome "${parent}" \
    --blocks "${fixture_dir}/${blocks}" \
    --output-prefix "${out_dir}/audit/${label}"
done

hal-multisynteny run-tree \
  --hal "${hal}" \
  --tree "${fixture_dir}/tree.nwk" \
  --node-map "${fixture_dir}/node-map.tsv" \
  --leaf-blocks "${fixture_dir}/leaf-blocks-traversal" \
  --output-dir "${out_dir}/run-a" \
  --min-block-length 1 \
  --backend hal

hal-multisynteny run-tree \
  --hal "${hal}" \
  --tree "${fixture_dir}/tree.nwk" \
  --node-map "${fixture_dir}/node-map.tsv" \
  --leaf-blocks "${fixture_dir}/leaf-blocks-traversal" \
  --output-dir "${out_dir}/run-a" \
  --min-block-length 1 \
  --backend hal \
  --resume

hal-multisynteny run-tree \
  --hal "${hal}" \
  --tree "${fixture_dir}/tree.nwk" \
  --node-map "${fixture_dir}/node-map.tsv" \
  --leaf-blocks "${fixture_dir}/leaf-blocks-traversal" \
  --output-dir "${out_dir}/run-a" \
  --min-block-length 1 \
  --backend hal \
  --resume \
  --force-node ancAB

hal-multisynteny run-tree \
  --hal "${hal}" \
  --tree "${fixture_dir}/tree.nwk" \
  --node-map "${fixture_dir}/node-map.tsv" \
  --leaf-blocks "${fixture_dir}/leaf-blocks-traversal" \
  --output-dir "${out_dir}/run-b" \
  --min-block-length 1 \
  --backend hal

cmp "${out_dir}/run-a/nodes/root/blocks.tsv" "${out_dir}/run-b/nodes/root/blocks.tsv"
cmp "${out_dir}/run-a/nodes/root/node_occurrences.tsv" "${out_dir}/run-b/nodes/root/node_occurrences.tsv"
cmp "${out_dir}/run-a/nodes/root/leaf_occurrences.tsv" "${out_dir}/run-b/nodes/root/leaf_occurrences.tsv"
