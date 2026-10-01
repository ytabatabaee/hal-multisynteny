#!/usr/bin/env bash
set -euo pipefail

out_dir=${1:-tiny-hal-output}
repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixture_dir="${repo_dir}/tests/fixtures/tiny_hal"
hal="${fixture_dir}/build/tiny.hal"

echo "HAL tool versions:"
for tool in maf2hal halStats halLiftover; do
  command -v "${tool}"
  "${tool}" --version 2>&1 | head -n 1 || true
done

"${repo_dir}/scripts/build_tiny_hal_fixture.sh"

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

python "${repo_dir}/scripts/validate_tiny_hal_outputs.py" \
  --fixture-dir "${fixture_dir}" \
  --output-dir "${out_dir}"

hal-multisynteny run-tree \
  --hal "${hal}" \
  --tree "${fixture_dir}/tree.nwk" \
  --node-map "${fixture_dir}/node-map.tsv" \
  --leaf-blocks "${fixture_dir}/leaf-blocks-traversal" \
  --output-dir "${out_dir}/run-a" \
  --min-block-length 1 \
  --backend hal \
  --resume

python - "${out_dir}/run-a/run-summary.json" <<'PY'
import json
import sys
summary = json.loads(open(sys.argv[1], encoding="utf-8").read())
assert summary["reused_nodes"] == ["ancAB", "ancCD", "root"], summary
print("resume reused:", ",".join(summary["reused_nodes"]))
PY

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

python - "${out_dir}/run-a/run-summary.json" <<'PY'
import json
import sys
summary = json.loads(open(sys.argv[1], encoding="utf-8").read())
assert summary["recomputed_nodes"] == ["ancAB", "root"], summary
assert summary["reused_nodes"] == ["ancCD"], summary
print("force-node recomputed:", ",".join(summary["recomputed_nodes"]), "reused:", ",".join(summary["reused_nodes"]))
PY

hal-multisynteny run-tree \
  --hal "${hal}" \
  --tree "${fixture_dir}/tree.nwk" \
  --node-map "${fixture_dir}/node-map.tsv" \
  --leaf-blocks "${fixture_dir}/leaf-blocks-traversal" \
  --output-dir "${out_dir}/run-b" \
  --min-block-length 1 \
  --backend hal

meaningful=(
  blocks.tsv
  node_occurrences.tsv
  leaf_occurrences.tsv
  left_edge_runs.tsv
  right_edge_runs.tsv
  left_unmapped_edge_evidence.tsv
  right_unmapped_edge_evidence.tsv
  provenance.tsv
  conflicts.tsv
  summary.json
)
for node in ancAB ancCD root; do
  for name in "${meaningful[@]}"; do
    a="${out_dir}/run-a/nodes/${node}/${name}"
    b="${out_dir}/run-b/nodes/${node}/${name}"
    if [[ -e "${a}" || -e "${b}" ]]; then
      cmp "${a}" "${b}"
    fi
  done
done

echo "tiny HAL integration passed: audits, semantic checkpoints, resume, force-node, and byte-stable scientific outputs"
