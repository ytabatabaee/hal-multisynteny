#!/usr/bin/env bash
set -euo pipefail

HAL="${HAL:-data/cache/palaeognathae/PalaeognathaeAnc0.hal}"
INDEX_DIR="${INDEX_DIR:-/tmp/pal-indexes}"
OUT="${OUT:-/tmp/anchor-30mb}"
PROFILE="${PROFILE:-/tmp/anchor-30mb-profile.json}"
CHROM="${CHROM:-PalaeognathaeAnc3refChr301}"
END="${END:-30374959}"

mkdir -p "$OUT"
env PATH="$PWD/.hal-env/bin:$PATH" \
  .venv/bin/hal-multisynteny run-tree \
  --backend hal-native-segments --leaf-init native \
  --hal "$HAL" --native-index-dir "$INDEX_DIR" \
  --tree experiments/downward-pilot/tree.nwk \
  --node-map experiments/downward-pilot/node-map.tsv \
  --output-dir "$OUT" --reconciliation-policy anchor-chaining \
  --gap-tolerant --max-merge-gap 1000 --max-child-gap 1000 \
  --min-block-length 1 --compact-ancestry --profile \
  --profile-output "$PROFILE" --region-node PalaeognathaeAnc3 \
  --region-chrom "$CHROM" --region-start 0 --region-end "$END"
