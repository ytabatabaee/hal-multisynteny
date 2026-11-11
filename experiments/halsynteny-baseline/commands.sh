#!/usr/bin/env bash
set -euo pipefail

HAL="${HAL:-data/cache/palaeognathae/PalaeognathaeAnc0.hal}"
OUT="${OUT:-halsynteny-output}"
mkdir -p "$OUT"

run_edge() {
  local child="$1"; local tmp="$OUT/${child}.psl.tmp"; local final="$OUT/${child}.psl"
  if [[ -s "$final" ]]; then echo "reusing complete $final"; return; fi
  rm -f "$tmp"
  /usr/bin/time -v env PATH="${HAL_ENV:-$PWD/.hal-env/bin}:$PATH" \
    halSynteny "$HAL" "$tmp" --queryGenome "$child" \
    --targetGenome PalaeognathaeAnc3 --minBlockSize 1 --maxAnchorDistance 1000
  [[ -s "$tmp" ]] || { echo "halSynteny produced no complete PSL for $child" >&2; exit 2; }
  mv "$tmp" "$final"
}

run_edge GCF_036417845.1
run_edge GCF_036370855.1

for child in GCF_036417845.1 GCF_036370855.1; do
  for end in 100000 1000000 10000000; do
    PYTHONPATH=src python scripts/summarize_halsynteny_psl.py \
      --psl "$OUT/$child.psl" --child "$child" --parent PalaeognathaeAnc3 \
      --chrom PalaeognathaeAnc3refChr301 --start 0 --end "$end" \
      --output "$OUT/${child}.${end}.json"
  done
done
