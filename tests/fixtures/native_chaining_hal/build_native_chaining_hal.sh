#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
build="${1:-$root/build}"
for tool in maf2hal halStats halAlignedExtract; do
  command -v "$tool" >/dev/null || { echo "missing $tool" >&2; exit 2; }
done
rm -rf "$build"
mkdir -p "$build"

cat > "$build/root.maf" <<'MAF'
##maf version=1 scoring=native-chaining

a score=0
s Root.chrRoot 0 80 + 80 AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT
s Anc1.chrAnc1 0 80 + 80 AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT
s Anc2.chrAnc2 0 80 + 80 AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT
MAF

cat > "$build/anc1.maf" <<'MAF'
##maf version=1 scoring=native-chaining

a score=0
s Anc1.chrAnc1 0 20 + 80 AAAAAAAAAACCCCCCCCCC
s A.chrA 0 20 + 80 AAAAAAAAAACCCCCCCCCC
s B.chrB 0 20 + 80 AAAAAAAAAACCCCCCCCCC

a score=0
s Anc1.chrAnc1 20 20 + 80 GGGGGGGGGGTTTTTTTTTT
s A.chrA 20 20 + 80 GGGGGGGGGGTTTTTTTTTT
s B.chrB 20 20 + 80 GGGGGGGGGGTTTTTTTTTT

a score=0
s Anc1.chrAnc1 40 20 + 80 AAAAAAAAAACCCCCCCCCC
s A.chrA 40 20 - 80 GGGGGGGGGGTTTTTTTTTT
s B.chrB 40 20 + 80 AAAAAAAAAACCCCCCCCCC

a score=0
s Anc1.chrAnc1 60 20 + 80 GGGGGGGGGGTTTTTTTTTT
s A.chrA 60 20 + 80 GGGGGGGGGGTTTTTTTTTT
s B.chrB 60 20 + 80 GGGGGGGGGGTTTTTTTTTT
MAF

cat > "$build/anc2.maf" <<'MAF'
##maf version=1 scoring=native-chaining

a score=0
s Anc2.chrAnc2 0 20 + 80 AAAAAAAAAACCCCCCCCCC
s C.chrC 0 20 + 80 AAAAAAAAAACCCCCCCCCC
s D.chrD 0 20 + 80 AAAAAAAAAACCCCCCCCCC

a score=0
s Anc2.chrAnc2 20 20 + 80 GGGGGGGGGGTTTTTTTTTT
s C.chrC 20 20 + 80 GGGGGGGGGGTTTTTTTTTT
s D.chrD 20 20 + 80 GGGGGGGGGGTTTTTTTTTT

a score=0
s Anc2.chrAnc2 40 20 + 80 AAAAAAAAAACCCCCCCCCC
s C.chrC 40 20 + 80 AAAAAAAAAACCCCCCCCCC
s D.chrD 40 20 + 80 AAAAAAAAAACCCCCCCCCC

a score=0
s Anc2.chrAnc2 60 20 + 80 GGGGGGGGGGTTTTTTTTTT
s C.chrC 60 20 + 80 GGGGGGGGGGTTTTTTTTTT
s D.chrD2 60 20 + 80 GGGGGGGGGGTTTTTTTTTT
MAF

maf2hal "$build/root.maf" "$build/test.hal" --refGenome Root --targetGenomes Anc1,Anc2
maf2hal "$build/anc1.maf" "$build/test.hal" --append --refGenome Anc1 --targetGenomes A,B
maf2hal "$build/anc2.maf" "$build/test.hal" --append --refGenome Anc2 --targetGenomes C,D
halStats "$build/test.hal" --tree > "$build/tree.nwk"
halStats "$build/test.hal" --genomes > "$build/genomes.txt"
echo "wrote $build/test.hal"
