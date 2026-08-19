#!/bin/bash
# Matrix assembly: no trailing tab in header, __ rows excluded, columns aligned to samples.
set -euo pipefail
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
OUT="$TMP/counts"; mkdir -p "$OUT"
SAMPLES=(S1 S2)
printf 'G1\t10\nG2\t20\n__no_feature\t5\n' > "$OUT/S1.gene.counts"
printf 'G1\t30\nG2\t40\n__no_feature\t7\n' > "$OUT/S2.gene.counts"

MATRIX="$TMP/m.tsv"
{ printf 'Gene'; for s in "${SAMPLES[@]}"; do printf '\t%s' "$s"; done; printf '\n'; } > "$MATRIX"
CMD="paste <(grep -v '^__' ${OUT}/${SAMPLES[0]}.gene.counts | cut -f1)"
for s in "${SAMPLES[@]}"; do CMD="$CMD <(grep -v '^__' ${OUT}/${s}.gene.counts | cut -f2)"; done
eval "$CMD" >> "$MATRIX"

head -1 "$MATRIX" | grep -q $'\t$' && { echo "FAIL: trailing tab in header"; exit 1; }
[ "$(head -1 "$MATRIX")" = "$(printf 'Gene\tS1\tS2')" ] || { echo "FAIL: header"; exit 1; }
grep -q '__no_feature' "$MATRIX" && { echo "FAIL: __ row leaked into matrix"; exit 1; }
[ "$(wc -l < "$MATRIX")" -eq 3 ] || { echo "FAIL: expected 1 header + 2 genes"; exit 1; }
[ "$(sed -n '2p' "$MATRIX")" = "$(printf 'G1\t10\t30')" ] || { echo "FAIL: row G1"; exit 1; }
echo "OK: header without trailing tab, __ rows dropped, 2 genes x 2 samples"
