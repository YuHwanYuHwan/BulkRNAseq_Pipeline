#!/bin/bash
# FastQC.sh <group_dir>
#   FastQC on every raw FASTQ in the group -> Processed/<project>/<group>/Fastqc_result/
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,3p' "$0"; exit 1; }
init_group "$1"

OUT="${PROC_DIR}/Fastqc_result"
mkdir -p "$OUT"

n=0
while IFS=$'\t' read -r sample r1 r2; do
    if is_done "$OUT" "$sample"; then
        echo "[SKIP] $sample"
        continue
    fi
    echo "[FQC ] $sample"
    "$FASTQC_BIN" --threads "$THREADS" --outdir "$OUT" ${r1//,/ } ${r2//,/ }
    mark_done "$OUT" "$sample"
    n=$((n+1))
done < <(list_samples)

echo "[DONE] FastQC $n samples -> $OUT"
