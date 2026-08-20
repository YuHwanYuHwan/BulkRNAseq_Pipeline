#!/bin/bash
# Trimming.sh <group_dir>
#   Adapter removal with cutadapt -> Processed/<project>/<group>/AdapterTrimming_result/
#   The adapter is looked up in AdapterSequenceList.csv by the adapter_kit set in group.conf.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,4p' "$0"; exit 1; }
init_group "$1"

KIT="${adapter_kit:-Illumina_universal}"
LIST="$(dirname "${BASH_SOURCE[0]}")/AdapterSequenceList.csv"
read -r A1 A2 < <(awk -F, -v k="$KIT" '$1==k { print $2, $3 }' "$LIST")
[ -n "${A1:-}" ] || { echo "[ERROR] kit '$KIT' not in $LIST - add a row or fix group.conf" >&2; exit 1; }
echo "[KIT ] $KIT  R1=$A1  R2=${A2:-(none)}"

OUT="${PROC_DIR}/AdapterTrimming_result"
mkdir -p "$OUT"

n=0 skip=0
while IFS=$'	' read -r sample r1 r2; do
    if is_done "$OUT" "$sample"; then echo "[SKIP] $sample"; skip=$((skip+1)); continue; fi
    echo "[TRIM] $sample"
    # A sample split over several runs is concatenated first - joined gzip stays valid gzip.
    M1="${OUT}/${sample}_1.merged.fastq.gz"; in1="$r1"
    [[ "$r1" == *,* ]] && { cat ${r1//,/ } > "$M1"; in1="$M1"; }
    if [ -n "$r2" ]; then
        M2="${OUT}/${sample}_2.merged.fastq.gz"; in2="$r2"
        [[ "$r2" == *,* ]] && { cat ${r2//,/ } > "$M2"; in2="$M2"; }
        cutadapt -j "$THREADS" -a "$A1" -A "$A2" --pair-filter=any --minimum-length 20 \
            -o "${OUT}/${sample}_1_trimmed.fastq.gz" -p "${OUT}/${sample}_2_trimmed.fastq.gz" \
            "$in1" "$in2" > "${OUT}/${sample}.cutadapt.log"
    else
        cutadapt -j "$THREADS" -a "$A1" --minimum-length 20 \
            -o "${OUT}/${sample}_trimmed.fastq.gz" \
            "$in1" > "${OUT}/${sample}.cutadapt.log"
    fi
    rm -f "${OUT}/${sample}"_[12].merged.fastq.gz
    mark_done "$OUT" "$sample"
    n=$((n+1))
done < <(list_samples)

# The counter reports every sample in the group, not only the ones processed this time -
# a resumed run that skips all of them still has to read as success.
SKIPMSG=""
[ "$skip" -gt 0 ] && SKIPMSG="  ($skip already done)"
echo "[DONE] cutadapt $((n+skip)) samples -> $OUT$SKIPMSG"
