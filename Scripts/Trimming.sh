#!/bin/bash
# Trimming.sh <group_dir>
#   cutadapt 로 어댑터 제거 → Processed/<project>/<group>/AdapterTrimming_result/
#   어댑터는 group.conf 의 adapter_kit 으로 AdapterSequenceList.csv 에서 조회한다.
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

# run 이 여러 개면 cat 으로 합친다 (gzip 은 이어붙여도 유효한 스트림)
merge_runs() {   # $1=콤마목록 $2=출력경로 → 실제 사용할 파일 경로를 stdout 으로
    local list="$1" out="$2"
    if [[ "$list" == *,* ]]; then
        [ -s "$out" ] || cat ${list//,/ } > "$out"
        echo "$out"
    else
        echo "$list"
    fi
}

n=0
while IFS=$'	' read -r sample r1 r2; do
    if is_done "$OUT" "$sample"; then echo "[SKIP] $sample"; continue; fi
    echo "[TRIM] $sample"
    in1=$(merge_runs "$r1" "${OUT}/${sample}_1.merged.fastq.gz")
    if [ -n "$r2" ]; then
        in2=$(merge_runs "$r2" "${OUT}/${sample}_2.merged.fastq.gz")
        cutadapt -a "$A1" -A "$A2" --pair-filter=any --minimum-length 20             -o "${OUT}/${sample}_1_trimmed.fastq.gz" -p "${OUT}/${sample}_2_trimmed.fastq.gz"             "$in1" "$in2" > "${OUT}/${sample}.cutadapt.log"
    else
        cutadapt -a "$A1" --minimum-length 20             -o "${OUT}/${sample}_trimmed.fastq.gz"             "$in1" > "${OUT}/${sample}.cutadapt.log"
    fi
    rm -f "${OUT}/${sample}"_[12].merged.fastq.gz
    mark_done "$OUT" "$sample"
    n=$((n+1))
done < <(list_samples)

echo "[DONE] cutadapt $n samples -> $OUT"
