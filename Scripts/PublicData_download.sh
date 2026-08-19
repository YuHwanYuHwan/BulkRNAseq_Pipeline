#!/bin/bash
# PublicData_download.sh <group_dir> <SRR ...>
#   SRR 목록을 받아 FASTQ + runinfo 메타데이터를 group 폴더에 내려받는다.
#   같은 샘플(SampleName)에 run 이 여러 개면 하위 폴더로 묶어 merge 대상임을 표현한다.
#
#   bash PublicData_download.sh rawData/MSK_Map/sarcopenia SRR27743399 SRR27743400
#   bash PublicData_download.sh rawData/MSK_Map/sarcopenia srr_list.txt
set -euo pipefail

RUNINFO_URL="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc="

[ $# -ge 2 ] || { sed -n '2,7p' "$0"; exit 1; }

GROUP_DIR="$1"; shift
# 인자가 파일 하나면 그 안의 SRR 목록을 쓴다
if [ $# -eq 1 ] && [ -f "$1" ]; then
    mapfile -t ACCS < <(grep -oE 'SRR[0-9]+' "$1")
else
    ACCS=("$@")
fi
[ ${#ACCS[@]} -gt 0 ] || { echo "[ERROR] SRR 없음"; exit 1; }

mkdir -p "$GROUP_DIR"
META="${GROUP_DIR}/metadata.csv"

# ── 1. 메타데이터 ────────────────────────────────────────────────────────────
curl -sf "${RUNINFO_URL}${ACCS[0]}" > "$META"
for acc in "${ACCS[@]:1}"; do
    curl -sf "${RUNINFO_URL}${acc}" | tail -n +2 >> "$META"
done
[ -s "$META" ] || { echo "[ERROR] runinfo 조회 실패"; exit 1; }
echo "[META] $META  ($(( $(wc -l < "$META") - 1 )) runs)"

# ── 2. acc -> SampleName map, one pass. Samples with >1 run go into a subfolder.
declare -A SAMPLE_OF MULTI count
while IFS=$'	' read -r acc smp; do
    SAMPLE_OF[$acc]="$smp"
    n=$(( ${count[$smp]:-0} + 1 )); count[$smp]=$n
    if [ "$n" -gt 1 ]; then MULTI[$smp]=1; fi
done < <(awk -F, '
    NR==1 { for (i=1;i<=NF;i++) if ($i=="SampleName") c=i; next }
    c { print $1 "	" $c }' "$META")

# ── 3. 다운로드 ──────────────────────────────────────────────────────────────
for acc in "${ACCS[@]}"; do
    sample="${SAMPLE_OF[$acc]:-$acc}"
    dest="${GROUP_DIR}"
    [ -n "${MULTI[$sample]:-}" ] && dest="${GROUP_DIR}/${sample}"   # multi-run -> merge folder
    mkdir -p "$dest"

    if [ -f "${dest}/.${acc}.done" ]; then
        echo "[SKIP] $acc"
        continue
    fi
    echo "[GET ] $acc -> $dest"
    prefetch --output-directory "$dest" "$acc" >/dev/null
    fasterq-dump --split-3 --outdir "$dest" "${dest}/${acc}" >/dev/null
    rm -rf "${dest:?}/${acc}"
    gzip -f "${dest}"/${acc}*.fastq
    touch "${dest}/.${acc}.done"
done

echo "[DONE] ${#ACCS[@]} runs -> $GROUP_DIR"
echo "       다음: group.conf 에 species 를 적고 run_stage1.sh 실행"
