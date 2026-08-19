#!/bin/bash
# list_samples 규칙 self-check: 평면 파일 = 개별 샘플, 하위 폴더 = merge 대상
set -euo pipefail
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
ROOT="$TMP/repo"; mkdir -p "$ROOT/Scripts/lib" "$ROOT/rawData/P/g"
cp "$(dirname "$0")/../lib/common.sh" "$ROOT/Scripts/lib/"

G="$ROOT/rawData/P/g"
touch "$G/SRR001_1.fastq.gz" "$G/SRR001_2.fastq.gz"   # paired 평면
touch "$G/SRR009.fastq.gz"                            # single 평면
mkdir -p "$G/GSM_B"
touch "$G/GSM_B/SRR002_1.fastq.gz" "$G/GSM_B/SRR002_2.fastq.gz" \
      "$G/GSM_B/SRR003_1.fastq.gz" "$G/GSM_B/SRR003_2.fastq.gz"

source "$ROOT/Scripts/lib/common.sh"
init_group "$G"
OUTPUT=$(list_samples | sort)

grep -q "^GSM_B" <<< "$OUTPUT" || { echo "FAIL: GSM_B 폴더 미인식"; exit 1; }
[ "$(grep -c '^GSM_B' <<< "$OUTPUT")" -eq 1 ] || { echo "FAIL: GSM_B 는 1행이어야 함"; exit 1; }
grep "^GSM_B" <<< "$OUTPUT" | grep -q "SRR002_1.*,.*SRR003_1" || { echo "FAIL: GSM_B R1 2개 병합 안 됨"; exit 1; }
grep -q "^SRR001" <<< "$OUTPUT" || { echo "FAIL: SRR001 미인식"; exit 1; }
grep "^SRR001" <<< "$OUTPUT" | grep -q "SRR001_2.fastq.gz" || { echo "FAIL: SRR001 R2 미연결"; exit 1; }
grep -q "^SRR009" <<< "$OUTPUT" || { echo "FAIL: single-end 미인식"; exit 1; }
[ "$(grep -c . <<< "$OUTPUT")" -eq 3 ] || { echo "FAIL: 샘플 수 = $(grep -c . <<< "$OUTPUT") (기대 3)"; echo "$OUTPUT"; exit 1; }

[ "$PROC_DIR" = "$ROOT/Processed/P/g" ] || { echo "FAIL: PROC_DIR=$PROC_DIR"; exit 1; }
echo "OK: 평면 paired 1 + single 1 + 폴더 merge 1 = 3 samples, 경로 미러링 정상"
