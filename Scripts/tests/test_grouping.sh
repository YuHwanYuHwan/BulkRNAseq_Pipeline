#!/bin/bash
# acc -> SampleName mapping: runs sharing a SampleName must be flagged for a merge folder.
set -euo pipefail
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
META="$TMP/metadata.csv"
cat > "$META" <<'CSV'
Run,LibraryLayout,SampleName,ScientificName
SRR001,PAIRED,GSM_A,Homo sapiens
SRR002,PAIRED,GSM_B,Homo sapiens
SRR003,PAIRED,GSM_B,Homo sapiens
CSV

declare -A SAMPLE_OF MULTI count
while IFS=$'	' read -r acc smp; do
    SAMPLE_OF[$acc]="$smp"
    n=$(( ${count[$smp]:-0} + 1 )); count[$smp]=$n
    if [ "$n" -gt 1 ]; then MULTI[$smp]=1; fi
done < <(awk -F, '
    NR==1 { for (i=1;i<=NF;i++) if ($i=="SampleName") c=i; next }
    c { print $1 "\t" $c }' "$META")

[ "${SAMPLE_OF[SRR003]}" = "GSM_B" ]  || { echo "FAIL: SRR003 -> ${SAMPLE_OF[SRR003]:-}"; exit 1; }
[ -n "${MULTI[GSM_B]:-}" ]            || { echo "FAIL: GSM_B not flagged multi-run"; exit 1; }
[ -z "${MULTI[GSM_A]:-}" ]            || { echo "FAIL: GSM_A wrongly flagged multi-run"; exit 1; }
[ -z "${SAMPLE_OF[SRR999]:-}" ]       || { echo "FAIL: unknown acc must be empty"; exit 1; }
echo "OK: GSM_B(2 runs) -> merge folder, GSM_A(1 run) -> flat"
