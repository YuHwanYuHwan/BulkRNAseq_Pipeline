#!/bin/bash
# Alignment.sh <group_dir>
#   STAR alignment of trimmed reads -> Processed/<project>/<group>/Alignment_result/<sample>/
#   sjdbOverhang is derived from the data (max read length - 1); the matching index is
#   built on demand by RefIndexing.sh and reused forever.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,5p' "$0"; exit 1; }
init_group "$1"

: "${species:?group.conf must define species (e.g. Homo_sapiens)}"

TRIM="${PROC_DIR}/AdapterTrimming_result"
OUT="${PROC_DIR}/Alignment_result"
mkdir -p "$OUT"
[ -d "$TRIM" ] || { echo "[ERROR] run Trimming.sh first" >&2; exit 1; }

# Overhang from the longest trimmed read in the group. One value per group keeps every
# sample on the same index, so counts stay comparable within the group.
MAXLEN=0
while IFS=$'\t' read -r sample _ _; do
    for f in "${TRIM}/${sample}"_1_trimmed.fastq.gz "${TRIM}/${sample}"_trimmed.fastq.gz; do
        [ -e "$f" ] || continue
        L=$(max_read_length "$f")
        [ "$L" -gt "$MAXLEN" ] && MAXLEN=$L
        break
    done
done < <(list_samples)
[ "$MAXLEN" -gt 0 ] || { echo "[ERROR] could not read trimmed FASTQ" >&2; exit 1; }
OVERHANG=$((MAXLEN - 1))
echo "[LEN ] max trimmed read = ${MAXLEN}bp -> sjdbOverhang=${OVERHANG}"
echo "$OVERHANG" > "${PROC_DIR}/.overhang"   # picked up by MultiQC.sh for the report

bash "$(dirname "${BASH_SOURCE[0]}")/RefIndexing.sh" "$species" "$OVERHANG"
IDX="${REF_ROOT}/${species}/index/overhang${OVERHANG}"
# The index was built from the GTF, so the annotation is already in it. Passing --sjdbGTFfile
# here as well makes STAR re-insert the same junctions for every sample - minutes and extra
# memory each time, for a splice database it already has.

n=0 skip=0
while IFS=$'\t' read -r sample _ _; do
    SDIR="${OUT}/${sample}"
    if is_done "$OUT" "$sample"; then echo "[SKIP] $sample"; skip=$((skip+1)); continue; fi
    # STAR refuses to start if _STARtmp from a killed run is left behind
    rm -rf "${SDIR}/_STARtmp"
    mkdir -p "$SDIR"

    R1="${TRIM}/${sample}_1_trimmed.fastq.gz"
    R2="${TRIM}/${sample}_2_trimmed.fastq.gz"
    [ -e "$R1" ] || R1="${TRIM}/${sample}_trimmed.fastq.gz"
    [ -e "$R2" ] || R2=""

    echo "[STAR] $sample"
    STAR --genomeDir "$IDX" \
        --readFilesIn "$R1" ${R2:+"$R2"} \
        --readFilesCommand zcat \
        --runThreadN "$THREADS" \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNoverLmax 0.03 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "${SDIR}/${sample}"

    BAM="${SDIR}/${sample}Aligned.sortedByCoord.out.bam"
    # A node failure can leave a 0-byte BAM; .done alone would hide it from the next run
    [ -s "$BAM" ] && [ -f "${SDIR}/${sample}Log.final.out" ] \
        || { echo "[ERROR] $sample produced no valid BAM" >&2; exit 1; }

    mark_done "$OUT" "$sample"
    n=$((n+1))
done < <(list_samples)

# The counter reports every sample in the group, not only the ones processed this time -
# a resumed run that skips all of them still has to read as success.
SKIPMSG=""
[ "$skip" -gt 0 ] && SKIPMSG="  ($skip already done)"
echo "[DONE] STAR $((n+skip)) samples -> $OUT$SKIPMSG"
echo "       next: probe_strandedness.sh $GROUP_DIR"
