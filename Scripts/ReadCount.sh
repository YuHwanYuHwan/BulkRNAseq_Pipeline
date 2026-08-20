#!/bin/bash
# ReadCount.sh <group_dir>
#   HTSeq gene counts per sample, then a merged count matrix for the group.
#   This is the first step of stage 2: it refuses to run until strandedness is set,
#   because a wrong value silently halves the counts instead of raising an error.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,5p' "$0"; exit 1; }
init_group "$1"

: "${species:?group.conf must define species}"
case "${strandedness:-}" in
    no|yes|reverse) ;;
    *) echo "[ERROR] group.conf strandedness must be no|yes|reverse (got '${strandedness:-empty}')" >&2
       echo "        run: bash Scripts/probe_strandedness.sh ${GROUP_DIR}" >&2; exit 1 ;;
esac

ALIGN="${PROC_DIR}/Alignment_result"
OUT="${PROC_DIR}/HTseqCount_result"
mkdir -p "$OUT"
GTF="$(ls "${REF_ROOT}/${species}"/*.gtf | head -1)"
echo "[HTSEQ] strandedness=$strandedness  gtf=$(basename "$GTF")"

n=0 skip=0
SAMPLES=()
while IFS=$'\t' read -r sample _ _; do
    SAMPLES+=("$sample")
    CNT="${OUT}/${sample}.gene.counts"
    if is_done "$OUT" "$sample"; then echo "[SKIP] $sample"; skip=$((skip+1)); continue; fi
    BAM="${ALIGN}/${sample}/${sample}Aligned.sortedByCoord.out.bam"
    [ -s "$BAM" ] || { echo "[ERROR] missing BAM for $sample" >&2; exit 1; }

    echo "[CNT ] $sample"
    htseq-count -r pos -s "$strandedness" "$BAM" "$GTF" > "$CNT"

    # Wrong strandedness produces no error, only a collapsed assignment rate. Warn loudly.
    awk -F'\t' -v s="$sample" '
        /^__no_feature/ { nf=$2 } { t+=$2 }
        END { if (t>0 && 100*nf/t > 50)
                 printf "[WARN] %s: __no_feature %.1f%% - strandedness may be wrong\n", s, 100*nf/t }' "$CNT"

    mark_done "$OUT" "$sample"
    n=$((n+1))
done < <(list_samples)

# ── merge into one matrix ────────────────────────────────────────────────────
MATRIX="${OUT_DIR}/${GROUP}_count_matrix.tsv"

# Header starts with a leading "Gene" and uses leading tabs, never a trailing one -
# a trailing tab makes R read.table() add an empty column and turn everything character.
{ printf 'Gene'; printf '\t%s' "${SAMPLES[@]}"; printf '\n'; } > "$MATRIX"

# Gene column from the first sample, then one count column per sample.
CMD="paste <(grep -v '^__' ${OUT}/${SAMPLES[0]}.gene.counts | cut -f1)"
for s in "${SAMPLES[@]}"; do
    CMD="$CMD <(grep -v '^__' ${OUT}/${s}.gene.counts | cut -f2)"
done
eval "$CMD" >> "$MATRIX"

# The counter reports every sample in the group, not only the ones processed this time -
# a resumed run that skips all of them still has to read as success.
SKIPMSG=""
[ "$skip" -gt 0 ] && SKIPMSG="  ($skip already done)"
echo "[DONE] HTSeq $((n+skip)) samples$SKIPMSG"
echo "       matrix: $MATRIX  ($(( $(wc -l < "$MATRIX") - 1 )) genes x ${#SAMPLES[@]} samples)"
