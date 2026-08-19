#!/bin/bash
# probe_strandedness.sh <group_dir>
#   Run HTSeq with -s reverse on ONE sample and report __no_feature.
#   Nothing is decided automatically - you read the table and write the value into group.conf.
#
#   A single reverse run separates all three cases:
#     __no_feature low   (~10-20%)  -> reverse
#     __no_feature ~50%             -> no        (unstranded: half the reads sit on the other strand)
#     __no_feature very high (80%+) -> yes       (forward: nearly nothing is assigned)
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,9p' "$0"; exit 1; }
init_group "$1"

: "${species:?group.conf must define species}"
ALIGN="${PROC_DIR}/Alignment_result"
GTF="$(ls "${REF_ROOT}/${species}"/*.gtf 2>/dev/null | head -1)"
[ -n "$GTF" ] || { echo "[ERROR] no GTF for $species" >&2; exit 1; }

read -r SAMPLE _ _ < <(list_samples)
BAM="${ALIGN}/${SAMPLE}/${SAMPLE}Aligned.sortedByCoord.out.bam"
[ -s "$BAM" ] || { echo "[ERROR] no BAM for $SAMPLE - run Alignment.sh first" >&2; exit 1; }

echo "[PROBE] sample=$SAMPLE  -s reverse"
TMP="${PROC_DIR}/.strand_probe.counts"
htseq-count -r pos -s reverse "$BAM" "$GTF" > "$TMP"

# HTSeq writes __no_feature / __ambiguous / ... as the last lines; everything else is a gene
awk -F'\t' '
    /^__/ { special[$1]=$2; tot+=$2; next }
    { assigned+=$2; tot+=$2 }
    END {
        nf = special["__no_feature"]+0
        printf "\n  total reads counted : %d\n", tot
        printf "  assigned to genes   : %d (%.1f%%)\n", assigned, 100*assigned/tot
        printf "  __no_feature        : %d (%.1f%%)\n", nf, 100*nf/tot
        printf "  __ambiguous         : %d\n\n", special["__ambiguous"]+0
        r = 100*nf/tot
        if (r < 35)      verdict = "reverse   (most reads assigned)"
        else if (r < 65) verdict = "no        (about half assigned -> unstranded)"
        else             verdict = "yes       (almost nothing assigned -> forward)"
        printf "  --> likely strandedness : %s\n", verdict
    }' "$TMP"
rm -f "$TMP"

cat <<MSG

  Thresholds are heuristic. Borderline cases are yours to judge.
  Write the value into ${GROUP_DIR}/group.conf :

      strandedness = no | yes | reverse

  Then run: bash Scripts/ReadCount.sh ${GROUP_DIR}
MSG
