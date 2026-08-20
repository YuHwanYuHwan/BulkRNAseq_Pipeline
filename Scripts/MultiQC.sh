#!/bin/bash
# MultiQC.sh <group_dir>
#   Aggregate QC with MultiQC and write pipeline_report.txt - the versions, parameters
#   and QC numbers you need for a paper Methods section, plus a draft paragraph.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,4p' "$0"; exit 1; }
init_group "$1"

multiqc "$PROC_DIR" -o "$OUT_DIR" -n "${GROUP}_multiqc_report" -f -q

GTF="$(ls "${REF_ROOT}/${species}"/*.gtf | head -1)"
FA="$(ls "${REF_ROOT}/${species}"/*.dna.*.fa 2>/dev/null | head -1)"
OVERHANG="$(cat "${PROC_DIR}/.overhang" 2>/dev/null || echo '?')"
KIT="${adapter_kit:-Illumina_universal}"
read -r A1 A2 < <(awk -F, -v k="$KIT" '$1==k { print $2, $3 }' "$(dirname "${BASH_SOURCE[0]}")/AdapterSequenceList.csv")
SAMPLE_TSV=$(list_samples)
NSAMP=$(wc -l <<< "$SAMPLE_TSV")
LAYOUT=$(awk -F'	' 'NR==1 { print ($3=="" ? "single-end" : "paired-end") }' <<< "$SAMPLE_TSV")
# Ensembl release is read from the GTF filename (Homo_sapiens.GRCh38.113.gtf) so the
# report reflects what is actually on disk, not what someone meant to download.
RELEASE=$(basename "$GTF" | sed -E 's/.*\.([0-9]+)\.gtf/\1/')
# Genome build from the FASTA name, and the species with a space, so the Methods paragraph
# reads as prose rather than as a directory name.
BUILD=$(basename "${FA:-}" | cut -d. -f2)
SPECIES_PROSE="${species//_/ }"
READLEN=$((OVERHANG + 1))
case "$strandedness" in
    no)      STRAND_PROSE="unstranded mode" ;;
    reverse) STRAND_PROSE="reverse-stranded mode" ;;
    yes)     STRAND_PROSE="forward-stranded mode" ;;
esac

# STAR unique mapping rate and sequencing depth across samples
MAPSTAT=$(find "${PROC_DIR}/Alignment_result" -name '*Log.final.out' -exec \
    grep -H 'Uniquely mapped reads %' {} + 2>/dev/null | sed 's/.*|\s*//; s/%//' |
    awk '{ s+=$1; if (min=="" || $1<min) min=$1 } END { if (NR) printf "mean %.1f%%  min %.1f%%", s/NR, min }')
DEPTH=$(find "${PROC_DIR}/Alignment_result" -name '*Log.final.out' -exec \
    grep -H 'Number of input reads' {} + 2>/dev/null | sed 's/.*|\s*//' |
    awk '{ s+=$1; if (min=="" || $1<min) min=$1; if ($1>max) max=$1 }
         END { if (NR) printf "mean %.1fM  range %.1f-%.1fM", s/NR/1e6, min/1e6, max/1e6 }')
# One awk over every counts file; FILENAME keys the per-sample totals.
NFSTAT=$(awk -F'\t' '/^__no_feature/ { nf[FILENAME]=$2 } { t[FILENAME]+=$2 }
    END { for (f in t) if (t[f]>0) { r=100*nf[f]/t[f]; s+=r; c++; if (r>max) max=r }
          if (c) printf "mean %.1f%%  max %.1f%%", s/c, max }' \
    $(find "${PROC_DIR}/HTseqCount_result" -name '*.gene.counts') 2>/dev/null)

# Short forms for the prose paragraph
DEPTH_PROSE=$(sed -E 's/mean ([0-9.]+M).*/an average of \1/' <<< "${DEPTH:-n\/a}")
MAP_PROSE=$(sed -E 's/mean ([0-9.]+%).*/\1/' <<< "${MAPSTAT:-n\/a}")
STAR_VER=$(tool_version STAR | sed 's/^STAR /v/')

REPORT="${OUT_DIR}/${GROUP}_pipeline_report.txt"
cat > "$REPORT" <<TXT
================================================================
 Bulk RNA-seq Preprocessing Report
 group   : ${GROUP}
 run     : $(date '+%F %H:%M')
 samples : ${NSAMP}  (${LAYOUT}, ${READLEN} bp)
================================================================

[ Reference ]
  species       ${SPECIES_PROSE}  (${BUILD:-?})
  genome FASTA  $(basename "${FA:-n/a}")
  annotation    $(basename "$GTF")   (Ensembl release ${RELEASE})
  STAR index    overhang ${OVERHANG}

[ Software ]
  $(tool_version fastqc)
  $(tool_version cutadapt)
  $(tool_version STAR)
  $(tool_version htseq)
  $(tool_version multiqc)
  $(tool_version R)

[ Parameters ]
  adapter kit   ${KIT}
    R1          ${A1}
    R2          ${A2:-(single-end)}
  cutadapt      --minimum-length 20 --pair-filter=any
  STAR          --outFilterMismatchNoverLmax 0.03 --outFilterMultimapNmax 10
                --outSAMtype BAM SortedByCoordinate
                (annotation and sjdbOverhang ${OVERHANG} come from the index)
  HTSeq         -r pos -s ${strandedness}
  normalization edgeR TMM -> CPM

[ QC summary ]
  input reads       ${DEPTH:-n/a}   per sample, after trimming
  uniquely mapped   ${MAPSTAT:-n/a}
  __no_feature      ${NFSTAT:-n/a}
  MultiQC           ${GROUP}_multiqc_report.html

[ Methods draft ]
  Sequencing was ${LAYOUT} at ${READLEN} bp, with ${DEPTH_PROSE} reads per sample. Raw
  reads were assessed with FastQC and adapters were removed with cutadapt (minimum
  length 20 bp). Trimmed reads were aligned to the ${SPECIES_PROSE} ${BUILD:-} reference
  genome with STAR ${STAR_VER} (sjdbOverhang ${OVERHANG}, mismatch rate <= 0.03, up to 10
  multimapping loci) using Ensembl release ${RELEASE} annotation; ${MAP_PROSE} of reads
  mapped uniquely. Gene-level counts were obtained with HTSeq-count in ${STRAND_PROSE} and
  normalized to counts per million using the TMM method in edgeR.
TXT

echo "[DONE] $REPORT"
