#!/bin/bash
# setup.sh [--create-env]
#   First-time setup and environment check. Safe to re-run - it only reports what is
#   missing unless --create-env is given.
#
#   bash setup.sh                # check only
#   bash setup.sh --create-env   # also create the conda env if it does not exist
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

[ "${1:-}" = "--create-env" ] && CREATE=1 || CREATE=0
FAIL=0
ok()   { printf '  [ OK ] %s\n' "$1"; }
miss() { printf '  [MISS] %s\n' "$1"; FAIL=1; }
warn() { printf '  [WARN] %s\n' "$1"; }

echo "== config =="
if [ -f config.sh ]; then ok "config.sh"
else
    # Template lives here, not in a separate .example file - one place to keep in sync.
    cat > config.sh <<'CFG'
# config.sh - per-machine settings. Gitignored.

# Conda environment holding the tools. `setup.sh --create-env` creates it, every script
# activates it. Set it empty if the tools are already on PATH some other way.
CONDA_ENV="rnaseq-preproc"

# Only if fastqc is not on PATH (e.g. a manually downloaded copy)
# FASTQC_BIN="/home/user/FastQC/fastqc"
CFG
    ok "config.sh created - edit it if tools are not on PATH"
fi
source config.sh 2>/dev/null || true

echo "== conda env =="
if [ -n "${CONDA_ENV:-}" ] && ! conda env list 2>/dev/null | grep -qE "^${CONDA_ENV}\s"; then
    if [ "$CREATE" = 1 ]; then
        echo "  creating '$CONDA_ENV' ..."
        # Channel order matters: bioconda requires conda-forge to take priority.
        conda create -y -n "$CONDA_ENV" -c conda-forge -c bioconda \
            sra-tools fastqc cutadapt star htseq multiqc bioconductor-edger pigz
    else
        miss "conda env '$CONDA_ENV' not found - rerun with --create-env"
    fi
fi
if [ -n "${CONDA_ENV:-}" ] && command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null && conda activate "$CONDA_ENV"
fi

echo "== tools =="
for t in prefetch fasterq-dump cutadapt STAR htseq-count multiqc Rscript "${FASTQC_BIN:-fastqc}"; do
    command -v "$t" >/dev/null 2>&1 && ok "$(basename "$t")" \
        || miss "$(basename "$t")$([ "$t" = "${FASTQC_BIN:-fastqc}" ] && echo '  (set FASTQC_BIN in config.sh)')"
done
Rscript -e 'quit(status = !requireNamespace("edgeR", quietly=TRUE))' 2>/dev/null && ok "R edgeR" || miss "R package edgeR"

echo "== reference genomes =="
REF_ROOT=reference_Genomes
shopt -s nullglob
found=0
for d in "$REF_ROOT"/*/; do
    sp="$(basename "$d")"
    fa=("$d"*.dna.*.fa); gtf=("$d"*.gtf)
    if [ ${#fa[@]} -gt 0 ] && [ ${#gtf[@]} -gt 0 ]; then
        ok "$sp  ($(basename "${gtf[0]}"))"; found=1
    else
        miss "$sp  missing $([ ${#fa[@]} -eq 0 ] && echo 'genome FASTA') $([ ${#gtf[@]} -eq 0 ] && echo 'GTF')"
    fi
done
[ "$found" = 1 ] || miss "no genome under $REF_ROOT - see README section 1"

echo "== disk =="
AVAIL=$(df -BG --output=avail . 2>/dev/null | tail -1 | tr -dc '0-9')
[ -n "$AVAIL" ] && { [ "$AVAIL" -ge 100 ] && ok "${AVAIL}G free" \
    || warn "${AVAIL}G free - a STAR index alone needs ~30G"; }

echo "== self-check =="
if out=$(bash Scripts/lib/selfcheck.sh 2>&1); then ok "logic $out"; else miss "logic - $out"; fi

echo
if [ "$FAIL" = 0 ]; then
    echo "Ready. Next: put FASTQ under rawData/<project>/<group>/ and write group.conf (species)."
else
    echo "Some items are missing - fix them before running the pipeline."
    exit 1
fi
