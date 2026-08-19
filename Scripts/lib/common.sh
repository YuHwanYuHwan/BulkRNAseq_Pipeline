# common.sh - shared helpers sourced by every step script.
# Contract: each step takes ONE group directory. Directory layout is the source of truth.
#   subdirectory present -> directory name = sample_id, FASTQs inside are merged
#   flat FASTQ files     -> filename stem (minus _1/_2) = sample_id

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
[ -f "${PIPELINE_ROOT}/config.sh" ] && source "${PIPELINE_ROOT}/config.sh"
REF_ROOT="${PIPELINE_ROOT}/reference_Genomes"

# Tool paths: use PATH if available, else *_BIN from config.sh
if [ -n "${CONDA_ENV:-}" ] && ! command -v fastqc >/dev/null 2>&1 && command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate "$CONDA_ENV"
fi
FASTQC_BIN="${FASTQC_BIN:-fastqc}"

# Cores. Inside a SLURM job SLURM_CPUS_PER_TASK is what was actually reserved; outside one,
# nproc reports the CPUs this process may use. A fixed default would either waste a big
# compute node or oversubscribe a shared login node.
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}"

# group dir -> project/group relative path. Processed/ and Output/ mirror the same layout.
init_group() {
    GROUP_DIR="$(cd "$1" && pwd)"
    REL="${GROUP_DIR#*/rawData/}"
    [ "$REL" != "$GROUP_DIR" ] || { echo "[ERROR] not under rawData/: $1" >&2; exit 1; }
    PROC_DIR="${PIPELINE_ROOT}/Processed/${REL}"
    OUT_DIR="${PIPELINE_ROOT}/Output/${REL}"
    GROUP="$(basename "$REL")"
    # group.conf is written by hand, so tolerate spaces around '=' and quotes.
    # Sourcing it directly would fail on "species = Homo_sapiens".
    CONF="${GROUP_DIR}/group.conf"
    if [ -f "$CONF" ]; then
        while IFS= read -r line; do
            line="${line%%#*}"
            [[ "$line" == *=* ]] || continue
            local k="${line%%=*}" v="${line#*=}"
            k="$(echo "$k" | tr -d '[:space:]')"
            v="$(echo "$v" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//; s/^"\(.*\)"$/\1/; s/^'"'"'\(.*\)'"'"'$/\1/')"
            [ -n "$k" ] && printf -v "$k" '%s' "$v"
        done < "$CONF"
    fi
    mkdir -p "$PROC_DIR" "$OUT_DIR"
}

# Sample list. One line per sample: sample_id<TAB>R1(comma-joined)<TAB>R2(comma-joined, empty if single-end)
list_samples() {
    local d f f2 stem s
    for d in "$GROUP_DIR"/*/; do
        [ -d "$d" ] || continue
        s="$(basename "$d")"
        printf '%s\t%s\t%s\n' "$s" \
            "$(ls "$d"*_1.fastq* "$d"*_1.fq* 2>/dev/null | paste -sd,)" \
            "$(ls "$d"*_2.fastq* "$d"*_2.fq* 2>/dev/null | paste -sd,)"
    done
    for f in "$GROUP_DIR"/*.fastq.gz "$GROUP_DIR"/*.fastq "$GROUP_DIR"/*.fq.gz "$GROUP_DIR"/*.fq; do
        [ -e "$f" ] || continue
        stem="$(basename "$f")"; stem="${stem%%.*}"
        case "$stem" in
            *_2) continue ;;                  # R2 is handled together with its R1
            # Only a _1 stem can have a mate. Without this branch the substitution below is a
            # no-op and ls returns the file itself, pairing a single-end read with itself.
            *_1) printf '%s\t%s\t%s\n' "${stem%_1}" "$f" "$(ls "${f/_1./_2.}" 2>/dev/null || true)"
                 continue ;;
        esac
        printf '%s\t%s\t\n' "$stem" "$f"
    done
}

is_done()   { [ -f "${1}/.${2}.done" ]; }
mark_done() { touch "${1}/.${2}.done"; }

# Max read length over the first N reads of a FASTQ (line 2 of every 4).
# Trimming makes read length variable, so the MAX must be used - sjdbOverhang has to
# accommodate the longest read that can span a junction.
max_read_length() {   # $1=fastq(.gz)  [$2=reads to scan, default 10000]
    local f="$1" n="${2:-10000}"
    zcat -f "$f" 2>/dev/null | awk -v n="$n" 'NR%4==2 { if (length($0)>m) m=length($0); c++ }
                                              c>=n { exit } END { print m+0 }'
}

# Tool versions, one line each -> material for pipeline_report.txt
tool_version() {
    case "$1" in
        fastqc)   "$FASTQC_BIN" --version 2>&1 | head -1 ;;
        cutadapt) echo "cutadapt $(cutadapt --version 2>&1 | head -1)" ;;
        STAR)     STAR --version 2>&1 | head -1 | sed 's/^/STAR /' ;;
        htseq)    htseq-count --version 2>&1 | head -1 | sed 's/^/HTSeq /' ;;
        multiqc)  multiqc --version 2>&1 | head -1 ;;
        R)        R --version 2>&1 | head -1 ;;
    esac
}
