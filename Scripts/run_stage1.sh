#!/bin/bash
# run_stage1.sh <group_dir>
#   FastQC -> Trimming -> Alignment -> strandedness probe.
#   Stops at the probe: you decide the strandedness, then run stage 2.
#   Works both as `bash run_stage1.sh <dir>` and `sbatch run_stage1.sh <dir>`.
#SBATCH --job-name=rnaseq_s1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --output=logs/stage1_%j.out
set -euo pipefail
S="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ $# -eq 1 ] || { sed -n '2,5p' "$0"; exit 1; }

bash "$S/FastQC.sh"    "$1"
bash "$S/Trimming.sh"  "$1"
bash "$S/Alignment.sh" "$1"
bash "$S/probe_strandedness.sh" "$1"
