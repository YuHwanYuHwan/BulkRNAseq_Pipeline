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
[ $# -eq 1 ] || { sed -n "2,5p" "$0"; exit 1; }
# sbatch copies this file to a spool directory, so its own path says nothing about where the
# repository is. The group directory does: it always lives under <repo>/rawData/.
GROUP_ABS="$(cd "$1" && pwd)"
S="${GROUP_ABS%%/rawData/*}/Scripts"
[ -d "$S" ] || { echo "[ERROR] $1 is not under a pipeline rawData/ directory" >&2; exit 1; }

bash "$S/FastQC.sh"    "$1"
bash "$S/Trimming.sh"  "$1"
bash "$S/Alignment.sh" "$1"
bash "$S/probe_strandedness.sh" "$1"
