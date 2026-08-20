#!/bin/bash
# run_stage2.sh <group_dir>
#   HTSeq -> count matrix -> CPM -> MultiQC + pipeline_report.txt
#   Refuses to start unless group.conf has strandedness (ReadCount.sh enforces it).
#SBATCH --job-name=rnaseq_s2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
# HTSeq counts one BAM at a time in a single process, so extra cores buy nothing here.
#SBATCH --output=logs/stage2_%j.out
set -euo pipefail
[ $# -eq 1 ] || { sed -n "2,4p" "$0"; exit 1; }
# sbatch copies this file to a spool directory, so its own path says nothing about where the
# repository is. The group directory does: it always lives under <repo>/rawData/.
GROUP_ABS="$(cd "$1" && pwd)"
S="${GROUP_ABS%%/rawData/*}/Scripts"
[ -d "$S" ] || { echo "[ERROR] $1 is not under a pipeline rawData/ directory" >&2; exit 1; }

bash "$S/ReadCount.sh" "$1"
bash "$S/CalcCPM.sh"   "$1"
bash "$S/MultiQC.sh"   "$1"
