#!/bin/bash
# list_samples.sh <group_dir>
#   Print what the pipeline sees in a group: sample <TAB> R1 <TAB> R2.
#   Worth running before a job that takes hours - a sample missing here stays missing.
NO_STEP_LOG=1
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,4p' "$0"; exit 1; }
init_group "$1"
list_samples
