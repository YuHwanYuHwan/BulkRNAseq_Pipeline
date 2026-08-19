#!/bin/bash
# group.conf parsing: spaces around '=' , quotes, comments, blank lines must all work.
set -euo pipefail
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
ROOT="$TMP/repo"; mkdir -p "$ROOT/Scripts/lib" "$ROOT/rawData/P/g"
cp "$(dirname "$0")/../lib/common.sh" "$ROOT/Scripts/lib/"

cat > "$ROOT/rawData/P/g/group.conf" <<'CONF'
# comment line
species       = Homo_sapiens
adapter_kit="Illumina_TruSeq"      # trailing comment
strandedness  =

CONF

source "$ROOT/Scripts/lib/common.sh"
init_group "$ROOT/rawData/P/g"

[ "${species:-}"      = "Homo_sapiens"    ] || { echo "FAIL species='${species:-}'"; exit 1; }
[ "${adapter_kit:-}"  = "Illumina_TruSeq" ] || { echo "FAIL adapter_kit='${adapter_kit:-}'"; exit 1; }
[ -z "${strandedness:-}" ]                  || { echo "FAIL strandedness should be empty"; exit 1; }
echo "OK: spaces / quotes / comments / empty value all parsed"
