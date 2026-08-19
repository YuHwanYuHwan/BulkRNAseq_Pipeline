#!/bin/bash
# max_read_length self-check: variable-length reads must yield the MAX, not mean/first.
set -euo pipefail
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
ROOT="$TMP/repo"; mkdir -p "$ROOT/Scripts/lib"
cp "$(dirname "$0")/../lib/common.sh" "$ROOT/Scripts/lib/"
source "$ROOT/Scripts/lib/common.sh"

FQ="$TMP/t.fastq"
{ printf '@r1\n%s\n+\n%s\n' "$(printf 'A%.0s' {1..100})" "$(printf 'I%.0s' {1..100})"
  printf '@r2\n%s\n+\n%s\n' "$(printf 'A%.0s' {1..150})" "$(printf 'I%.0s' {1..150})"
  printf '@r3\n%s\n+\n%s\n' "$(printf 'A%.0s' {1..75})"  "$(printf 'I%.0s' {1..75})"; } > "$FQ"

L=$(max_read_length "$FQ")
[ "$L" -eq 150 ] || { echo "FAIL: max_read_length=$L (expected 150)"; exit 1; }

gzip -c "$FQ" > "$FQ.gz"
LZ=$(max_read_length "$FQ.gz")
[ "$LZ" -eq 150 ] || { echo "FAIL: gzip max_read_length=$LZ (expected 150)"; exit 1; }

echo "OK: max read 150 (not 100/75/108avg) -> sjdbOverhang 149, gzip handled"
