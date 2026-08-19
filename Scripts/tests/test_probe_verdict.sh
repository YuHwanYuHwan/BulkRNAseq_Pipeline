#!/bin/bash
# The verdict thresholds must map __no_feature fraction to the right library type.
set -euo pipefail
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

verdict() {   # $1 = __no_feature count, total 1000
    printf 'GENE1\t%d\n__no_feature\t%d\n__ambiguous\t0\n' "$((1000-$1))" "$1" > "$TMP/c"
    awk -F'\t' '/^__/ { s[$1]=$2; t+=$2; next } { a+=$2; t+=$2 }
        END { r=100*(s["__no_feature"]+0)/t
              if (r<35) print "reverse"; else if (r<65) print "no"; else print "yes" }' "$TMP/c"
}

[ "$(verdict 150)" = "reverse" ] || { echo "FAIL: 15% -> $(verdict 150)"; exit 1; }
[ "$(verdict 500)" = "no" ]      || { echo "FAIL: 50% -> $(verdict 500)"; exit 1; }
[ "$(verdict 850)" = "yes" ]     || { echo "FAIL: 85% -> $(verdict 850)"; exit 1; }
echo "OK: 15%->reverse, 50%->no, 85%->yes"
