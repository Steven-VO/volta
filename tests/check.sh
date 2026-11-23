#!/bin/sh
# tests/check.sh -- run volta on parameter sets with known answers
#
# Copyright (C) 2026 Steven Van Overberghe
# This file is part of volta, distributed under the GNU General Public License,
# version 3 or later; see the file LICENSE.
#
# usage: tests/check.sh [-j JOBS]
#
# Compiles the binaries of tests/known_counts.txt (JOBS at a time, default: one
# per CPU; they are cached like any other volta binary, so a second run is
# fast), runs every case and compares the number of graphs with the known one.
# Then a few checks of the output itself: graph6 lines, --output, --part, -v.
# Exit status 0 iff everything passed.

set -u
here=$(cd "$(dirname "$0")" && pwd)
VOLTA=${VOLTA:-$here/../volta}
CASES=$here/known_counts.txt

JOBS=
case ${1:-} in
    -j) JOBS=${2:-} ;;
    -j*) JOBS=${1#-j} ;;
    '') ;;
    *) echo "usage: $0 [-j JOBS]" >&2; exit 2 ;;
esac
if [ -z "$JOBS" ]; then
    JOBS=$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 2)
fi

TMP=$(mktemp -d "${TMPDIR:-/tmp}/volta-check.XXXXXX") || exit 2
trap 'rm -rf "$TMP"' EXIT INT TERM

# The cases, without comments: "EXPECTED ARGS...".
sed -e 's/#.*//' -e 's/[[:space:]]*$//' "$CASES" | grep -v '^$' > "$TMP/cases"
total=$(wc -l < "$TMP/cases" | tr -d ' ')

# 1. Compile every parameter set once, in parallel (-B: build only).
echo "compiling the test binaries ($JOBS at a time; cached for later runs) ..."
{ cut -d ' ' -f 2- "$TMP/cases" | sed 's/^ *//; s/ *-t [0-9]*//'
  echo "6,6,3 -f C4 -6"; echo "6,6,3 -f C4 -6 -p 0/5"; echo "63 -D 2 -6"; } | sort -u > "$TMP/sets"
# shellcheck disable=SC2016
if ! xargs -P "$JOBS" -I '{}' sh -c \
        '"$0" {} -B > /dev/null 2> "$1/build.$$.log" || { cat "$1/build.$$.log"; echo "check: cannot compile: volta {}"; exit 1; }' \
        "$VOLTA" "$TMP" < "$TMP/sets" >&2; then
    exit 1
fi

# 2. Run them.
pass=0; fail=0
while read -r expected args; do
    # shellcheck disable=SC2086  # args is a word list
    got=$("$VOLTA" $args 2> "$TMP/err" | sed -n 's/^graphs (N-classes) = //p')
    if [ "$got" = "$expected" ]; then
        pass=$((pass + 1))
        printf 'ok    %-36s %s\n' "$args" "$got"
    else
        fail=$((fail + 1))
        printf 'FAIL  %-36s expected %s, got %s\n' "$args" "$expected" "${got:-nothing}"
        sed 's/^/      /' "$TMP/err"
    fi
done < "$TMP/cases"

# 3. The output, not only the counts.
check() {   # check DESCRIPTION CONDITION...
    _d=$1; shift
    if "$@"; then pass=$((pass + 1)); printf 'ok    %s\n' "$_d"
    else fail=$((fail + 1)); printf 'FAIL  %s\n' "$_d"; fi
    total=$((total + 1))
}
"$VOLTA" 6,6,3 -f C4 -6 2> /dev/null > "$TMP/a.g6"
check "graph6 on stdout: 139 lines for 6,6,3 -f C4" \
    test "$(wc -l < "$TMP/a.g6" | tr -d ' ')" = 139
check "graph6 lines have the right header and length (n = 15)" \
    test "$(grep -c '^N.\{18\}$' "$TMP/a.g6")" = 139
"$VOLTA" 6,6,3 -f C4 -o "$TMP/b.g6" > /dev/null 2>&1
check "--output writes the same graphs" cmp -s "$TMP/a.g6" "$TMP/b.g6"
"$VOLTA" 6,6,3 -f C4 -6 -t 4 2> /dev/null | sort > "$TMP/c.g6"
sort "$TMP/a.g6" > "$TMP/a.sorted"
check "4 threads give the same graphs" cmp -s "$TMP/a.sorted" "$TMP/c.g6"
: > "$TMP/parts.g6"
for r in 0 1 2 3 4; do
    "$VOLTA" 6,6,3 -f C4 -6 -p "$r/5" 2> /dev/null >> "$TMP/parts.g6"
done
sort "$TMP/parts.g6" > "$TMP/parts.sorted"
check "--part 0/5 ... 4/5 together give every graph exactly once" cmp -s "$TMP/a.sorted" "$TMP/parts.sorted"
"$VOLTA" 63 -D 2 -6 2> /dev/null > "$TMP/d.g6"
check "graph6 with the 4-byte header: 6 lines for 63 -D 2 (n = 63)" \
    test "$(grep -c '^~??~.\{326\}$' "$TMP/d.g6")" = 6 -a "$(wc -l < "$TMP/d.g6" | tr -d ' ')" = 6
"$VOLTA" 6,6,3 -f C4 > "$TMP/short.txt" 2>&1
"$VOLTA" 6,6,3 -f C4 -v > "$TMP/long.txt" 2>&1
check "the default summary is short; -v adds the statistics and the counts by edges" \
    test "$(grep -c 'nodes visited\|by number of edges' "$TMP/short.txt")" = 0 \
      -a "$(grep -c 'nodes visited\|by number of edges' "$TMP/long.txt")" = 2

echo
if [ $fail = 0 ]; then
    echo "all $total tests passed"
    exit 0
fi
echo "$fail of $total tests FAILED"
exit 1
