#!/bin/sh
# SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
# SPDX-License-Identifier: GPL-3.0-or-later
# Runs location_pca on a list of instances and collects every result line in
# one tab-separated file. Instances are independent processes, so several
# can run at the same time (-j). Keep -j 1 when measuring running times.
#
# Usage:
#   scripts/run_experiments.sh [-j JOBS] [-o OUTPUT] [-i "AC AL ..."] [-- SOLVER OPTIONS]
#
# Examples:
#   # Table 3, heuristic setting: PCA, alpha0 = 1/16, TL1 = 500 s, 15 runs per instance
#   scripts/run_experiments.sh -o results/pca.tsv -- --runs 15 --seed 1
#   # Table 3, exact setting: PCA, alpha0 = 1, TL1 = 1200 s
#   scripts/run_experiments.sh -o results/pca_exact.tsv -- --alpha0 1 --iteration-time-limit 1200
#   # Table 2: greedy framework with each cleaning procedure, 30 runs
#   for c in none random greedy exact; do
#     scripts/run_experiments.sh -j 4 -o results/greedy_$c.tsv -- -a greedy --cleaning $c --runs 30
#   done
set -eu

root=$(cd "$(dirname "$0")/.." && pwd)
program="$root/bin/location_pca"
jobs=1
output="results/results.tsv"
instances="AC AL AM AP BA CE ES GOeDF MA MG MS MT PA PB PE PI PR RJ RN RO RR RS SC SE SP TO"

usage() {
    sed -n '2,/^[^#]/{/^# SPDX-/d;/^#/p;}' "$0" | sed 's/^# \{0,1\}//'
    exit "$1"
}

while getopts "j:o:i:h" flag; do
    case $flag in
    j) jobs=$OPTARG ;;
    o) output=$OPTARG ;;
    i) instances=$OPTARG ;;
    h) usage 0 ;;
    *) usage 2 ;;
    esac
done
shift $((OPTIND - 1))
[ "${1:-}" = "--" ] && shift

if [ ! -x "$program" ]; then
    echo "run_experiments: $program not found; run 'make' first" >&2
    exit 1
fi

mkdir -p "$(dirname "$output")"
partial=$(mktemp -d)
trap 'rm -rf "$partial"' EXIT

# One process per instance; each writes its own file, merged at the end.
for name in $instances; do
    echo "$name"
done | xargs -P "$jobs" -I{} sh -c '
    echo "running {}" >&2
    "$0" --quiet --no-header "$@" "'"$root"'/instances/{}" > "'"$partial"'/{}.tsv"
' "$program" "$@"

{
    printf 'instance\trun\tseed\talgorithm\tcost\tlower_bound\tgap_pct\tfeasible\tpaths\t'
    printf 'iterations\tmip_solves\ttime_total_s\ttime_mip_s\ttime_cut_s\n'
    for name in $instances; do
        cat "$partial/$name.tsv"
    done
} > "$output"
echo "results written to $output" >&2
