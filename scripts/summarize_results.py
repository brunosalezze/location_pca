#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
# SPDX-License-Identifier: GPL-3.0-or-later
"""Summarize location_pca result files per instance (standard library only).

For every instance prints: number of runs, best / average / worst cost, the
relative standard deviation of the cost, the best lower bound, the gap between
the best cost and the best lower bound, the average running time and its
relative standard deviation, and how many runs were infeasible (should be 0).

Usage:
    scripts/summarize_results.py results/pca.tsv [more.tsv ...] [--markdown]
"""

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict


def read_rows(paths):
    """Yield the result rows (as dicts) of all given TSV files."""
    for path in paths:
        with open(path, newline="") as handle:
            yield from csv.DictReader(handle, delimiter="\t")


def relative_std(values):
    """Population standard deviation as a percentage of the mean (0 if undefined)."""
    mean = statistics.fmean(values)
    return 100.0 * statistics.pstdev(values) / mean if mean else 0.0


def summarize(rows):
    """Group rows by instance and compute the summary statistics."""
    groups = defaultdict(list)
    for row in rows:
        groups[row["instance"]].append(row)

    summary = []
    for name in sorted(groups):
        runs = groups[name]
        costs = [int(r["cost"]) for r in runs]
        times = [float(r["time_total_s"]) for r in runs]
        bounds = [float(r["lower_bound"]) for r in runs if r["lower_bound"] != "NA"]
        best_bound = max(bounds) if bounds else math.nan
        best = min(costs)
        gap = 100.0 * (best - best_bound) / best_bound if bounds and best_bound > 0 else math.nan
        summary.append({
            "instance": name,
            "runs": len(runs),
            "best": best,
            "average": statistics.fmean(costs),
            "worst": max(costs),
            "cost_rsd_pct": relative_std(costs),
            "best_lb": best_bound,
            "gap_pct": gap,
            "time_avg_s": statistics.fmean(times),
            "time_rsd_pct": relative_std(times),
            "infeasible": sum(1 for r in runs if r["feasible"] != "1"),
        })
    return summary


def format_value(value):
    if isinstance(value, float):
        return "NA" if math.isnan(value) else f"{value:.2f}"
    return str(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+", help="TSV files written by location_pca")
    parser.add_argument("--markdown", action="store_true", help="print a Markdown table")
    args = parser.parse_args()

    summary = summarize(read_rows(args.files))
    if not summary:
        sys.exit("no result rows found")
    columns = list(summary[0].keys())
    if args.markdown:
        print("| " + " | ".join(columns) + " |")
        print("|" + "---|" * len(columns))
        for entry in summary:
            print("| " + " | ".join(format_value(entry[c]) for c in columns) + " |")
    else:
        print("\t".join(columns))
        for entry in summary:
            print("\t".join(format_value(entry[c]) for c in columns))


if __name__ == "__main__":
    main()
