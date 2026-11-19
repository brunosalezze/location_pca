/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file report.h
 * @brief Everything the program prints or writes: instance summary, result
 *        table (TSV) and solution files.
 */
#ifndef LOCATION_PCA_REPORT_H
#define LOCATION_PCA_REPORT_H

#include "algorithms.h"
#include "graph.h"
#include "options.h"

#include <stdint.h>
#include <stdio.h>

/**
 * @brief Prints the size of the instance and what the preprocessing did.
 *
 * @param stream  Destination (stdout for --info, stderr for logging).
 * @param graph   Preprocessed network.
 * @param name    Instance name.
 */
void report_instance(FILE *stream, const Graph *graph, const char *name);

/**
 * @brief Prints the method and its parameters (to the log, at LOG_INFO).
 *
 * @param options Parsed options.
 * @param solver  Solver description, or nullptr if no MIP solver is used.
 */
void report_configuration(const Options *options, const char *solver);

/**
 * @brief Prints the header line of the result table.
 *
 * Columns: instance, run, seed, algorithm, cost, lower_bound, gap_pct,
 * feasible, paths, iterations, mip_solves, time_total_s, time_mip_s, time_cut_s.
 *
 * @param stream Destination.
 */
void report_table_header(FILE *stream);

/**
 * @brief Prints one row of the result table.
 *
 * @param stream   Destination.
 * @param options  Parsed options (instance name, algorithm).
 * @param run      Run index (0-based).
 * @param seed     Seed of the run.
 * @param result   Outcome of the run.
 */
void report_table_row(FILE *stream, const Options *options, int32_t run, uint64_t seed,
                      const RunResult *result);

/**
 * @brief Writes the selected edges of a run to "<dir>/<instance>_<algorithm>_seed<seed>.txt".
 *
 * The file starts with '#' comment lines describing the run, followed by one
 * edge identifier (from the input edge file) per line. The directory is
 * created if needed.
 *
 * @param options Parsed options (solution_dir must be set).
 * @param graph   Network (for edge labels).
 * @param seed    Seed of the run.
 * @param result  Outcome of the run.
 */
void report_write_solution(const Options *options, const Graph *graph, uint64_t seed,
                           const RunResult *result);

#endif /* LOCATION_PCA_REPORT_H */
