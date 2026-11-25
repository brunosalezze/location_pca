/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file options.h
 * @brief Command-line options.
 */
#ifndef LOCATION_PCA_OPTIONS_H
#define LOCATION_PCA_OPTIONS_H

#include "algorithms.h"
#include "util.h"

#include <stdint.h>

/** Program version, reported by --version. */
#define LOCATION_PCA_VERSION "1.0.0"

/** Everything the user can configure from the command line. */
typedef struct Options {
    char *edges_path;         /**< Edge file. */
    char *od_path;            /**< O-D node file. */
    char *instance_name;      /**< Name used in the output. */
    SolverParams params;      /**< Algorithm and its parameters. */
    int threads;              /**< MIP solver threads. */
    int32_t runs;             /**< Number of independent runs. */
    uint64_t seed;            /**< Seed of the first run (run k uses seed + k). */
    char *solution_dir;       /**< Directory for solution files, or nullptr. */
    bool print_header;        /**< Print the TSV header line. */
    bool info_only;           /**< Only print instance statistics. */
    bool contract_degree_two; /**< Preprocessing step 3 (see graph_load()). */
    bool planarity_only;      /**< With -a planar: only test planarity. */
    bool planar_paper_only;   /**< With -a planar: only the Klein-Marx method. */
    LogLevel verbosity;       /**< Logging level. */
} Options;

/**
 * @brief Parses the command line into @p options.
 *
 * Inner workings: options are read with getopt_long(). Values given
 * explicitly override the defaults of the selected algorithm, whatever their
 * position relative to --algorithm. The positional arguments are either an
 * instance directory (containing edges.txt and od_nodes.txt) or an edge file
 * followed by an O-D file.
 *
 * Prints the help or version text and exits for --help / --version. Prints
 * an error and exits with status 2 for invalid arguments.
 *
 * @param argc    Argument count from main().
 * @param argv    Argument vector from main().
 * @param options Output; release with options_free().
 */
void options_parse(int argc, char **argv, Options *options);

/**
 * @brief Releases the strings owned by @p options.
 *
 * @param options Options filled by options_parse().
 */
void options_free(Options *options);

#endif /* LOCATION_PCA_OPTIONS_H */
