/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file main.c
 * @brief Entry point: reads the instance, runs the chosen method, prints results.
 *
 * Program flow:
 *  1. parse the command line (options.c);
 *  2. load and preprocess the network (graph_load.c);
 *  3. start the MIP solver once, if the method needs it (exact_gurobi.c);
 *  4. for each independent run: solve (algorithms.c) and report (report.c).
 */
#include "algorithms.h"
#include "exact.h"
#include "graph.h"
#include "options.h"
#include "planar/entry.h"
#include "report.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Tells whether the configured method calls the MIP solver.
 *
 * @param params Method parameters.
 * @return false only for the greedy framework without exact cleaning.
 */
static bool needs_mip_solver(const SolverParams *params)
{
    return params->algorithm != ALGORITHM_GREEDY || params->cleaning == CLEANING_EXACT;
}

/**
 * @brief Program entry point.
 *
 * @param argc Argument count.
 * @param argv Arguments (see --help).
 * @return EXIT_SUCCESS, or EXIT_FAILURE if any run ended without a feasible cover.
 */
int main(int argc, char **argv)
{
    Options options;
    options_parse(argc, argv, &options);
    log_set_level(options.verbosity);

    if (options.params.algorithm == ALGORITHM_PLANAR) {
        /* The planar algorithm works on the untouched input and never calls
           the MIP solver; see planar/entry.c. */
        int status = planar_run_all(&options);
        options_free(&options);
        return status;
    }

    Graph *graph = graph_load(options.edges_path, options.od_path, options.contract_degree_two);
    if (options.info_only) {
        report_instance(stdout, graph, options.instance_name);
        graph_free(graph);
        options_free(&options);
        return EXIT_SUCCESS;
    }
    if (log_get_level() >= LOG_INFO) {
        report_instance(stderr, graph, options.instance_name);
    }

    ExactSolver *solver = nullptr;
    if (needs_mip_solver(&options.params)) {
        solver = exact_solver_create(options.threads, log_get_level() >= LOG_TRACE);
    }
    report_configuration(&options, solver != nullptr ? exact_solver_description() : nullptr);

    if (options.print_header) {
        report_table_header(stdout);
    }
    int status = EXIT_SUCCESS;
    for (int32_t run = 0; run < options.runs; run++) {
        uint64_t seed = options.seed + (uint64_t)run;
        RunResult result;
        run_algorithm(graph, &options.params, solver, seed, &result);

        report_table_row(stdout, &options, run, seed, &result);
        log_message(LOG_INFO, "run %d (seed %llu): %d counting stations%s, %.2fs", run,
                    (unsigned long long)seed, result.cost,
                    result.connected_pairs == 0 ? "" : " (INFEASIBLE)", result.time_total);
        if (result.connected_pairs != 0) {
            status = EXIT_FAILURE; /* should never happen; reported loudly if it does */
        }
        if (options.solution_dir != nullptr) {
            report_write_solution(&options, graph, seed, &result);
        }
        run_result_free(&result);
    }

    exact_solver_free(solver);
    graph_free(graph);
    options_free(&options);
    return status;
}
