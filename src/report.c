/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file report.c
 * @brief Output formatting: instance summary, TSV result table, solution files.
 */
#define _POSIX_C_SOURCE 200809L

#include "report.h"

#include "util.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/**
 * @brief Name of a spread criterion as used on the command line.
 *
 * @param spread Spread criterion.
 * @return Static string.
 */
static const char *spread_name(SpreadCriterion spread)
{
    return spread == SPREAD_MLC ? "mlc" : "random";
}

/**
 * @brief Name of a cleaning method as used on the command line.
 *
 * @param cleaning Cleaning method.
 * @return Static string.
 */
static const char *cleaning_name(CleaningMethod cleaning)
{
    switch (cleaning) {
    case CLEANING_NONE:
        return "none";
    case CLEANING_RANDOM:
        return "random";
    case CLEANING_GREEDY:
        return "greedy";
    case CLEANING_EXACT:
        return "exact";
    }
    return "?";
}

void report_instance(FILE *stream, const Graph *graph, const char *name)
{
    const GraphStats *stats = &graph->stats;
    fprintf(stream,
            "instance %s: %d nodes, %d edges, %d O-D nodes (%" PRId64 " pairs), %d locked edges\n",
            name, graph->num_nodes, graph->num_edges, graph->num_od_nodes,
            graph_num_od_pairs(graph), graph->num_locked_edges);
    fprintf(stream,
            "  preprocessing: %d edge lines; removed %d self-loops, %d parallel edges, "
            "%d ring edges; contracted %d degree-2 nodes\n",
            stats->raw_edges, stats->self_loops_removed, stats->parallel_edges_removed,
            stats->ring_edges_removed, stats->nodes_contracted);
    if (stats->duplicate_od_entries > 0 || stats->isolated_od_nodes > 0) {
        fprintf(stream, "  note: %d duplicate O-D entries ignored, %d O-D nodes without edges\n",
                stats->duplicate_od_entries, stats->isolated_od_nodes);
    }
}

void report_configuration(const Options *options, const char *solver)
{
    const SolverParams *params = &options->params;
    char alpha[48] = "";
    if (params->algorithm == ALGORITHM_PCA) {
        snprintf(alpha, sizeof alpha, ", alpha0 %g", params->alpha0);
    }
    log_message(
        LOG_INFO, "method: %s%s, spread %s, cleaning %s, paths per O-D node %d then %d (0 = all)",
        algorithm_name(params->algorithm), alpha, spread_name(params->spread),
        cleaning_name(params->cleaning), params->initial_paths_per_root, params->paths_per_root);
    char threads[32] = "";
    if (solver != nullptr) {
        snprintf(threads, sizeof threads, " (%d thread%s)", options->threads,
                 options->threads == 1 ? "" : "s");
    }
    log_message(LOG_INFO, "limits: run %gs, exact covering %gs, exact cleaning %gs; solver: %s%s",
                params->time_limit_total, params->time_limit_iteration, params->time_limit_cleaning,
                solver != nullptr ? solver : "not used", threads);
}

void report_table_header(FILE *stream)
{
    fputs("instance\trun\tseed\talgorithm\tcost\tlower_bound\tgap_pct\tfeasible\tpaths\t"
          "iterations\tmip_solves\ttime_total_s\ttime_mip_s\ttime_cut_s\n",
          stream);
}

void report_table_row(FILE *stream, const Options *options, int32_t run, uint64_t seed,
                      const RunResult *result)
{
    char bound[32] = "NA";
    char gap[32] = "NA";
    if (!isnan(result->lower_bound)) {
        snprintf(bound, sizeof bound, "%.0f", result->lower_bound);
        if (result->lower_bound > 0) {
            snprintf(gap, sizeof gap, "%.2f",
                     100.0 * (result->cost - result->lower_bound) / result->lower_bound);
        }
    }
    fprintf(stream, "%s\t%d\t%" PRIu64 "\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n",
            options->instance_name, run, seed, algorithm_name(options->params.algorithm),
            result->cost, bound, gap, result->connected_pairs == 0, result->num_paths,
            result->iterations, result->exact_solves, result->time_total, result->time_solver,
            result->time_separation);
    fflush(stream);
}

void report_write_solution(const Options *options, const Graph *graph, uint64_t seed,
                           const RunResult *result)
{
    if (mkdir(options->solution_dir, 0777) != 0 && errno != EEXIST) {
        fatal("cannot create directory '%s': %s", options->solution_dir, strerror(errno));
    }
    const char *algorithm = algorithm_name(options->params.algorithm);
    size_t size =
        strlen(options->solution_dir) + strlen(options->instance_name) + strlen(algorithm) + 48;
    char *path = xmalloc(size);
    snprintf(path, size, "%s/%s_%s_seed%" PRIu64 ".txt", options->solution_dir,
             options->instance_name, algorithm, seed);

    FILE *file = fopen(path, "w");
    if (file == nullptr) {
        fatal("cannot write '%s': %s", path, strerror(errno));
    }
    fprintf(file, "# location_pca %s solution\n", LOCATION_PCA_VERSION);
    fprintf(file, "# instance: %s\n# algorithm: %s\n# seed: %" PRIu64 "\n", options->instance_name,
            algorithm, seed);
    fprintf(file, "# counting stations: %d\n# feasible: %s\n", result->cost,
            result->connected_pairs == 0 ? "yes" : "no");
    if (!isnan(result->lower_bound)) {
        fprintf(file, "# lower bound: %.0f\n", result->lower_bound);
    }
    fputs("# One edge identifier (first column of the edge file) per line.\n", file);
    for (int32_t e = 0; e < graph->num_edges; e++) {
        if (result->selected[e]) {
            fprintf(file, "%s\n", graph->edge_label[e]);
        }
    }
    if (fclose(file) != 0) {
        fatal("error while writing '%s'", path);
    }
    free(path);
}
