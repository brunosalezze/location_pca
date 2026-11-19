/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file cleaning.c
 * @brief Implementation of the RC, GC and EC cleaning procedures.
 */
#include "cleaning.h"

#include "greedy.h"
#include "util.h"

#include <stdlib.h>

/**
 * @brief RC: removes random redundant edges one by one until none is left.
 *
 * @param cover     Cover to clean.
 * @param redundant Scratch array with room for num_edges entries.
 * @param rng       Random generator.
 */
static void clean_random(Cover *cover, int32_t *redundant, Rng *rng)
{
    for (;;) {
        int32_t count = cover_collect_redundant(cover, redundant);
        if (count == 0) {
            return;
        }
        cover_remove_edge(cover, redundant[rng_below(rng, (uint32_t)count)]);
    }
}

/**
 * @brief EC: removes all of R_C, then completes the cover with the MIP solver.
 *
 * @param cover      Cover to clean.
 * @param redundant  Edges of R_C.
 * @param count      Size of R_C (> 0).
 * @param solver     MIP solver.
 * @param time_limit Solver time limit in seconds.
 * @return Outcome of the solve (runtime and whether the solver was called).
 */
static ExactResult clean_exact(Cover *cover, const int32_t *redundant, int32_t count,
                               ExactSolver *solver, double time_limit)
{
    uint8_t *warm_start = xcalloc((size_t)cover->graph->num_edges, sizeof *warm_start);
    for (int32_t i = 0; i < count; i++) {
        warm_start[redundant[i]] = 1;
        cover_remove_edge(cover, redundant[i]);
    }
    ExactResult result = exact_cover(solver, cover, time_limit, warm_start, 0.0);
    if (!result.has_solution) {
        for (int32_t i = 0; i < count; i++) {
            cover_add_edge(cover, redundant[i]);
        }
    }
    free(warm_start);
    return result;
}

CleaningStats clean_cover(Cover *cover, CleaningMethod method, ExactSolver *solver,
                          double time_limit, Rng *rng)
{
    CleaningStats stats = {0};
    if (method == CLEANING_NONE) {
        return stats;
    }
    int32_t *redundant = xmalloc_array((size_t)cover->graph->num_edges, sizeof *redundant);

    if (method == CLEANING_RANDOM) {
        clean_random(cover, redundant, rng);
    } else {
        int32_t count = cover_collect_redundant(cover, redundant);
        if (count > 0) {
            if (method == CLEANING_GREEDY) {
                for (int32_t i = 0; i < count; i++) {
                    cover_remove_edge(cover, redundant[i]);
                }
                greedy_cover(cover, 0, rng);
            } else {
                ExactResult solve = clean_exact(cover, redundant, count, solver, time_limit);
                stats.solver_time = solve.runtime;
                stats.solver_calls = solve.solver_calls;
            }
        }
    }
    free(redundant);
    return stats;
}
