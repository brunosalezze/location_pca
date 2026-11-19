/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file greedy.c
 * @brief Implementation of the greedy covering heuristic.
 */
#include "greedy.h"

#include "util.h"

#include <stdlib.h>

void greedy_cover(Cover *cover, int32_t max_uncovered, Rng *rng)
{
    if (cover->num_uncovered <= max_uncovered) {
        return;
    }
    int32_t num_edges = cover->graph->num_edges;
    int32_t *best_edges = xmalloc_array((size_t)num_edges, sizeof *best_edges);

    while (cover->num_uncovered > max_uncovered) {
        int32_t best_value = 0;
        uint32_t ties = 0;
        for (int32_t e = 0; e < num_edges; e++) {
            int32_t value = cover->uncovered_paths_through[e];
            if (value > best_value) {
                best_value = value;
                ties = 0;
                best_edges[ties++] = e;
            } else if (value == best_value && value > 0) {
                best_edges[ties++] = e;
            }
        }
        /* An uncovered path always has an unselected edge with value >= 1. */
        if (ties == 0) {
            fatal("internal error: greedy_cover found no useful edge");
        }
        cover_add_edge(cover, best_edges[rng_below(rng, ties)]);
    }
    free(best_edges);
}
