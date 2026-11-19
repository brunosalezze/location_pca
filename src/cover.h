/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file cover.h
 * @brief A cover C ⊆ A (set of edges with counting stations) and its coverage counters.
 *
 * A path is *covered* when at least one of its edges is selected. To make the
 * covering heuristics fast, the structure keeps, for the current path pool P1:
 *  - path_hits[p]: how many selected edges path p contains (0 = uncovered);
 *  - uncovered_paths_through[e]: how many uncovered paths contain edge e
 *    (|P*_e| in Algorithms 2 and 3 of the paper);
 *  - num_uncovered: how many paths of P1 are uncovered (|P*|).
 * Adding or removing an edge updates the counters incrementally.
 *
 * The counters refer to the paths that existed at the last cover_reset() or
 * cover_sync(). Call cover_sync() after adding paths to the pool before
 * modifying the cover again.
 *
 * Locked edges (graph.h) are part of every cover: cover_reset() adds them and
 * the cleaning procedures never remove them.
 */
#ifndef LOCATION_PCA_COVER_H
#define LOCATION_PCA_COVER_H

#include "graph.h"
#include "path_pool.h"

#include <stddef.h>
#include <stdint.h>

/** Selected edges plus coverage counters relative to a path pool. */
typedef struct Cover {
    const Graph *graph;   /**< Network (not owned). */
    const PathPool *pool; /**< Path set the counters refer to (not owned). */

    uint8_t *selected;    /**< [num_edges] 1 if the edge has a counting station. */
    int32_t num_selected; /**< |C|, the objective value. */

    int32_t *uncovered_paths_through; /**< [num_edges] Uncovered paths containing each edge. */
    int32_t *path_hits;               /**< [paths] Selected edges on each path. */
    size_t path_hits_capacity;        /**< Allocated length of path_hits. */
    int32_t num_paths;                /**< Number of pool paths the counters account for. */
    int32_t num_uncovered;            /**< Paths with path_hits == 0. */
} Cover;

/**
 * @brief Creates an empty cover (no edge selected, no path known).
 *
 * @param graph Network; must outlive the cover.
 * @param pool  Path pool; must outlive the cover.
 * @return Newly allocated cover; call cover_reset() before using it.
 */
[[nodiscard]] Cover *cover_create(const Graph *graph, const PathPool *pool);

/**
 * @brief Releases a cover.
 *
 * @param cover Cover to release, or nullptr.
 */
void cover_free(Cover *cover);

/**
 * @brief Empties the cover, re-synchronizes it with the whole pool, then adds the locked edges.
 *
 * Equivalent to clearSolution() in the original code. Cost: O(|A| + |P1|)
 * plus the paths through locked edges.
 *
 * @param cover Cover to reset.
 */
void cover_reset(Cover *cover);

/**
 * @brief Accounts for the paths appended to the pool since the last reset/sync.
 *
 * Inner workings: for each new path, counts its selected edges. If there
 * are none, the path is uncovered and the counters of its edges increase.
 *
 * @param cover Cover to update.
 */
void cover_sync(Cover *cover);

/**
 * @brief Selects @p edge (adds a counting station).
 *
 * Inner workings: every path through the edge gains one hit. A path that was
 * uncovered becomes covered, so each of its edges loses one uncovered path
 * and num_uncovered decreases.
 *
 * @param cover The cover (must be synchronized with the pool).
 * @param edge  An edge that is not selected yet.
 */
void cover_add_edge(Cover *cover, int32_t edge);

/**
 * @brief Unselects @p edge (exact inverse of cover_add_edge()).
 *
 * @param cover The cover (must be synchronized with the pool).
 * @param edge  A selected edge.
 */
void cover_remove_edge(Cover *cover, int32_t edge);

/**
 * @brief Lists the redundant edges of the cover (the set R_C).
 *
 * R_C is used by the Cleaning procedures: a selected, non-locked edge is
 * redundant when every path of P1 through it contains at least one other
 * selected edge, so removing that edge alone keeps every path covered.
 *
 * @param cover The cover.
 * @param out   Output array with room for num_edges entries.
 * @return Number of edges written to @p out, in increasing index order.
 */
int32_t cover_collect_redundant(const Cover *cover, int32_t *out);

/**
 * @brief Counts the O-D pairs still connected when the selected edges are removed.
 *
 * This is an independent feasibility check that does not rely on the path pool:
 * a cover is feasible for the TCLP if and only if no two O-D nodes lie in the
 * same connected component of G \ C. Uses union-find, O(|A| α(|N|)).
 *
 * @param graph    The network.
 * @param selected [num_edges] Selection flags.
 * @return Number of O-D pairs joined by some path avoiding the selected edges
 *         (0 means the cover is feasible).
 */
int64_t cover_count_connected_od_pairs(const Graph *graph, const uint8_t *selected);

#endif /* LOCATION_PCA_COVER_H */
