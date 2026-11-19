/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file separation.h
 * @brief Cut phase (Section 4.2 and Algorithm 5 of the paper): finding O-D paths
 *        that avoid every counting station of the current cover.
 *
 * For each O-D node i (the *root*), a search tree is grown in G' = (N, A \ C).
 * Starting from the root, one frontier edge (an edge from a visited node to an
 * unvisited node) is added to the tree at each step. Which one is decided by the
 * *spread criterion*. Whenever the tree reaches an O-D partner j of i (see
 * graph_od_partner_count()), the tree path i -> j is uncovered by C and is
 * added to the path pool as a new covering constraint. When no root reaches
 * any of its partners, C separates every O-D pair and is feasible.
 *
 * Only *direct* paths are generated: the tree does not grow beyond an O-D
 * node, so the intermediate nodes of every path are intersections. This is
 * enough to decide feasibility. A path from i to j through another O-D node k
 * starts with a direct path from i to the first O-D node it meets; if the long
 * path avoids the cover, so does that first part, which is itself an uncovered
 * O-D path. As a covering constraint the long path is implied by the short one
 * (its edge set is a superset), so P1 gets fewer and stronger rows, and the
 * covering model splits into independent blocks, one per region of
 * intersections enclosed by O-D nodes. The program used for the article let
 * the tree grow through O-D nodes; README.md compares the two.
 */
#ifndef LOCATION_PCA_SEPARATION_H
#define LOCATION_PCA_SEPARATION_H

#include "cover.h"
#include "graph.h"
#include "path_pool.h"
#include "rng.h"

#include <stdint.h>

/** How the next frontier edge is chosen while growing a search tree. */
typedef enum SpreadCriterion {
    /**
     * Minimum Link Count (MLC, the paper's default): the frontier edge with
     * the fewest links from the root is chosen, so the tree is a breadth-first
     * search tree and paths have as few links as possible. Ties are broken in
     * favour of the edge used by the fewest stored paths, then uniformly at random.
     */
    SPREAD_MLC,
    /**
     * Non-Discriminated (ND): depth is ignored. The frontier edge used by the
     * fewest stored paths is chosen, ties uniformly at random (this is how
     * the original implementation realized "select next node randomly").
     */
    SPREAD_RANDOM,
} SpreadCriterion;

/** Reusable scratch memory for the Cut phase (one per run). */
typedef struct Separator Separator;

/**
 * @brief Allocates the scratch buffers needed to search @p graph.
 *
 * @param graph Network; must outlive the separator.
 * @return Newly allocated separator; release it with separator_free().
 */
[[nodiscard]] Separator *separator_create(const Graph *graph);

/**
 * @brief Releases a separator.
 *
 * @param separator Separator to release, or nullptr.
 */
void separator_free(Separator *separator);

/**
 * @brief Runs the Cut phase once: adds to @p pool new O-D paths avoiding @p cover.
 *
 * Inner workings, for each O-D root in O-D list order:
 *  1. the root is visited and its non-selected incident edges become the frontier;
 *  2. repeatedly, an edge is chosen from the frontier according to @p criterion,
 *     its unvisited endpoint v joins the tree, the frontier edges that also lead
 *     to v are dropped and v's own non-selected edges towards unvisited nodes
 *     are added, unless v is an O-D node (direct paths, see above);
 *  3. if v is an O-D partner of the root, the tree path root -> v is stored;
 *     the search of this root stops after @p max_paths_per_root paths or when
 *     the frontier is empty.
 * Because the edge-usage counts change as paths are added, later trees tend to
 * use different edges ("spreading" the paths over the network).
 *
 * Complexity: O(|T| · |N| · F), where F is the frontier size, usually small.
 *
 * @param separator          Scratch memory.
 * @param cover              Current cover C; only its selection flags are read.
 * @param pool               Path pool that receives the new paths.
 * @param criterion          Spread criterion.
 * @param max_paths_per_root Maximum number of new paths per root (0 = no limit).
 * @param rng                Random generator used to break ties.
 * @return Number of paths added. 0 means the cover separates every O-D pair.
 */
int32_t separator_run(Separator *separator, const Cover *cover, PathPool *pool,
                      SpreadCriterion criterion, int32_t max_paths_per_root, Rng *rng);

#endif /* LOCATION_PCA_SEPARATION_H */
