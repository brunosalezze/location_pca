/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file path_pool.h
 * @brief The growing set P1 of O-D paths, i.e. the rows of the set covering model.
 *
 * Every path found by the Cut phase (separation.h) is appended here and never
 * removed during a run. A path is stored as the list of its edge indices, and
 * for every edge the pool also keeps the list of paths that contain it. Both
 * directions are needed:
 *  - path -> edges: to build the covering constraints and to update coverage counters;
 *  - edge -> paths: to know which rows a column (edge) covers.
 *
 * Storage is contiguous and grows on demand. The original code instead used
 * pre-allocated pools of 10 million path-tree nodes and 50 million list nodes
 * (about 3 GB reserved at start-up, whatever the instance size).
 */
#ifndef LOCATION_PCA_PATH_POOL_H
#define LOCATION_PCA_PATH_POOL_H

#include "util.h"

#include <stddef.h>
#include <stdint.h>

/** Set of O-D paths with an edge -> paths index. */
typedef struct PathPool {
    int32_t num_paths;       /**< Number of stored paths. */
    int32_t num_graph_edges; /**< Number of edges of the underlying graph. */

    size_t *path_start;   /**< [num_paths + 1] Offsets of each path in path_edges. */
    size_t path_capacity; /**< Allocated length of path_start, minus one. */

    int32_t *path_edges;        /**< Edge indices of all paths, concatenated. */
    size_t path_edges_count;    /**< Used length of path_edges. */
    size_t path_edges_capacity; /**< Allocated length of path_edges. */

    IntVec *paths_through_edge; /**< [num_graph_edges] Indices of the paths containing each edge. */
} PathPool;

/**
 * @brief Creates an empty pool for a graph with @p num_graph_edges edges.
 *
 * @param num_graph_edges Number of edges of the graph.
 * @return Newly allocated pool; release it with path_pool_free().
 */
[[nodiscard]] PathPool *path_pool_create(int32_t num_graph_edges);

/**
 * @brief Releases a pool.
 *
 * @param pool Pool to release, or nullptr.
 */
void path_pool_free(PathPool *pool);

/**
 * @brief Appends a path.
 *
 * @param pool   The pool.
 * @param edges  Edge indices of the path (any order; no repeated edge).
 * @param length Number of edges, at least 1.
 * @return Index of the new path.
 */
int32_t path_pool_add(PathPool *pool, const int32_t *edges, int32_t length);

/**
 * @brief Number of edges of path @p path.
 *
 * @param pool The pool.
 * @param path Path index.
 * @return Path length in edges.
 */
static inline int32_t path_pool_length(const PathPool *pool, int32_t path)
{
    return (int32_t)(pool->path_start[path + 1] - pool->path_start[path]);
}

/**
 * @brief Edges of path @p path.
 *
 * @param pool The pool.
 * @param path Path index.
 * @return Pointer to path_pool_length() edge indices, valid until the next path_pool_add().
 */
static inline const int32_t *path_pool_edges(const PathPool *pool, int32_t path)
{
    return pool->path_edges + pool->path_start[path];
}

/**
 * @brief Number of stored paths that contain @p edge.
 *
 * This is |P1_a| in the paper's notation and is used by the spread criteria
 * to prefer edges that are still rarely used.
 *
 * @param pool The pool.
 * @param edge Edge index.
 * @return Number of paths through the edge.
 */
static inline int32_t path_pool_edge_usage(const PathPool *pool, int32_t edge)
{
    return (int32_t)pool->paths_through_edge[edge].count;
}

#endif /* LOCATION_PCA_PATH_POOL_H */
