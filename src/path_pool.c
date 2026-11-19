/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file path_pool.c
 * @brief Implementation of the O-D path set.
 */
#include "path_pool.h"

#include <stdlib.h>
#include <string.h>

PathPool *path_pool_create(int32_t num_graph_edges)
{
    PathPool *pool = xcalloc(1, sizeof *pool);
    pool->num_graph_edges = num_graph_edges;
    pool->path_capacity = 1024;
    pool->path_start = xmalloc_array(pool->path_capacity + 1, sizeof *pool->path_start);
    pool->path_start[0] = 0;
    pool->path_edges_capacity = 16384;
    pool->path_edges = xmalloc_array(pool->path_edges_capacity, sizeof *pool->path_edges);
    pool->paths_through_edge = xcalloc((size_t)num_graph_edges, sizeof *pool->paths_through_edge);
    return pool;
}

void path_pool_free(PathPool *pool)
{
    if (pool == nullptr) {
        return;
    }
    for (int32_t e = 0; e < pool->num_graph_edges; e++) {
        intvec_free(&pool->paths_through_edge[e]);
    }
    free(pool->paths_through_edge);
    free(pool->path_start);
    free(pool->path_edges);
    free(pool);
}

int32_t path_pool_add(PathPool *pool, const int32_t *edges, int32_t length)
{
    if ((size_t)pool->num_paths == pool->path_capacity) {
        pool->path_capacity *= 2;
        pool->path_start =
            xrealloc_array(pool->path_start, pool->path_capacity + 1, sizeof *pool->path_start);
    }
    size_t needed = pool->path_edges_count + (size_t)length;
    if (needed > pool->path_edges_capacity) {
        while (pool->path_edges_capacity < needed) {
            pool->path_edges_capacity *= 2;
        }
        pool->path_edges =
            xrealloc_array(pool->path_edges, pool->path_edges_capacity, sizeof *pool->path_edges);
    }

    int32_t path = pool->num_paths++;
    memcpy(pool->path_edges + pool->path_edges_count, edges, (size_t)length * sizeof *edges);
    pool->path_edges_count = needed;
    pool->path_start[path + 1] = needed;
    for (int32_t i = 0; i < length; i++) {
        intvec_push(&pool->paths_through_edge[edges[i]], path);
    }
    return path;
}
