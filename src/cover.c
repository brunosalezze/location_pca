/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file cover.c
 * @brief Implementation of the cover and its incremental coverage counters.
 */
#include "cover.h"

#include "util.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

Cover *cover_create(const Graph *graph, const PathPool *pool)
{
    Cover *cover = xcalloc(1, sizeof *cover);
    cover->graph = graph;
    cover->pool = pool;
    cover->selected = xcalloc((size_t)graph->num_edges, sizeof *cover->selected);
    cover->uncovered_paths_through =
        xcalloc((size_t)graph->num_edges, sizeof *cover->uncovered_paths_through);
    cover->path_hits_capacity = 1024;
    cover->path_hits = xcalloc(cover->path_hits_capacity, sizeof *cover->path_hits);
    return cover;
}

void cover_free(Cover *cover)
{
    if (cover == nullptr) {
        return;
    }
    free(cover->selected);
    free(cover->uncovered_paths_through);
    free(cover->path_hits);
    free(cover);
}

/**
 * @brief Grows path_hits so that it can hold a counter for every pool path.
 *
 * @param cover The cover.
 */
static void reserve_path_counters(Cover *cover)
{
    size_t needed = (size_t)cover->pool->num_paths;
    if (needed <= cover->path_hits_capacity) {
        return;
    }
    size_t capacity = cover->path_hits_capacity;
    while (capacity < needed) {
        capacity *= 2;
    }
    cover->path_hits = xrealloc_array(cover->path_hits, capacity, sizeof *cover->path_hits);
    cover->path_hits_capacity = capacity;
}

void cover_reset(Cover *cover)
{
    const Graph *graph = cover->graph;
    const PathPool *pool = cover->pool;
    reserve_path_counters(cover);

    memset(cover->selected, 0, (size_t)graph->num_edges * sizeof *cover->selected);
    memset(cover->path_hits, 0, (size_t)pool->num_paths * sizeof *cover->path_hits);
    for (int32_t e = 0; e < graph->num_edges; e++) {
        cover->uncovered_paths_through[e] = path_pool_edge_usage(pool, e);
    }
    cover->num_selected = 0;
    cover->num_paths = pool->num_paths;
    cover->num_uncovered = pool->num_paths;

    for (int32_t i = 0; i < graph->num_locked_edges; i++) {
        cover_add_edge(cover, graph->locked_edges[i]);
    }
}

void cover_sync(Cover *cover)
{
    const PathPool *pool = cover->pool;
    reserve_path_counters(cover);
    for (int32_t p = cover->num_paths; p < pool->num_paths; p++) {
        const int32_t *edges = path_pool_edges(pool, p);
        int32_t length = path_pool_length(pool, p);
        int32_t hits = 0;
        for (int32_t i = 0; i < length; i++) {
            hits += cover->selected[edges[i]];
        }
        cover->path_hits[p] = hits;
        if (hits == 0) {
            cover->num_uncovered++;
            for (int32_t i = 0; i < length; i++) {
                cover->uncovered_paths_through[edges[i]]++;
            }
        }
    }
    cover->num_paths = pool->num_paths;
}

void cover_add_edge(Cover *cover, int32_t edge)
{
    assert(!cover->selected[edge]);
    assert(cover->num_paths == cover->pool->num_paths);
    const PathPool *pool = cover->pool;
    cover->selected[edge] = 1;
    cover->num_selected++;

    const IntVec *paths = &pool->paths_through_edge[edge];
    for (size_t i = 0; i < paths->count; i++) {
        int32_t path = paths->items[i];
        if (cover->path_hits[path]++ == 0) {
            const int32_t *edges = path_pool_edges(pool, path);
            int32_t length = path_pool_length(pool, path);
            for (int32_t k = 0; k < length; k++) {
                cover->uncovered_paths_through[edges[k]]--;
            }
            cover->num_uncovered--;
        }
    }
}

void cover_remove_edge(Cover *cover, int32_t edge)
{
    assert(cover->selected[edge]);
    assert(cover->num_paths == cover->pool->num_paths);
    const PathPool *pool = cover->pool;
    cover->selected[edge] = 0;
    cover->num_selected--;

    const IntVec *paths = &pool->paths_through_edge[edge];
    for (size_t i = 0; i < paths->count; i++) {
        int32_t path = paths->items[i];
        if (--cover->path_hits[path] == 0) {
            const int32_t *edges = path_pool_edges(pool, path);
            int32_t length = path_pool_length(pool, path);
            for (int32_t k = 0; k < length; k++) {
                cover->uncovered_paths_through[edges[k]]++;
            }
            cover->num_uncovered++;
        }
    }
}

/**
 * @brief Tells whether removing @p edge keeps every path of P1 covered.
 *
 * Inner workings: the edge is redundant when no path through it has a hit
 * count of 1, i.e. no path depends on this edge alone.
 *
 * @param cover The cover.
 * @param edge  Any edge index.
 * @return true if @p edge is selected, not locked and redundant.
 */
static bool edge_is_redundant(const Cover *cover, int32_t edge)
{
    if (!cover->selected[edge] || cover->graph->edge_is_locked[edge]) {
        return false;
    }
    const IntVec *paths = &cover->pool->paths_through_edge[edge];
    for (size_t i = 0; i < paths->count; i++) {
        if (cover->path_hits[paths->items[i]] == 1) {
            return false; /* this edge is the only one covering that path */
        }
    }
    return true;
}

int32_t cover_collect_redundant(const Cover *cover, int32_t *out)
{
    int32_t count = 0;
    for (int32_t e = 0; e < cover->graph->num_edges; e++) {
        if (edge_is_redundant(cover, e)) {
            out[count++] = e;
        }
    }
    return count;
}

/**
 * @brief Union-find "find" with path halving.
 *
 * @param parent Parent array of the disjoint-set forest.
 * @param node   Element to look up.
 * @return Representative of the set containing @p node.
 */
static int32_t find_root(int32_t *parent, int32_t node)
{
    while (parent[node] != node) {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    return node;
}

int64_t cover_count_connected_od_pairs(const Graph *graph, const uint8_t *selected)
{
    int32_t *parent = xmalloc_array((size_t)graph->num_nodes, sizeof *parent);
    for (int32_t v = 0; v < graph->num_nodes; v++) {
        parent[v] = v;
    }
    for (int32_t e = 0; e < graph->num_edges; e++) {
        if (!selected[e]) {
            int32_t a = find_root(parent, graph->edge_node_a[e]);
            int32_t b = find_root(parent, graph->edge_node_b[e]);
            if (a != b) {
                parent[a] = b;
            }
        }
    }
    int32_t *od_in_component = xcalloc((size_t)graph->num_nodes, sizeof *od_in_component);
    for (int32_t rank = 0; rank < graph->num_od_nodes; rank++) {
        od_in_component[find_root(parent, graph->od_nodes[rank])]++;
    }
    int64_t connected_pairs = 0;
    for (int32_t v = 0; v < graph->num_nodes; v++) {
        int64_t count = od_in_component[v];
        connected_pairs += count * (count - 1) / 2;
    }
    free(parent);
    free(od_in_component);
    return connected_pairs;
}
