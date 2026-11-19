/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file separation.c
 * @brief Implementation of the Cut phase (search trees with spread criteria).
 *
 * Differences from the original implementation (same behaviour, faster):
 *  - visited flags use a "generation stamp", so they are not cleared for every root;
 *  - membership of a node in the root's partner list is an O(1) rank test
 *    instead of a scan over the partner list;
 *  - the frontier edge is selected with a linear scan instead of sorting
 *    the whole frontier twice per step;
 *  - paths are stored as flat edge lists instead of persistent tree nodes.
 * One difference in behaviour: the tree stops at O-D nodes, so only direct
 * paths are generated (see separation.h).
 */
#include "separation.h"

#include "util.h"

#include <stdlib.h>
#include <string.h>

struct Separator {
    const Graph *graph;

    uint32_t stamp;          /**< Identifier of the current search tree. */
    uint32_t *visited_stamp; /**< [num_nodes] == stamp iff the node is in the current tree. */
    int32_t *parent_edge;    /**< [num_nodes] Tree edge leading to each visited node. */
    int32_t *depth;          /**< [num_nodes] Number of links from the root. */

    /* Frontier, as three parallel arrays (each edge appears at most once). */
    int32_t *frontier_edge;  /**< [num_edges] Candidate edge.                    */
    int32_t *frontier_node;  /**< [num_edges] Its endpoint outside the tree.     */
    int32_t *frontier_depth; /**< [num_edges] Depth that endpoint would receive. */
    int32_t frontier_size;   /**< Number of candidates.                          */

    int32_t *path_buffer; /**< [num_nodes] Scratch for extracting a tree path. */
};

Separator *separator_create(const Graph *graph)
{
    size_t nodes = (size_t)graph->num_nodes;
    size_t edges = (size_t)graph->num_edges;
    Separator *separator = xcalloc(1, sizeof *separator);
    separator->graph = graph;
    separator->visited_stamp = xcalloc(nodes, sizeof *separator->visited_stamp);
    separator->parent_edge = xmalloc_array(nodes, sizeof *separator->parent_edge);
    separator->depth = xmalloc_array(nodes, sizeof *separator->depth);
    separator->frontier_edge = xmalloc_array(edges, sizeof *separator->frontier_edge);
    separator->frontier_node = xmalloc_array(edges, sizeof *separator->frontier_node);
    separator->frontier_depth = xmalloc_array(edges, sizeof *separator->frontier_depth);
    separator->path_buffer = xmalloc_array(nodes, sizeof *separator->path_buffer);
    return separator;
}

void separator_free(Separator *separator)
{
    if (separator == nullptr) {
        return;
    }
    free(separator->visited_stamp);
    free(separator->parent_edge);
    free(separator->depth);
    free(separator->frontier_edge);
    free(separator->frontier_node);
    free(separator->frontier_depth);
    free(separator->path_buffer);
    free(separator);
}

/**
 * @brief Starts a new, empty search tree (invalidates all visited flags in O(1)).
 *
 * @param separator Scratch memory.
 */
static void start_new_tree(Separator *separator)
{
    separator->stamp++;
    if (separator->stamp == 0) { /* counter wrapped around: clear for real */
        memset(separator->visited_stamp, 0,
               (size_t)separator->graph->num_nodes * sizeof *separator->visited_stamp);
        separator->stamp = 1;
    }
    separator->frontier_size = 0;
}

/**
 * @brief Tells whether @p node belongs to the current search tree.
 *
 * @param separator Scratch memory.
 * @param node      Node index.
 * @return true if visited in the current tree.
 */
static inline bool is_visited(const Separator *separator, int32_t node)
{
    return separator->visited_stamp[node] == separator->stamp;
}

/**
 * @brief Adds to the frontier every usable edge leaving the tree node @p node.
 *
 * An edge is usable when it has no counting station and its other endpoint
 * is not in the tree yet.
 *
 * @param separator Scratch memory.
 * @param cover     Current cover (selection flags).
 * @param node      A node that was just added to the tree.
 */
static void expand_node(Separator *separator, const Cover *cover, int32_t node)
{
    const Graph *graph = separator->graph;
    int32_t next_depth = separator->depth[node] + 1;
    for (int32_t i = graph->adjacency_start[node]; i < graph->adjacency_start[node + 1]; i++) {
        int32_t edge = graph->adjacency_edges[i];
        if (cover->selected[edge]) {
            continue;
        }
        int32_t far = graph_other_endpoint(graph, edge, node);
        if (is_visited(separator, far)) {
            continue;
        }
        int32_t slot = separator->frontier_size++;
        separator->frontier_edge[slot] = edge;
        separator->frontier_node[slot] = far;
        separator->frontier_depth[slot] = next_depth;
    }
}

/**
 * @brief Priority of a frontier candidate under a spread criterion (lower is better).
 *
 * @param separator Scratch memory.
 * @param pool      Path pool (for edge usage counts).
 * @param criterion Spread criterion.
 * @param slot      Frontier position.
 * @return (depth, usage) packed into one integer for SPREAD_MLC, usage alone for SPREAD_RANDOM.
 */
static inline uint64_t candidate_key(const Separator *separator, const PathPool *pool,
                                     SpreadCriterion criterion, int32_t slot)
{
    uint64_t usage = (uint64_t)path_pool_edge_usage(pool, separator->frontier_edge[slot]);
    if (criterion == SPREAD_RANDOM) {
        return usage;
    }
    return ((uint64_t)separator->frontier_depth[slot] << 32) | usage;
}

/**
 * @brief Chooses the next frontier position according to the spread criterion.
 *
 * Inner workings: a first scan finds the best (lowest) key and how many
 * candidates share it; a second scan returns one of those ties chosen
 * uniformly at random.
 *
 * @param separator Scratch memory (frontier must be non-empty).
 * @param pool      Path pool (for edge usage counts).
 * @param criterion Spread criterion.
 * @param rng       Random generator.
 * @return Frontier position of the chosen candidate.
 */
static int32_t choose_candidate(const Separator *separator, const PathPool *pool,
                                SpreadCriterion criterion, Rng *rng)
{
    uint64_t best_key = UINT64_MAX;
    uint32_t ties = 0;
    for (int32_t slot = 0; slot < separator->frontier_size; slot++) {
        uint64_t key = candidate_key(separator, pool, criterion, slot);
        if (key < best_key) {
            best_key = key;
            ties = 1;
        } else if (key == best_key) {
            ties++;
        }
    }
    uint32_t chosen = rng_below(rng, ties);
    for (int32_t slot = 0; slot < separator->frontier_size; slot++) {
        if (candidate_key(separator, pool, criterion, slot) == best_key && chosen-- == 0) {
            return slot;
        }
    }
    return 0; /* unreachable */
}

/**
 * @brief Moves the candidate at @p slot into the tree.
 *
 * Inner workings: records the tree edge and depth of its outer endpoint v,
 * marks v visited, then removes from the frontier every candidate leading to
 * v (they would now join two tree nodes).
 *
 * @param separator Scratch memory.
 * @param slot      Frontier position of the chosen candidate.
 * @return The node that joined the tree.
 */
static int32_t attach_candidate(Separator *separator, int32_t slot)
{
    int32_t node = separator->frontier_node[slot];
    separator->visited_stamp[node] = separator->stamp;
    separator->parent_edge[node] = separator->frontier_edge[slot];
    separator->depth[node] = separator->frontier_depth[slot];

    for (int32_t i = 0; i < separator->frontier_size;) {
        if (separator->frontier_node[i] == node) {
            int32_t last = --separator->frontier_size;
            separator->frontier_edge[i] = separator->frontier_edge[last];
            separator->frontier_node[i] = separator->frontier_node[last];
            separator->frontier_depth[i] = separator->frontier_depth[last];
        } else {
            i++;
        }
    }
    return node;
}

/**
 * @brief Writes the edges of the tree path from @p node back to @p root.
 *
 * @param separator Scratch memory; the edges go to path_buffer.
 * @param node      Tree node where the path ends.
 * @param root      Root of the current tree.
 * @return Number of edges in the path.
 */
static int32_t extract_tree_path(Separator *separator, int32_t node, int32_t root)
{
    int32_t length = 0;
    while (node != root) {
        int32_t edge = separator->parent_edge[node];
        separator->path_buffer[length++] = edge;
        node = graph_other_endpoint(separator->graph, edge, node);
    }
    return length;
}

int32_t separator_run(Separator *separator, const Cover *cover, PathPool *pool,
                      SpreadCriterion criterion, int32_t max_paths_per_root, Rng *rng)
{
    const Graph *graph = separator->graph;
    int32_t new_paths = 0;

    for (int32_t rank = 0; rank < graph->num_od_nodes; rank++) {
        int32_t paths_left = graph_od_partner_count(graph, rank);
        if (max_paths_per_root > 0 && max_paths_per_root < paths_left) {
            paths_left = max_paths_per_root;
        }
        if (paths_left <= 0) {
            continue;
        }

        int32_t root = graph->od_nodes[rank];
        start_new_tree(separator);
        separator->visited_stamp[root] = separator->stamp;
        separator->depth[root] = 0;
        expand_node(separator, cover, root);

        while (separator->frontier_size > 0) {
            int32_t slot = choose_candidate(separator, pool, criterion, rng);
            int32_t node = attach_candidate(separator, slot);
            if (graph->od_rank[node] < 0) { /* direct paths: do not grow past O-D nodes */
                expand_node(separator, cover, node);
            }

            if (graph_is_od_partner(graph, rank, node)) {
                int32_t length = extract_tree_path(separator, node, root);
                path_pool_add(pool, separator->path_buffer, length);
                new_paths++;
                if (--paths_left == 0) {
                    break;
                }
            }
        }
    }
    return new_paths;
}
