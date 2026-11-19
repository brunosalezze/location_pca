/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file graph.h
 * @brief Road network G = (N, A) with its set T of origin-destination (O-D) nodes.
 *
 * Notation follows Section 3 of the paper: N is the set of nodes, A the set of
 * (bidirectional) road links, called *edges* in the code, and T ⊆ N the set of O-D
 * nodes (municipalities). Every unordered pair of distinct O-D nodes is an O-D
 * pair w ∈ W that must be separated by the counting stations.
 *
 * Nodes and edges are identified by dense indices (0 .. num_nodes-1 and
 * 0 .. num_edges-1). The labels from the input files are kept only to write
 * solutions. The adjacency is stored in CSR form: the edges incident to node v
 * are adjacency_edges[adjacency_start[v] .. adjacency_start[v+1]-1].
 *
 * A graph is created by graph_load(), which also preprocesses the network
 * (see graph_load() for the exact rules), and released with graph_free().
 */
#ifndef LOCATION_PCA_GRAPH_H
#define LOCATION_PCA_GRAPH_H

#include <stdint.h>

/** Counters describing what the preprocessing in graph_load() changed. */
typedef struct GraphStats {
    int32_t raw_edges;              /**< Edge lines read from the edge file.             */
    int32_t self_loops_removed;     /**< Edges whose two endpoints are the same node.    */
    int32_t parallel_edges_removed; /**< Duplicate links joining an already-linked pair. */
    int32_t nodes_contracted;       /**< Degree-2 intersections merged into one link.    */
    int32_t ring_edges_removed;     /**< Links of dead-end rings deleted by contraction. */
    int32_t duplicate_od_entries;   /**< Repeated lines in the O-D file (ignored).       */
    int32_t isolated_od_nodes;      /**< O-D nodes that touch no edge at all.            */
} GraphStats;

/** Road network plus O-D nodes. All arrays are owned by the structure. */
typedef struct Graph {
    int32_t num_nodes;    /**< |N| after preprocessing. */
    int32_t num_edges;    /**< |A| after preprocessing. */
    int32_t num_od_nodes; /**< |T|.                     */

    /* --- edges ----------------------------------------------------------- */
    int32_t *edge_node_a; /**< [num_edges] First endpoint of each edge.  */
    int32_t *edge_node_b; /**< [num_edges] Second endpoint of each edge. */
    char **edge_label;    /**< [num_edges] Edge identifier from the input file. */

    /**
     * [num_edges] 1 if both endpoints are O-D nodes. Such an edge is by itself a
     * path between an O-D pair, so every feasible cover must contain it.
     */
    uint8_t *edge_is_locked;
    int32_t *locked_edges;    /**< [num_locked_edges] Indices of the locked edges. */
    int32_t num_locked_edges; /**< Number of locked edges. */

    /* --- nodes ----------------------------------------------------------- */
    char **node_label; /**< [num_nodes] Node identifier from the input files. */
    int32_t *od_rank;  /**< [num_nodes] Position of the node in the O-D list, or -1. */
    int32_t *od_nodes; /**< [num_od_nodes] Node index of each O-D node, in file order. */

    /* --- CSR adjacency ---------------------------------------------------- */
    int32_t *adjacency_start; /**< [num_nodes + 1] Offsets into adjacency_edges. */
    int32_t *adjacency_edges; /**< [2 * num_edges] Incident edge indices.       */

    GraphStats stats; /**< What the preprocessing did (for reporting). */
} Graph;

/**
 * @brief Reads an instance from disk and preprocesses it.
 *
 * Input files (plain text, whitespace separated, `#` starts a comment line;
 * Windows or Unix line endings, a UTF-8 byte order mark and a missing final
 * newline are all fine):
 *  - @p edges_path: one road link per line, `edge_id node_a node_b [lanes]`.
 *    Extra columns such as the number of lanes are ignored.
 *  - @p od_path: one O-D node label per line. Duplicate lines are ignored.
 *
 * Preprocessing, in order (each step keeps the optimal number of counting
 * stations unchanged, except step 2 which is a modelling decision of the
 * paper's instances):
 *  1. Self-loops are dropped: they cannot lie on a simple O-D path.
 *  2. Parallel duplicates (several links joining the same two nodes in the
 *     raw data, e.g. the two carriageways of one road) are merged. The first
 *     link in file order is kept.
 *  3. If @p contract_degree_two is true, every intersection (non-O-D node) of
 *     degree 2 is removed and its two links are merged into one: a path through
 *     one of them always uses the other, so they are interchangeable for
 *     covering. When this closes a dead-end ring (both links lead back to the
 *     same node) the ring is deleted, since no simple O-D path can use it.
 *
 * @param edges_path          Path of the edge file.
 * @param od_path             Path of the O-D node file.
 * @param contract_degree_two Apply step 3 (the paper's setting is true).
 * @return Newly allocated graph; release it with graph_free(). Malformed or
 *         unreadable files stop the program with an error message.
 */
[[nodiscard]] Graph *graph_load(const char *edges_path, const char *od_path,
                                bool contract_degree_two);

/**
 * @brief Releases a graph and all its arrays.
 *
 * @param graph Graph returned by graph_load(), or nullptr.
 */
void graph_free(Graph *graph);

/**
 * @brief Returns the endpoint of @p edge that is not @p node.
 *
 * @param graph The graph.
 * @param edge  Edge index.
 * @param node  One endpoint of @p edge.
 * @return The opposite endpoint.
 */
static inline int32_t graph_other_endpoint(const Graph *graph, int32_t edge, int32_t node)
{
    int32_t a = graph->edge_node_a[edge];
    return a != node ? a : graph->edge_node_b[edge];
}

/**
 * @brief Number of O-D partners whose paths are searched from O-D node @p rank.
 *
 * Each unordered O-D pair {i, j} is assigned to exactly one of its two nodes,
 * which becomes its *root*. Instead of the textbook rule "i < j", which
 * makes the first O-D nodes do almost all the work, pairs are assigned
 * cyclically: root i is responsible for the next k(i) nodes in the O-D list,
 * i+1, ..., i+k(i) (modulo |T|), with
 *  - k(i) = (|T|-1)/2                         if |T| is odd;
 *  - k(i) = |T|/2 for i < |T|/2, |T|/2 - 1 otherwise, if |T| is even.
 * Every root therefore searches for about |T|/2 partners and every pair is
 * covered exactly once (same rule as the original implementation).
 *
 * @param graph The graph.
 * @param rank  Position of the root in the O-D list (0 .. num_od_nodes-1).
 * @return k(rank).
 */
int32_t graph_od_partner_count(const Graph *graph, int32_t rank);

/**
 * @brief Tells whether @p node is one of the O-D partners assigned to root @p root_rank.
 *
 * Constant time: computes the cyclic distance between the two O-D ranks and
 * compares it with graph_od_partner_count().
 *
 * @param graph     The graph.
 * @param root_rank O-D rank of the root.
 * @param node      Any node index.
 * @return true if {root, node} is an O-D pair assigned to the root.
 */
static inline bool graph_is_od_partner(const Graph *graph, int32_t root_rank, int32_t node)
{
    int32_t rank = graph->od_rank[node];
    if (rank < 0) {
        return false;
    }
    int32_t distance = rank - root_rank;
    if (distance < 0) {
        distance += graph->num_od_nodes;
    }
    return distance >= 1 && distance <= graph_od_partner_count(graph, root_rank);
}

/**
 * @brief Number of O-D pairs |W| = |T| (|T| - 1) / 2.
 *
 * @param graph The graph.
 * @return The number of unordered pairs of distinct O-D nodes.
 */
int64_t graph_num_od_pairs(const Graph *graph);

#endif /* LOCATION_PCA_GRAPH_H */
