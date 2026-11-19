/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file graph.c
 * @brief Queries on the road network and its destruction.
 *
 * Loading and preprocessing live in graph_load.c.
 */
#include "graph.h"

#include <stdlib.h>

int32_t graph_od_partner_count(const Graph *graph, int32_t rank)
{
    int32_t count = graph->num_od_nodes;
    if (count % 2 == 1) {
        return (count - 1) / 2;
    }
    return rank < count / 2 ? count / 2 : count / 2 - 1;
}

int64_t graph_num_od_pairs(const Graph *graph)
{
    int64_t count = graph->num_od_nodes;
    return count * (count - 1) / 2;
}

void graph_free(Graph *graph)
{
    if (graph == nullptr) {
        return;
    }
    for (int32_t e = 0; e < graph->num_edges; e++) {
        free(graph->edge_label[e]);
    }
    for (int32_t v = 0; v < graph->num_nodes; v++) {
        free(graph->node_label[v]);
    }
    free(graph->edge_node_a);
    free(graph->edge_node_b);
    free(graph->edge_label);
    free(graph->edge_is_locked);
    free(graph->locked_edges);
    free(graph->node_label);
    free(graph->od_rank);
    free(graph->od_nodes);
    free(graph->adjacency_start);
    free(graph->adjacency_edges);
    free(graph);
}
