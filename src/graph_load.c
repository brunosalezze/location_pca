/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file graph_load.c
 * @brief Reading instance files and preprocessing the road network.
 *
 * The work is done on a mutable "builder" representation (per-node incident
 * lists that can be edited) and the result is then compacted into the
 * immutable CSR Graph used by the algorithms.
 */
#include "graph.h"

#include "strmap.h"
#include "util.h"

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** Mutable network used while reading and preprocessing an instance. */
typedef struct Builder {
    /* nodes */
    StrMap node_index;   /**< label -> node index */
    char **node_label;   /**< [node_capacity] */
    uint8_t *node_is_od; /**< [node_capacity] */
    uint8_t *node_alive; /**< [num_nodes] (allocated after reading) */
    IntVec *incident;    /**< [num_nodes] incident edges (allocated after reading) */
    int32_t num_nodes;
    size_t node_capacity;

    /* edges */
    char **edge_label;   /**< [edge_capacity] */
    int32_t *edge_a;     /**< [edge_capacity] */
    int32_t *edge_b;     /**< [edge_capacity] */
    uint8_t *edge_alive; /**< [num_edges] (allocated after reading) */
    int32_t num_edges;
    size_t edge_capacity;

    IntVec od_list; /**< O-D node indices in file order, without duplicates. */
    GraphStats stats;
} Builder;

/* ------------------------------------------------------------------------- */
/* Text parsing helpers                                                      */
/* ------------------------------------------------------------------------- */

/**
 * @brief Reads a whole text file into a NUL-terminated buffer.
 *
 * A UTF-8 byte order mark at the start of the file (written by some Windows
 * editors) is dropped, otherwise it would become part of the first label.
 *
 * @param path File to read.
 * @return Newly allocated buffer (free() it). Stops the program if the file
 *         cannot be read or contains a NUL byte (e.g. a UTF-16 file), which
 *         would otherwise silently cut the input short.
 */
static char *read_text_file(const char *path)
{
    FILE *file = fopen(path, "rb");
    if (file == nullptr) {
        fatal("cannot open '%s': %s", path, strerror(errno));
    }
    size_t capacity = 1 << 16;
    size_t length = 0;
    char *buffer = xmalloc(capacity);
    for (;;) {
        if (capacity - length < 2) {
            capacity *= 2;
            buffer = xrealloc_array(buffer, capacity, 1);
        }
        size_t wanted = capacity - length - 1;
        size_t read = fread(buffer + length, 1, wanted, file);
        length += read;
        if (read < wanted) { /* end of file or read error: never read past either */
            break;
        }
    }
    if (ferror(file)) {
        fatal("error while reading '%s'", path);
    }
    fclose(file);
    buffer[length] = '\0';
    if (memchr(buffer, '\0', length) != nullptr) {
        fatal("'%s' contains a NUL byte: it is not a plain-text (ASCII or UTF-8) file", path);
    }
    static const char utf8_bom[] = "\xEF\xBB\xBF";
    if (length >= 3 && memcmp(buffer, utf8_bom, 3) == 0) {
        memmove(buffer, buffer + 3, length - 3 + 1); /* + 1: the terminating NUL */
    }
    return buffer;
}

/**
 * @brief Cuts the next line out of a text buffer, in place.
 *
 * @param cursor In: start of the remaining text. Out: start of the following line.
 * @return The line (without its '\n'), or nullptr when the text is exhausted.
 */
static char *next_line(char **cursor)
{
    char *line = *cursor;
    if (*line == '\0') {
        return nullptr;
    }
    char *end = strchr(line, '\n');
    if (end != nullptr) {
        *end = '\0';
        *cursor = end + 1;
    } else {
        *cursor = line + strlen(line);
    }
    return line;
}

/**
 * @brief Cuts the next whitespace-separated token out of a line, in place.
 *
 * Any whitespace (spaces, tabs and the carriage return of Windows line endings)
 * separates tokens.
 *
 * @param cursor In: position in the line. Out: position after the token.
 * @return The token, or nullptr if the line has no more tokens.
 */
static char *next_token(char **cursor)
{
    char *p = *cursor;
    while (*p != '\0' && isspace((unsigned char)*p)) {
        p++;
    }
    if (*p == '\0') {
        *cursor = p;
        return nullptr;
    }
    char *token = p;
    while (*p != '\0' && !isspace((unsigned char)*p)) {
        p++;
    }
    if (*p != '\0') {
        *p++ = '\0';
    }
    *cursor = p;
    return token;
}

/* ------------------------------------------------------------------------- */
/* Builder: reading                                                          */
/* ------------------------------------------------------------------------- */

/**
 * @brief Returns the index of the node labelled @p label, creating it if needed.
 *
 * @param builder Network under construction.
 * @param label   Node label as read from a file.
 * @return Node index.
 */
static int32_t builder_node(Builder *builder, const char *label)
{
    int32_t node = strmap_get(&builder->node_index, label);
    if (node >= 0) {
        return node;
    }
    if ((size_t)builder->num_nodes == builder->node_capacity) {
        builder->node_capacity = builder->node_capacity > 0 ? 2 * builder->node_capacity : 256;
        builder->node_label = xrealloc_array(builder->node_label, builder->node_capacity,
                                             sizeof *builder->node_label);
        builder->node_is_od = xrealloc_array(builder->node_is_od, builder->node_capacity,
                                             sizeof *builder->node_is_od);
    }
    node = builder->num_nodes++;
    builder->node_label[node] = xstrdup(label);
    builder->node_is_od[node] = 0;
    strmap_put(&builder->node_index, builder->node_label[node], node);
    return node;
}

/**
 * @brief Appends an edge between two existing nodes.
 *
 * @param builder Network under construction.
 * @param label   Edge label (copied).
 * @param a       First endpoint.
 * @param b       Second endpoint.
 */
static void builder_add_edge(Builder *builder, const char *label, int32_t a, int32_t b)
{
    if ((size_t)builder->num_edges == builder->edge_capacity) {
        builder->edge_capacity = builder->edge_capacity > 0 ? 2 * builder->edge_capacity : 256;
        builder->edge_label = xrealloc_array(builder->edge_label, builder->edge_capacity,
                                             sizeof *builder->edge_label);
        builder->edge_a =
            xrealloc_array(builder->edge_a, builder->edge_capacity, sizeof *builder->edge_a);
        builder->edge_b =
            xrealloc_array(builder->edge_b, builder->edge_capacity, sizeof *builder->edge_b);
    }
    int32_t edge = builder->num_edges++;
    builder->edge_label[edge] = xstrdup(label);
    builder->edge_a[edge] = a;
    builder->edge_b[edge] = b;
}

/**
 * @brief Parses the edge file: one `edge_id node_a node_b [lanes]` per line.
 *
 * Self-loops are counted and skipped right away (preprocessing step 1).
 *
 * @param builder Network under construction.
 * @param path    Edge file path (used for reading and in error messages).
 */
static void read_edge_file(Builder *builder, const char *path)
{
    char *text = read_text_file(path);
    char *cursor = text;
    int line_number = 0;
    for (char *line; (line = next_line(&cursor)) != nullptr;) {
        line_number++;
        char *field = line;
        const char *edge_label = next_token(&field);
        if (edge_label == nullptr || edge_label[0] == '#') {
            continue; /* blank or comment line */
        }
        const char *label_a = next_token(&field);
        const char *label_b = next_token(&field);
        if (label_b == nullptr) {
            fatal("%s:%d: expected 'edge_id node_a node_b [lanes]'", path, line_number);
        }
        builder->stats.raw_edges++;
        int32_t a = builder_node(builder, label_a);
        int32_t b = builder_node(builder, label_b);
        if (a == b) {
            builder->stats.self_loops_removed++;
            continue;
        }
        builder_add_edge(builder, edge_label, a, b);
    }
    free(text);
    if (builder->stats.raw_edges == 0) {
        fatal("'%s' contains no edges", path);
    }
}

/**
 * @brief Parses the O-D file: one node label per line.
 *
 * Labels that do not occur in the edge file create isolated nodes (such an
 * O-D node trivially needs no counting station). Repeated labels are counted
 * and ignored: the original code turned them into an O-D "pair" made of a
 * single node, which can never be separated.
 *
 * @param builder Network under construction.
 * @param path    O-D file path.
 */
static void read_od_file(Builder *builder, const char *path)
{
    char *text = read_text_file(path);
    char *cursor = text;
    for (char *line; (line = next_line(&cursor)) != nullptr;) {
        char *field = line;
        const char *label = next_token(&field);
        if (label == nullptr || label[0] == '#') {
            continue;
        }
        for (; label != nullptr; label = next_token(&field)) {
            int32_t node = builder_node(builder, label);
            if (builder->node_is_od[node]) {
                builder->stats.duplicate_od_entries++;
                continue;
            }
            builder->node_is_od[node] = 1;
            intvec_push(&builder->od_list, node);
        }
    }
    free(text);
    if (builder->od_list.count == 0) {
        fatal("'%s' contains no O-D nodes", path);
    }
}

/* ------------------------------------------------------------------------- */
/* Builder: preprocessing                                                    */
/* ------------------------------------------------------------------------- */

/** Sort key used to find parallel edges. */
typedef struct EdgeKey {
    int32_t low;  /**< Smaller endpoint index. */
    int32_t high; /**< Larger endpoint index.  */
    int32_t edge; /**< Edge index (tie-breaker: keeps file order). */
} EdgeKey;

/**
 * @brief qsort() comparator ordering EdgeKey by (low, high, edge).
 *
 * @param lhs Pointer to an EdgeKey.
 * @param rhs Pointer to an EdgeKey.
 * @return Negative, zero or positive, like strcmp().
 */
static int compare_edge_keys(const void *lhs, const void *rhs)
{
    const EdgeKey *a = lhs;
    const EdgeKey *b = rhs;
    if (a->low != b->low) {
        return a->low < b->low ? -1 : 1;
    }
    if (a->high != b->high) {
        return a->high < b->high ? -1 : 1;
    }
    return (a->edge > b->edge) - (a->edge < b->edge);
}

/**
 * @brief Preprocessing step 2: removes parallel duplicate edges.
 *
 * Inner workings: sorts the edges by their (unordered) pair of endpoints, so
 * duplicates become adjacent; within a group the edge that appears first in
 * the file is kept. O(|A| log |A|), instead of the original O(|A|^2) scan.
 *
 * @param builder Network under construction (edge_alive must be allocated).
 */
static void remove_parallel_edges(Builder *builder)
{
    EdgeKey *keys = xmalloc_array((size_t)builder->num_edges, sizeof *keys);
    for (int32_t e = 0; e < builder->num_edges; e++) {
        int32_t a = builder->edge_a[e];
        int32_t b = builder->edge_b[e];
        keys[e] = (EdgeKey){.low = a < b ? a : b, .high = a < b ? b : a, .edge = e};
    }
    qsort(keys, (size_t)builder->num_edges, sizeof *keys, compare_edge_keys);
    for (int32_t i = 1; i < builder->num_edges; i++) {
        if (keys[i].low == keys[i - 1].low && keys[i].high == keys[i - 1].high) {
            builder->edge_alive[keys[i].edge] = 0;
            builder->stats.parallel_edges_removed++;
        }
    }
    free(keys);
}

/**
 * @brief Removes the first occurrence of @p value from @p vec, keeping the order.
 *
 * @param vec   Vector to modify.
 * @param value Value to remove (must be present).
 */
static void intvec_remove_value(IntVec *vec, int32_t value)
{
    for (size_t i = 0; i < vec->count; i++) {
        if (vec->items[i] == value) {
            memmove(&vec->items[i], &vec->items[i + 1], (vec->count - i - 1) * sizeof *vec->items);
            vec->count--;
            return;
        }
    }
}

/**
 * @brief Replaces the first occurrence of @p old_value in @p vec by @p new_value.
 *
 * @param vec       Vector to modify.
 * @param old_value Value to look for (must be present).
 * @param new_value Replacement.
 */
static void intvec_replace_value(IntVec *vec, int32_t old_value, int32_t new_value)
{
    for (size_t i = 0; i < vec->count; i++) {
        if (vec->items[i] == old_value) {
            vec->items[i] = new_value;
            return;
        }
    }
}

/**
 * @brief Tells whether @p node is an intersection that step 3 can contract.
 *
 * @param builder Network under construction.
 * @param node    Node index.
 * @return true for a live, non-O-D node with exactly two incident edges.
 */
static bool is_contractible(const Builder *builder, int32_t node)
{
    return builder->node_alive[node] && !builder->node_is_od[node] &&
           builder->incident[node].count == 2;
}

/**
 * @brief Returns the endpoint of @p edge opposite to @p node in the builder.
 *
 * @param builder Network under construction.
 * @param edge    Edge index.
 * @param node    One endpoint of the edge.
 * @return The other endpoint.
 */
static int32_t builder_other_endpoint(const Builder *builder, int32_t edge, int32_t node)
{
    return builder->edge_a[edge] != node ? builder->edge_a[edge] : builder->edge_b[edge];
}

/**
 * @brief Preprocessing step 3: contracts every degree-2 intersection.
 *
 * Inner workings: for a contractible node v with incident edges e0 = {x, v}
 * and e1 = {v, y}:
 *  - if x != y, edge e0 is re-attached to become {x, y} (keeping e0's label),
 *    e1 is deleted and v disappears. Degrees of x and y do not change;
 *  - if x == y, v closes a dead-end ring hanging from x: both edges and v are
 *    deleted, which lowers the degree of x by two, so x is re-examined. (The
 *    original code turned such rings into self-loops instead. Harmless in the
 *    instances, but a self-loop on an O-D node would have been "locked" into
 *    every solution.)
 * A stack of candidate nodes makes the whole step linear in the graph size.
 *
 * @param builder Network under construction (incident lists must be built).
 */
static void contract_degree_two_nodes(Builder *builder)
{
    IntVec pending = {0};
    for (int32_t v = builder->num_nodes - 1; v >= 0; v--) {
        if (is_contractible(builder, v)) {
            intvec_push(&pending, v);
        }
    }
    while (pending.count > 0) {
        int32_t v = pending.items[--pending.count];
        if (!is_contractible(builder, v)) {
            continue;
        }
        int32_t kept = builder->incident[v].items[0];
        int32_t dropped = builder->incident[v].items[1];
        int32_t x = builder_other_endpoint(builder, kept, v);
        int32_t y = builder_other_endpoint(builder, dropped, v);

        builder->node_alive[v] = 0;
        builder->incident[v].count = 0;
        builder->edge_alive[dropped] = 0;
        builder->stats.nodes_contracted++;

        if (x == y) {
            intvec_remove_value(&builder->incident[x], kept);
            intvec_remove_value(&builder->incident[x], dropped);
            builder->edge_alive[kept] = 0;
            builder->stats.ring_edges_removed += 2;
            if (is_contractible(builder, x)) {
                intvec_push(&pending, x);
            }
        } else {
            if (builder->edge_a[kept] == v) {
                builder->edge_a[kept] = y;
            } else {
                builder->edge_b[kept] = y;
            }
            intvec_replace_value(&builder->incident[y], dropped, kept);
        }
    }
    intvec_free(&pending);
}

/* ------------------------------------------------------------------------- */
/* Builder -> Graph                                                          */
/* ------------------------------------------------------------------------- */

/**
 * @brief Converts the preprocessed builder into a compact, immutable Graph.
 *
 * Inner workings: live nodes and edges are renumbered densely, preserving
 * their relative order. Labels are moved (not copied) into the graph. The
 * CSR adjacency, the O-D ranks and the list of locked edges are built here.
 *
 * @param builder Preprocessed network. Its label arrays are emptied.
 * @return Newly allocated graph.
 */
static Graph *build_graph(Builder *builder)
{
    Graph *graph = xcalloc(1, sizeof *graph);
    graph->stats = builder->stats;

    int32_t *new_node = xmalloc_array((size_t)builder->num_nodes, sizeof *new_node);
    int32_t *new_edge = xmalloc_array((size_t)builder->num_edges, sizeof *new_edge);
    for (int32_t v = 0; v < builder->num_nodes; v++) {
        new_node[v] = builder->node_alive[v] ? graph->num_nodes++ : -1;
    }
    for (int32_t e = 0; e < builder->num_edges; e++) {
        new_edge[e] = builder->edge_alive[e] ? graph->num_edges++ : -1;
    }

    /* Edges. */
    size_t num_edges = (size_t)graph->num_edges;
    graph->edge_node_a = xmalloc_array(num_edges, sizeof *graph->edge_node_a);
    graph->edge_node_b = xmalloc_array(num_edges, sizeof *graph->edge_node_b);
    graph->edge_label = xmalloc_array(num_edges, sizeof *graph->edge_label);
    for (int32_t e = 0; e < builder->num_edges; e++) {
        if (new_edge[e] < 0) {
            free(builder->edge_label[e]);
        } else {
            graph->edge_node_a[new_edge[e]] = new_node[builder->edge_a[e]];
            graph->edge_node_b[new_edge[e]] = new_node[builder->edge_b[e]];
            graph->edge_label[new_edge[e]] = builder->edge_label[e];
        }
        builder->edge_label[e] = nullptr;
    }

    /* Nodes and O-D ranks. */
    size_t num_nodes = (size_t)graph->num_nodes;
    graph->node_label = xmalloc_array(num_nodes, sizeof *graph->node_label);
    graph->od_rank = xmalloc_array(num_nodes, sizeof *graph->od_rank);
    for (int32_t v = 0; v < builder->num_nodes; v++) {
        if (new_node[v] < 0) {
            free(builder->node_label[v]);
        } else {
            graph->node_label[new_node[v]] = builder->node_label[v];
            graph->od_rank[new_node[v]] = -1;
        }
        builder->node_label[v] = nullptr;
    }
    graph->num_od_nodes = (int32_t)builder->od_list.count;
    graph->od_nodes = xmalloc_array(builder->od_list.count, sizeof *graph->od_nodes);
    for (int32_t rank = 0; rank < graph->num_od_nodes; rank++) {
        int32_t node = new_node[builder->od_list.items[rank]];
        graph->od_nodes[rank] = node;
        graph->od_rank[node] = rank;
    }

    /* CSR adjacency (keeps each node's incident-edge order). */
    graph->adjacency_start = xmalloc_array(num_nodes + 1, sizeof *graph->adjacency_start);
    graph->adjacency_edges = xmalloc_array(2 * num_edges, sizeof *graph->adjacency_edges);
    int32_t offset = 0;
    for (int32_t v = 0; v < builder->num_nodes; v++) {
        if (new_node[v] < 0) {
            continue;
        }
        graph->adjacency_start[new_node[v]] = offset;
        const IntVec *incident = &builder->incident[v];
        for (size_t i = 0; i < incident->count; i++) {
            graph->adjacency_edges[offset++] = new_edge[incident->items[i]];
        }
        if (incident->count == 0 && builder->node_is_od[v]) {
            graph->stats.isolated_od_nodes++;
        }
    }
    graph->adjacency_start[num_nodes] = offset;

    /* Locked edges: both endpoints are O-D nodes. */
    graph->edge_is_locked = xcalloc(num_edges, sizeof *graph->edge_is_locked);
    graph->locked_edges = xmalloc_array(num_edges, sizeof *graph->locked_edges);
    for (int32_t e = 0; e < graph->num_edges; e++) {
        if (graph->od_rank[graph->edge_node_a[e]] >= 0 &&
            graph->od_rank[graph->edge_node_b[e]] >= 0) {
            graph->edge_is_locked[e] = 1;
            graph->locked_edges[graph->num_locked_edges++] = e;
        }
    }

    free(new_node);
    free(new_edge);
    return graph;
}

/**
 * @brief Releases everything the builder still owns.
 *
 * @param builder Builder to release.
 */
static void builder_free(Builder *builder)
{
    for (int32_t v = 0; v < builder->num_nodes; v++) {
        free(builder->node_label[v]); /* nullptr after build_graph() */
        if (builder->incident != nullptr) {
            intvec_free(&builder->incident[v]);
        }
    }
    for (int32_t e = 0; e < builder->num_edges; e++) {
        free(builder->edge_label[e]);
    }
    strmap_free(&builder->node_index);
    free(builder->node_label);
    free(builder->node_is_od);
    free(builder->node_alive);
    free(builder->incident);
    free(builder->edge_label);
    free(builder->edge_a);
    free(builder->edge_b);
    free(builder->edge_alive);
    intvec_free(&builder->od_list);
}

Graph *graph_load(const char *edges_path, const char *od_path, bool contract_degree_two)
{
    Builder builder = {0};
    strmap_init(&builder.node_index, 1024);

    read_edge_file(&builder, edges_path);
    read_od_file(&builder, od_path);

    /* Every node and edge starts alive. */
    builder.node_alive = xmalloc_array((size_t)builder.num_nodes, 1);
    memset(builder.node_alive, 1, (size_t)builder.num_nodes);
    builder.edge_alive = xmalloc_array((size_t)builder.num_edges, 1);
    memset(builder.edge_alive, 1, (size_t)builder.num_edges);

    remove_parallel_edges(&builder);

    builder.incident = xcalloc((size_t)builder.num_nodes, sizeof *builder.incident);
    for (int32_t e = 0; e < builder.num_edges; e++) {
        if (builder.edge_alive[e]) {
            intvec_push(&builder.incident[builder.edge_a[e]], e);
            intvec_push(&builder.incident[builder.edge_b[e]], e);
        }
    }

    if (contract_degree_two) {
        contract_degree_two_nodes(&builder);
    }

    Graph *graph = build_graph(&builder);
    builder_free(&builder);
    return graph;
}
