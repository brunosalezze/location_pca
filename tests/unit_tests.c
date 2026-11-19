/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file unit_tests.c
 * @brief Unit tests for the solver-independent core, built by `make test`.
 *
 * The tests run on small hand-written networks created in a temporary
 * directory, so no Gurobi and no repository instances are required (one
 * optional check reads instances/AC when it is available). exact_cover() is
 * provided by exact_bruteforce.c. The exact method is cross-checked against a
 * separate brute force that enumerates all edge subsets and tests feasibility
 * with cover_count_connected_od_pairs(), an independent union-find check.
 *
 * Everything runs with AddressSanitizer/UBSan and the project's warning set;
 * a failing check prints its location and the program exits with EXIT_FAILURE.
 */
#define _POSIX_C_SOURCE 200809L /* mkdtemp() */

#include "algorithms.h"
#include "cleaning.h"
#include "cover.h"
#include "exact.h"
#include "graph.h"
#include "greedy.h"
#include "path_pool.h"
#include "rng.h"
#include "separation.h"
#include "timer.h"
#include "util.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* The solver handle used whenever a method wants one: the brute-force
 * stand-in from exact_bruteforce.c. */
static ExactSolver *g_solver;

/* ------------------------------------------------------------------------- */
/* Check bookkeeping                                                         */
/* ------------------------------------------------------------------------- */

static int g_checks;
static int g_failures;

static void check_true(bool ok, const char *what, const char *file, int line)
{
    g_checks++;
    if (!ok) {
        g_failures++;
        fprintf(stderr, "FAILED %s:%d: %s\n", file, line, what);
    }
}

static void check_eq(long long got, long long want, const char *what, const char *file, int line)
{
    g_checks++;
    if (got != want) {
        g_failures++;
        fprintf(stderr, "FAILED %s:%d: %s (got %lld, want %lld)\n", file, line, what, got, want);
    }
}

#define CHECK(cond) check_true((cond), #cond, __FILE__, __LINE__)
#define CHECK_MSG(cond, msg) check_true((cond), (msg), __FILE__, __LINE__)
#define CHECK_EQ(got, want, msg) check_eq((long long)(got), (long long)(want), (msg), __FILE__, __LINE__)

/* ------------------------------------------------------------------------- */
/* Temporary instance files                                                  */
/* ------------------------------------------------------------------------- */

static char g_tmpdir[4096];

typedef struct FileNode {
    struct FileNode *next;
    char *name;
} FileNode;

static FileNode *g_files;

/** Records @p name so cleanup_tmpdir() can delete the file later. */
static void register_tmp_file(const char *name)
{
    FileNode *node = xmalloc(sizeof *node);
    node->name = xstrdup(name);
    node->next = g_files;
    g_files = node;
}

/** Newly allocated "<tmpdir>/<name>"; free() after use. */
static char *tmp_path(const char *name)
{
    size_t dir_len = strlen(g_tmpdir);
    size_t name_len = strlen(name);
    char *path = xmalloc_array(dir_len + name_len + 2, 1);
    memcpy(path, g_tmpdir, dir_len);
    path[dir_len] = '/';
    memcpy(path + dir_len + 1, name, name_len + 1);
    return path;
}

static void write_text_file(const char *name, const char *content)
{
    char *path = tmp_path(name);
    FILE *file = fopen(path, "wb");
    if (file == nullptr) {
        fprintf(stderr, "unit_tests: cannot write '%s'\n", path);
        exit(EXIT_FAILURE);
    }
    fwrite(content, 1, strlen(content), file);
    fclose(file);
    free(path);
    register_tmp_file(name);
}

/** Writes both instance files and loads the graph they describe. */
static Graph *write_and_load(const char *stem, const char *edges_text, const char *od_text,
                             bool contract_degree_two)
{
    char edges_name[512];
    char od_name[512];
    snprintf(edges_name, sizeof edges_name, "%s.edges.txt", stem);
    snprintf(od_name, sizeof od_name, "%s.od.txt", stem);
    char *edges_path = tmp_path(edges_name);
    char *od_path = tmp_path(od_name);
    write_text_file(edges_name, edges_text);
    write_text_file(od_name, od_text);
    Graph *graph = graph_load(edges_path, od_path, contract_degree_two);
    free(edges_path);
    free(od_path);
    return graph;
}

static void cleanup_tmpdir(void)
{
    for (FileNode *node = g_files; node != nullptr;) {
        char *path = tmp_path(node->name);
        remove(path);
        free(path);
        FileNode *next = node->next;
        free(node->name);
        free(node);
        node = next;
    }
    g_files = nullptr;
    rmdir(g_tmpdir);
}

/* ------------------------------------------------------------------------- */
/* Test networks                                                             */
/* ------------------------------------------------------------------------- */

/* Diamond: two O-D pairs {O1,O2} and {O2,O3} joined by two edge-disjoint
 * routes each (through m1/m2 and through m3/m4). No locked edges; the optimal
 * TCLP cover has 4 edges. Edge index = order of appearance: e1..e8 -> 0..7. */
#define DIAMOND_EDGES                                                                                     \
    "# diamond network\n"                                                                                 \
    "e1\tO1\tm1\t2\r\n"                                                                                  \
    "e2 O1 m2\n"                                                                                          \
    "e3 m1 O2 1\n"                                                                                        \
    "e4 m2 O2\n"                                                                                          \
    "e5 O2 m3\n"                                                                                          \
    "e6 m3 O3\n"                                                                                          \
    "e7 O2 m4\n"                                                                                          \
    "e8 m4 O3" /* no final newline */

static const char *const diamond_od = "O1\nO2\nO3\n";

/* Path A - x - B - y - C with A, B, C O-D nodes: the two intersections get
 * contracted, leaving the two locked edges l1 = {A,B} and l3 = {B,C}. */
#define LINE_EDGES "l1 A x\nl2 x B\nl3 B y\nl4 y C\n"
static const char *const line_od = "A\nB\nC\n";

/* Path A - c - B with a dead-end ring c - d1 - d2 - c hanging from c.
 * Contracting d1, d2 deletes the ring, contracting c then deletes l2. */
#define RING_EDGES "r1 A c\nr2 c B\nr3 c d1\nr4 d1 d2\nr5 d2 c\n"
static const char *const ring_od = "A\nB\n";

/* Parsing torture file: UTF-8 BOM, comment lines, tab and CRLF separators,
 * an extra lanes column, a parallel duplicate and a missing final newline.
 * The O-D file repeats a label and adds two labels that never occur in an
 * edge (isolated O-D nodes). */
#define PARSING_EDGES                                     \
    "\xEF\xBB\xBF# edges\r\n"                             \
    "e1 OD1 OD2 2\r\n"                                    \
    "\te2\tOD1\tX\t1\n"                                   \
    "# comment line\n"                                    \
    "e3 OD1 X\n"                                          \
    "e4 X Y\n"                                            \
    "e5 Y OD2" /* no final newline */
static const char *const parsing_od = "OD1\nOD2\nOD1\nOD9\nOD8";

/* ------------------------------------------------------------------------- */
/* Small helpers                                                             */
/* ------------------------------------------------------------------------- */

/** Adds the four diamond block paths p0..p3 and the long path p4. */
static void add_diamond_pool_paths(PathPool *pool)
{
    const int32_t p0[] = {0, 2};
    const int32_t p1[] = {1, 3};
    const int32_t p2[] = {4, 5};
    const int32_t p3[] = {6, 7};
    const int32_t p4[] = {0, 2, 4, 5};
    path_pool_add(pool, p0, 2);
    path_pool_add(pool, p1, 2);
    path_pool_add(pool, p2, 2);
    path_pool_add(pool, p3, 2);
    path_pool_add(pool, p4, 4);
}

static void select_edges(Cover *cover, const int32_t *edges, int32_t count)
{
    for (int32_t i = 0; i < count; i++) {
        cover_add_edge(cover, edges[i]);
    }
}

/**
 * @brief True optimum of the TCLP on @p graph, found by enumerating every
 *        edge subset and testing it with the union-find feasibility check.
 */
static int32_t brute_force_tclp_optimum(const Graph *graph)
{
    const int32_t n = graph->num_edges;
    CHECK_MSG(n <= 20, "brute-force optimum is only meant for small test graphs");
    uint8_t selected[20];
    int32_t best = n;
    const uint64_t limit = (uint64_t)1 << (uint64_t)n;
    for (uint64_t mask = 0; mask < limit; mask++) {
        int32_t bits = 0;
        for (int32_t e = 0; e < n; e++) {
            selected[e] = (uint8_t)((mask >> (uint64_t)e) & 1u);
            bits += selected[e];
        }
        if (bits >= best) {
            continue;
        }
        if (cover_count_connected_od_pairs(graph, selected) == 0) {
            best = bits;
        }
    }
    return best;
}

/** Number of unordered O-D pairs assigned to some root (must be |W|). */
static int32_t count_assigned_od_pairs(const Graph *graph)
{
    int32_t count = 0;
    for (int32_t r1 = 0; r1 < graph->num_od_nodes; r1++) {
        for (int32_t r2 = 0; r2 < graph->num_od_nodes; r2++) {
            if (r1 != r2 && graph_is_od_partner(graph, r1, graph->od_nodes[r2])) {
                count++;
            }
        }
    }
    return count;
}

static int32_t sum_od_partner_counts(const Graph *graph)
{
    int32_t sum = 0;
    for (int32_t r = 0; r < graph->num_od_nodes; r++) {
        sum += graph_od_partner_count(graph, r);
    }
    return sum;
}

/* ------------------------------------------------------------------------- */
/* util.h                                                                    */
/* ------------------------------------------------------------------------- */

static void test_util(void)
{
    log_set_level(LOG_DEBUG);
    CHECK_MSG(log_get_level() == LOG_DEBUG, "log_set_level round-trips");
    log_set_level(LOG_QUIET);

    IntVec vec = {0};
    for (int32_t i = 0; i < 5; i++) {
        intvec_push(&vec, 3 * i);
    }
    CHECK_EQ(vec.count, 5, "intvec_push appends");
    CHECK_MSG(vec.items[3] == 9, "intvec_push stores values");
    intvec_reserve(&vec, 64);
    CHECK_MSG(vec.capacity >= 64, "intvec_reserve grows capacity");
    CHECK_EQ(vec.count, 5, "intvec_reserve keeps the values");
    CHECK_MSG(vec.items[4] == 12, "intvec_reserve keeps the contents");
    intvec_free(&vec);
    CHECK_MSG(vec.count == 0 && vec.items == nullptr, "intvec_free resets");

    char *copy = xstrdup("location_pca");
    CHECK_MSG(strcmp(copy, "location_pca") == 0, "xstrdup copies the text");
    CHECK_MSG(copy != nullptr, "xstrdup allocates");
    free(copy);
    free(xmalloc(0)); /* documented as a valid, freeable pointer */
    int *array = xrealloc_array(nullptr, 4, sizeof *array);
    for (int i = 0; i < 4; i++) {
        array[i] = i;
    }
    CHECK_MSG(array[3] == 3, "xrealloc_array yields a usable array");
    free(array);

    double t0 = timer_now();
    volatile double sink = 0.0;
    for (int i = 0; i < 10000; i++) {
        sink += i;
    }
    CHECK_MSG(timer_now() >= t0, "timer_now is monotone");
}

/* ------------------------------------------------------------------------- */
/* rng.h                                                                     */
/* ------------------------------------------------------------------------- */

static void test_rng(void)
{
    Rng a, b, c;
    rng_seed(&a, 12345);
    rng_seed(&b, 12345);
    rng_seed(&c, 987654321);
    bool same = true, different = false;
    for (int i = 0; i < 16; i++) {
        uint32_t va = rng_next(&a), vb = rng_next(&b), vc = rng_next(&c);
        if (va != vb) {
            same = false;
        }
        if (va != vc) {
            different = true;
        }
    }
    CHECK_MSG(same, "equal seeds yield equal sequences");
    CHECK_MSG(different, "different seeds yield different sequences");

    rng_seed(&a, 777);
    bool seen[7] = {0};
    for (int i = 0; i < 5000; i++) {
        uint32_t v = rng_below(&a, 7);
        bool in_range = v < 7;
        CHECK_MSG(in_range, "rng_below stays below the bound");
        if (!in_range) {
            break;
        }
        seen[v] = true;
    }
    int seen_count = 0;
    for (int i = 0; i < 7; i++) {
        seen_count += seen[i];
    }
    CHECK_EQ(seen_count, 7, "rng_below hits every value of its range");
    for (int i = 0; i < 32; i++) {
        CHECK_EQ(rng_below(&a, 1), 0, "rng_below(x, 1) is always 0");
    }
}

/* ------------------------------------------------------------------------- */
/* path_pool.h                                                               */
/* ------------------------------------------------------------------------- */

static void test_path_pool(void)
{
    PathPool *pool = path_pool_create(4);
    const int32_t p0[] = {0, 1};
    const int32_t p1[] = {1, 2};
    const int32_t p2[] = {1};
    CHECK_EQ(path_pool_add(pool, p0, 2), 0, "first path index");
    CHECK_EQ(path_pool_add(pool, p1, 2), 1, "second path index");
    CHECK_EQ(path_pool_add(pool, p2, 1), 2, "third path index");
    CHECK_EQ(pool->num_paths, 3, "num_paths");
    CHECK_EQ(path_pool_length(pool, 0), 2, "path length");
    CHECK_EQ(path_pool_length(pool, 2), 1, "single-edge path length");
    CHECK_MSG(path_pool_edges(pool, 0)[0] == 0 && path_pool_edges(pool, 0)[1] == 1,
              "path edges are stored unchanged");
    CHECK_EQ(path_pool_edge_usage(pool, 0), 1, "edge usage");
    CHECK_EQ(path_pool_edge_usage(pool, 1), 3, "edge used by all three paths");
    CHECK_EQ(path_pool_edge_usage(pool, 3), 0, "unused edge");
    path_pool_free(pool);
}

/* ------------------------------------------------------------------------- */
/* graph.h / graph_load.c                                                    */
/* ------------------------------------------------------------------------- */

static void test_graph_parsing(void)
{
    Graph *g = write_and_load("parsing", PARSING_EDGES, parsing_od, false);

    CHECK_EQ(g->stats.raw_edges, 5, "all five edge lines counted");
    CHECK_EQ(g->stats.self_loops_removed, 0, "no self-loops in this file");
    CHECK_EQ(g->stats.parallel_edges_removed, 1, "e2/e3 duplicate merged");
    CHECK_EQ(g->stats.duplicate_od_entries, 1, "repeated O-D entry counted");
    CHECK_EQ(g->stats.isolated_od_nodes, 2, "OD8 and OD9 are isolated");
    CHECK_EQ(g->num_nodes, 6, "four network nodes + two isolated O-D nodes");
    CHECK_EQ(g->num_edges, 4, "one of the parallel pair dropped");
    CHECK_EQ(g->num_od_nodes, 4, "distinct O-D nodes");
    CHECK_EQ(graph_num_od_pairs(g), 6, "|W| = 4*3/2");

    /* Nodes keep file order of first appearance: OD1 OD2 X Y OD9 OD8. */
    CHECK_MSG(strcmp(g->node_label[0], "OD1") == 0, "BOM dropped from first label");
    CHECK_MSG(strcmp(g->node_label[3], "Y") == 0, "label after CRLF line intact");
    CHECK_MSG(strcmp(g->edge_label[0], "e1") == 0, "edge labels kept in file order");
    CHECK_MSG(strcmp(g->edge_label[3], "e5") == 0, "last line parsed without final newline");

    CHECK_EQ(g->od_rank[0], 0, "OD1 rank");
    CHECK_EQ(g->od_rank[2], -1, "X is not an O-D node");
    CHECK_EQ(g->od_rank[4], 2, "OD9 rank");
    CHECK_EQ(g->locked_edges[0], 0, "e1 is locked (both endpoints O-D)");
    CHECK_EQ(g->num_locked_edges, 1, "only e1 is locked");

    /* CSR adjacency: OD1 -> {e1, e2}, OD2 -> {e1, e5}. */
    CHECK_EQ(g->adjacency_start[1] - g->adjacency_start[0], 2, "degree of OD1");
    CHECK_EQ(g->adjacency_start[2] - g->adjacency_start[1], 2, "degree of OD2");
    CHECK_EQ(g->adjacency_edges[g->adjacency_start[1] + 1], 3, "OD2 is also incident to e5");

    /* Cyclic partner assignment must cover every pair exactly once. */
    CHECK_EQ(count_assigned_od_pairs(g), 6, "every O-D pair assigned exactly once");
    CHECK_EQ(sum_od_partner_counts(g), 6, "partner counts add up to |W|");
    CHECK_MSG(!graph_is_od_partner(g, 0, 0), "a root is not its own partner");

    graph_free(g);
}

static void test_graph_contraction(void)
{
    Graph *line = write_and_load("line", LINE_EDGES, line_od, true);
    CHECK_EQ(line->stats.nodes_contracted, 2, "x and y contracted");
    CHECK_EQ(line->num_nodes, 3, "A B C remain");
    CHECK_EQ(line->num_edges, 2, "two merged links remain");
    CHECK_MSG(strcmp(line->edge_label[0], "l1") == 0, "merged edge keeps the first label");
    CHECK_MSG(strcmp(line->edge_label[1], "l3") == 0, "second merged edge keeps its label");
    CHECK_MSG(graph_other_endpoint(line, 0, line->edge_node_a[0]) == line->od_nodes[1],
              "l1 now joins A and B");
    CHECK_EQ(line->num_locked_edges, 2, "both contracted links are locked");
    CHECK_EQ(brute_force_tclp_optimum(line), 2, "optimum of the path graph");
    graph_free(line);

    Graph *ring = write_and_load("ring", RING_EDGES, ring_od, true);
    CHECK_EQ(ring->stats.nodes_contracted, 3, "d1, d2 and c contracted");
    CHECK_EQ(ring->stats.ring_edges_removed, 2, "dead-end ring deleted");
    CHECK_EQ(ring->num_nodes, 2, "A and B remain");
    CHECK_EQ(ring->num_edges, 1, "a single link A-B remains");
    CHECK_EQ(ring->num_locked_edges, 1, "the link is locked");
    graph_free(ring);

    /* Without contraction the raw structure survives. */
    Graph *raw = write_and_load("line", LINE_EDGES, line_od, false);
    CHECK_EQ(raw->num_nodes, 5, "no node contracted");
    CHECK_EQ(raw->num_edges, 4, "no edge contracted");
    CHECK_EQ(raw->num_locked_edges, 0, "no edge has two O-D endpoints");
    CHECK_EQ(brute_force_tclp_optimum(raw), 2, "same optimum as contracted");
    graph_free(raw);
}

static void test_diamond_structure(void)
{
    Graph *g = write_and_load("diamond", DIAMOND_EDGES, diamond_od, false);
    CHECK_EQ(g->num_nodes, 7, "three O-D nodes and four intersections");
    CHECK_EQ(g->num_edges, 8, "eight links");
    CHECK_EQ(g->num_locked_edges, 0, "no O-D node pair is adjacent");
    CHECK_EQ(graph_num_od_pairs(g), 3, "three O-D pairs");
    CHECK_EQ(count_assigned_od_pairs(g), 3, "partner assignment is exact");
    CHECK_EQ(brute_force_tclp_optimum(g), 4, "reference optimum of the diamond");
    graph_free(g);
}

/* ------------------------------------------------------------------------- */
/* cover.h                                                                   */
/* ------------------------------------------------------------------------- */

static void test_cover_counters(void)
{
    Graph *g = write_and_load("diamond", DIAMOND_EDGES, diamond_od, false);
    PathPool *pool = path_pool_create(g->num_edges);
    add_diamond_pool_paths(pool);
    Cover *cover = cover_create(g, pool);

    cover_reset(cover);
    CHECK_EQ(cover->num_selected, 0, "empty cover");
    CHECK_EQ(cover->num_paths, 5, "counters know all pool paths");
    CHECK_EQ(cover->num_uncovered, 5, "all paths uncovered");
    CHECK_EQ(cover->uncovered_paths_through[0], 2, "two uncovered paths through e1");
    CHECK_EQ(cover->uncovered_paths_through[1], 1, "one uncovered path through e2");
    CHECK_EQ(cover->uncovered_paths_through[3], 1, "one uncovered path through e4");

    cover_add_edge(cover, 2); /* e3: covers p0 and p4 */
    CHECK_EQ(cover->num_uncovered, 3, "adding e3 covers two paths");
    CHECK_EQ(cover->path_hits[0], 1, "p0 hit once");
    CHECK_EQ(cover->path_hits[4], 1, "p4 hit once");
    CHECK_EQ(cover->uncovered_paths_through[0], 0, "e1 lies on no uncovered path now");
    CHECK_EQ(cover->uncovered_paths_through[1], 1, "e2 still on uncovered p1");

    cover_add_edge(cover, 1); /* e2 covers p1 */
    CHECK_EQ(cover->num_uncovered, 2, "e2 covers p1");
    CHECK_EQ(cover->uncovered_paths_through[3], 0, "e4 still unused");

    cover_remove_edge(cover, 2); /* exact inverse */
    CHECK_EQ(cover->num_uncovered, 4, "removing e3 uncovers p0 and p4");
    CHECK_EQ(cover->uncovered_paths_through[0], 2, "e1 counters restored");
    cover_remove_edge(cover, 1);
    CHECK_EQ(cover->num_uncovered, 5, "add/remove are inverse");

    /* Paths added after the reset are picked up by cover_sync(). */
    const int32_t late[] = {5};
    path_pool_add(pool, late, 1);
    cover_sync(cover);
    CHECK_EQ(cover->num_paths, 6, "cover_sync accounts for the new path");
    CHECK_EQ(cover->num_uncovered, 6, "the new path starts uncovered");
    CHECK_EQ(cover->uncovered_paths_through[5], 3, "p2, p4 and the late path run through e6");
    cover_add_edge(cover, 5);
    CHECK_EQ(cover->num_uncovered, 3, "e6 covers p2, p4 and the late path at once");

    cover_free(cover);
    path_pool_free(pool);
    graph_free(g);
}

static void test_cover_redundant_and_feasibility(void)
{
    Graph *g = write_and_load("diamond", DIAMOND_EDGES, diamond_od, false);
    PathPool *pool = path_pool_create(g->num_edges);
    add_diamond_pool_paths(pool);
    Cover *cover = cover_create(g, pool);
    cover_reset(cover);
    const int32_t cover_edges[] = {0, 1, 2, 4, 6}; /* e1 e2 e3 e5 e7 */
    select_edges(cover, cover_edges, 5);
    CHECK_EQ(cover->num_uncovered, 0, "the five-edge cover covers every pool path");

    int32_t redundant[16];
    CHECK_EQ(cover_collect_redundant(cover, redundant), 2, "R_C = {e1, e3}");
    CHECK_EQ(redundant[0], 0, "first redundant edge");
    CHECK_EQ(redundant[1], 2, "second redundant edge");

    /* Independent union-find feasibility check on the graph itself. */
    uint8_t empty[8] = {0};
    CHECK_EQ(cover_count_connected_od_pairs(g, empty), 3, "empty cover leaves all pairs");
    uint8_t block_cut[8] = {0, 0, 1, 1, 0, 1, 1, 0}; /* e3 e4 e6 e7 */
    CHECK_EQ(cover_count_connected_od_pairs(g, block_cut), 0, "a cut hits every route");
    CHECK_EQ(brute_force_tclp_optimum(g), 4, "so four edges separate all pairs");

    cover_free(cover);
    path_pool_free(pool);
    graph_free(g);
}

/* ------------------------------------------------------------------------- */
/* exact_bruteforce.c (the exact.h stand-in)                                 */
/* ------------------------------------------------------------------------- */

static void test_exact_stand_in(void)
{
    Graph *g = write_and_load("diamond", DIAMOND_EDGES, diamond_od, false);
    CHECK_MSG(exact_solver_description() != nullptr && exact_solver_description()[0] != '\0',
              "the stand-in describes itself");

    /* Four disjoint two-edge paths: one edge from each is needed, so the
     * unicost optimum is 4. */
    PathPool *disjoint = path_pool_create(g->num_edges);
    const int32_t q0[] = {0, 2};
    const int32_t q1[] = {1, 3};
    const int32_t q2[] = {4, 5};
    const int32_t q3[] = {6, 7};
    path_pool_add(disjoint, q0, 2);
    path_pool_add(disjoint, q1, 2);
    path_pool_add(disjoint, q2, 2);
    path_pool_add(disjoint, q3, 2);
    Cover *cover = cover_create(g, disjoint);
    cover_reset(cover);
    ExactResult res = exact_cover(g_solver, cover, 10.0, nullptr, 0.0);
    CHECK_MSG(res.has_solution, "solver finds a completion");
    CHECK_EQ(res.solver_calls, 1, "one model solved");
    CHECK_EQ((long long)res.objective_bound, 4, "optimum of four disjoint paths");
    CHECK_MSG(res.runtime >= 0.0, "runtime reported");
    CHECK_EQ(cover->num_selected, 4, "four edges added");
    CHECK_EQ(cover->num_uncovered, 0, "every path covered");

    /* Already covered -> solver not called at all. */
    ExactResult idle = exact_cover(g_solver, cover, 10.0, nullptr, 0.0);
    CHECK_EQ(idle.solver_calls, 0, "nothing to solve -> no solver call");
    CHECK_EQ((long long)idle.objective_bound, 0, "bound 0 when nothing is uncovered");
    CHECK_MSG(idle.has_solution, "trivially feasible");
    CHECK_EQ(cover->num_selected, 4, "cover untouched");

    /* known_bound = optimum: the same solution is still returned. */
    cover_reset(cover);
    ExactResult bounded = exact_cover(g_solver, cover, 10.0, nullptr, 4.0);
    CHECK_MSG(bounded.has_solution, "known bound still solved");
    CHECK_EQ((long long)bounded.objective_bound, 4, "solution size equals the known bound");
    CHECK_EQ(cover->num_uncovered, 0, "still a feasible completion");

    cover_free(cover);
    path_pool_free(disjoint);

    /* Overlapping paths {e2}, {e1,e2}, {e2,e3}: one edge (e2) is optimal. */
    PathPool *overlap = path_pool_create(g->num_edges);
    const int32_t r0[] = {1};
    const int32_t r1[] = {0, 1};
    const int32_t r2[] = {1, 2};
    path_pool_add(overlap, r0, 1);
    path_pool_add(overlap, r1, 2);
    path_pool_add(overlap, r2, 2);
    cover = cover_create(g, overlap);
    cover_reset(cover);
    res = exact_cover(g_solver, cover, 10.0, nullptr, 0.0);
    CHECK_EQ((long long)res.objective_bound, 1, "single shared edge is optimal");
    CHECK_EQ(cover->num_selected, 1, "one edge added");
    CHECK_MSG(cover->selected[1], "the shared edge e2 was picked");

    /* Warm start bigger than the optimum must not be preferred over it. */
    cover_reset(cover);
    uint8_t warm[8] = {0};
    warm[0] = 1;
    warm[1] = 1;
    warm[2] = 1;
    ExactResult warm_result = exact_cover(g_solver, cover, 10.0, warm, 0.0);
    CHECK_MSG(warm_result.has_solution, "warm start case solved");
    CHECK_EQ((long long)warm_result.objective_bound, 1, "optimum beats the warm start");
    cover_free(cover);
    path_pool_free(overlap);
    graph_free(g);
}

/* ------------------------------------------------------------------------- */
/* cleaning.h                                                                */
/* ------------------------------------------------------------------------- */

static void rebuild_redundant_cover(Cover *cover)
{
    cover_reset(cover);
    const int32_t cover_edges[] = {0, 1, 2, 4, 6};
    select_edges(cover, cover_edges, 5);
}

static void test_cleaning(void)
{
    Graph *g = write_and_load("diamond", DIAMOND_EDGES, diamond_od, false);
    PathPool *pool = path_pool_create(g->num_edges);
    add_diamond_pool_paths(pool);
    Cover *cover = cover_create(g, pool);
    Rng rng;
    rng_seed(&rng, 4242);

    rebuild_redundant_cover(cover);
    CleaningStats stats = clean_cover(cover, CLEANING_NONE, g_solver, 5.0, &rng);
    CHECK_EQ(cover->num_selected, 5, "CLEANING_NONE keeps the cover");
    CHECK_EQ(stats.solver_calls, 0, "and uses no solver");

    rebuild_redundant_cover(cover);
    stats = clean_cover(cover, CLEANING_RANDOM, nullptr, 0.0, &rng);
    CHECK_EQ(cover->num_selected, 4, "RC removes the redundant edge");
    CHECK_EQ(cover->num_uncovered, 0, "RC keeps every path covered");
    CHECK_EQ(stats.solver_calls, 0, "RC uses no solver");

    rebuild_redundant_cover(cover);
    stats = clean_cover(cover, CLEANING_GREEDY, nullptr, 0.0, &rng);
    CHECK_EQ(cover->num_selected, 4, "GC re-covers with one edge");
    CHECK_EQ(cover->num_uncovered, 0, "GC keeps coverage");

    rebuild_redundant_cover(cover);
    stats = clean_cover(cover, CLEANING_EXACT, g_solver, 5.0, &rng);
    CHECK_EQ(cover->num_selected, 4, "EC re-covers optimally");
    CHECK_EQ(cover->num_uncovered, 0, "EC keeps coverage");
    CHECK_EQ(stats.solver_calls, 1, "EC calls the solver once");
    CHECK_MSG(stats.solver_time >= 0.0, "EC reports solver time");

    cover_free(cover);
    path_pool_free(pool);

    /* Locked edges are never cleaned away. */
    Graph *line = write_and_load("line", LINE_EDGES, line_od, true);
    PathPool *line_pool = path_pool_create(line->num_edges);
    const int32_t path01[] = {0, 1};
    path_pool_add(line_pool, path01, 2);
    Cover *locked = cover_create(line, line_pool);
    cover_reset(locked);
    CHECK_EQ(locked->num_selected, 2, "reset selects the locked edges");
    CHECK_EQ(locked->num_uncovered, 0, "the locked edges cover the path");
    clean_cover(locked, CLEANING_RANDOM, nullptr, 0.0, &rng);
    CHECK_EQ(locked->num_selected, 2, "locked edges stay");
    clean_cover(locked, CLEANING_EXACT, g_solver, 5.0, &rng);
    CHECK_EQ(locked->num_selected, 2, "locked edges stay after EC");
    cover_free(locked);
    path_pool_free(line_pool);
    graph_free(line);
    graph_free(g);
}

/* ------------------------------------------------------------------------- */
/* separation.h                                                              */
/* ------------------------------------------------------------------------- */

static void test_separation(void)
{
    Graph *g = write_and_load("diamond", DIAMOND_EDGES, diamond_od, false);
    PathPool *pool = path_pool_create(g->num_edges);
    Cover *cover = cover_create(g, pool);
    Separator *sep = separator_create(g);
    Rng rng;
    rng_seed(&rng, 20261119);

    /* An empty cover leaves uncovered paths for the roots O1 and O2. */
    cover_reset(cover);
    int32_t found = separator_run(sep, cover, pool, SPREAD_MLC, 0, &rng);
    CHECK_MSG(found >= 1 && found <= 4, "first Cut phase finds O-D paths");
    cover_sync(cover);
    CHECK_EQ(cover->num_uncovered, found, "new paths enter uncovered");
    for (int32_t p = 0; p < pool->num_paths; p++) {
        CHECK_MSG(path_pool_length(pool, p) >= 1, "stored paths are non-empty");
    }

    /* Alternate greedy covering and Cut phases until the cover is feasible. */
    int32_t rounds = 0;
    while (found > 0 && rounds < 50) {
        greedy_cover(cover, 0, &rng);
        CHECK_EQ(cover->num_uncovered, 0, "greedy covers the pool after each Cut phase");
        found = separator_run(sep, cover, pool, SPREAD_MLC, 0, &rng);
        cover_sync(cover);
        rounds++;
    }
    CHECK_MSG(found == 0, "Cut phase ends when the cover is feasible");
    CHECK_EQ(cover_count_connected_od_pairs(g, cover->selected), 0,
             "the final cover separates all O-D pairs");

    /* One path per root at most. */
    PathPool *small_pool = path_pool_create(g->num_edges);
    Cover *small_cover = cover_create(g, small_pool);
    cover_reset(small_cover);
    found = separator_run(sep, small_cover, small_pool, SPREAD_MLC, 1, &rng);
    CHECK_MSG(found >= 1 && found <= 2, "max_paths_per_root is respected");
    cover_free(small_cover);
    path_pool_free(small_pool);

    /* The random spread criterion also runs and finds paths. */
    PathPool *random_pool = path_pool_create(g->num_edges);
    Cover *random_cover = cover_create(g, random_pool);
    cover_reset(random_cover);
    found = separator_run(sep, random_cover, random_pool, SPREAD_RANDOM, 0, &rng);
    CHECK_MSG(found >= 1 && found <= 4, "SPREAD_RANDOM finds uncovered paths too");
    cover_free(random_cover);
    path_pool_free(random_pool);

    separator_free(sep);
    cover_free(cover);
    path_pool_free(pool);
    graph_free(g);
}

/* ------------------------------------------------------------------------- */
/* algorithms.h                                                              */
/* ------------------------------------------------------------------------- */

static void test_algorithm_parameters(void)
{
    SolverParams params;
    solver_params_defaults(&params, ALGORITHM_PCA);
    CHECK_EQ(params.algorithm, ALGORITHM_PCA, "algorithm stored");
    CHECK_MSG(fabs(params.alpha0 - 1.0 / 16.0) < 1e-12, "PCA alpha0 default");
    CHECK_EQ(params.initial_paths_per_root, 0, "PCA uses all initial paths");
    CHECK_EQ(params.paths_per_root, 5, "PCA later paths per root");
    CHECK_EQ(params.spread, SPREAD_MLC, "MLC default");
    CHECK_EQ(params.cleaning, CLEANING_EXACT, "exact cleaning default");
    solver_params_defaults(&params, ALGORITHM_EXACT);
    CHECK_EQ(params.initial_paths_per_root, 20, "EXACT initial paths per root");
    CHECK_EQ(params.paths_per_root, 10, "EXACT paths per root");
    CHECK_MSG(strcmp(algorithm_name(ALGORITHM_PCA), "pca") == 0, "name of the PCA");
    CHECK_MSG(strcmp(algorithm_name(ALGORITHM_EXACT), "exact") == 0, "name of the exact method");
    CHECK_MSG(strcmp(algorithm_name(ALGORITHM_GREEDY), "greedy") == 0, "name of the greedy");
}

static void test_algorithms_on_diamond(void)
{
    Graph *g = write_and_load("diamond", DIAMOND_EDGES, diamond_od, false);
    const int32_t optimum = brute_force_tclp_optimum(g);
    CHECK_EQ(optimum, 4, "reference optimum");
    RunResult result;
    SolverParams params;

    /* The exact framework with the brute-force solver must prove optimality. */
    solver_params_defaults(&params, ALGORITHM_EXACT);
    params.time_limit_total = 60.0;
    run_algorithm(g, &params, g_solver, 42, &result);
    CHECK_EQ(result.connected_pairs, 0, "exact run produces a feasible cover");
    CHECK_EQ(result.cost, optimum, "exact run finds the optimum");
    CHECK_MSG(!isnan(result.lower_bound), "exact run proves a bound");
    CHECK_MSG(result.lower_bound >= optimum - 1e-9 && result.lower_bound <= result.cost + 1e-9,
              "bound sits between the optimum and the cost");
    CHECK_MSG(result.exact_solves >= 1, "the solver was used");
    CHECK_MSG(result.iterations >= 1 && result.num_paths >= 1, "framework iterations ran");
    CHECK_MSG(result.time_total >= 0.0 && result.time_solver >= 0.0,
              "timing statistics are non-negative");
    run_result_free(&result);
    CHECK_MSG(result.selected == nullptr, "run_result_free releases the selection");

    /* The PCA must also be feasible and never better than the optimum. */
    solver_params_defaults(&params, ALGORITHM_PCA);
    params.time_limit_total = 60.0;
    run_algorithm(g, &params, g_solver, 7, &result);
    CHECK_EQ(result.connected_pairs, 0, "PCA run is feasible");
    CHECK_MSG(result.cost >= optimum, "PCA cannot beat the optimum");
    CHECK_MSG(isnan(result.lower_bound) || (double)result.cost >= result.lower_bound,
              "PCA cost respects its bound");
    int32_t pca_cost = result.cost;
    run_result_free(&result);

    /* Same seed, same result: runs are reproducible. */
    run_algorithm(g, &params, g_solver, 7, &result);
    CHECK_EQ(result.cost, pca_cost, "same seed reproduces the same cost");
    run_result_free(&result);

    /* Pure greedy needs no solver at all. */
    solver_params_defaults(&params, ALGORITHM_GREEDY);
    params.cleaning = CLEANING_NONE;
    run_algorithm(g, &params, nullptr, 3, &result);
    CHECK_EQ(result.connected_pairs, 0, "greedy run is feasible");
    CHECK_MSG(result.cost >= optimum, "greedy cost above optimum");
    CHECK_MSG(isnan(result.lower_bound), "greedy run computes no bound");
    run_result_free(&result);
    graph_free(g);
}

static void test_algorithms_with_locked_edges(void)
{
    Graph *g = write_and_load("line", LINE_EDGES, line_od, true);
    const int32_t optimum = brute_force_tclp_optimum(g);
    CHECK_EQ(optimum, 2, "path graph needs its two (locked) links");
    SolverParams params;
    solver_params_defaults(&params, ALGORITHM_EXACT);
    params.time_limit_total = 60.0;
    RunResult result;
    run_algorithm(g, &params, g_solver, 1, &result);
    CHECK_EQ(result.connected_pairs, 0, "exact run on the locked graph is feasible");
    CHECK_EQ(result.cost, optimum, "exact run respects locked edges");
    CHECK_MSG(result.lower_bound >= optimum - 1e-9, "bound respects locked edges");
    run_result_free(&result);
    graph_free(g);
}

/* ------------------------------------------------------------------------- */
/* Real repository instance (skipped when not found)                          */
/* ------------------------------------------------------------------------- */

static void test_real_instance(void)
{
    const char *edges = "instances/AC/edges.txt";
    const char *ods = "instances/AC/od_nodes.txt";
    FILE *probe = fopen(edges, "rb");
    if (probe == nullptr) {
        fprintf(stderr, "note: instances/AC not found, skipping the real-instance check\n");
        return;
    }
    fclose(probe);
    Graph *g = graph_load(edges, ods, true);
    CHECK_EQ(g->stats.raw_edges, 114, "AC has 114 raw edge lines");
    CHECK_EQ(g->num_od_nodes, 20, "AC has 20 O-D nodes");
    CHECK_EQ(graph_num_od_pairs(g), 190, "AC has 20*19/2 O-D pairs");
    CHECK_MSG(g->num_edges > 0 && g->num_nodes > 0, "AC loaded successfully");
    CHECK_EQ(count_assigned_od_pairs(g), 190, "AC partner assignment is exact");
    graph_free(g);
}

/* ------------------------------------------------------------------------- */

int main(void);

int main(void)
{
    const char *base = getenv("TMPDIR");
    if (base == nullptr || base[0] == '\0') {
        base = "/tmp";
    }
    char tmpl[4096];
    snprintf(tmpl, sizeof tmpl, "%s/location_pca_tests_XXXXXX", base);
    if (mkdtemp(tmpl) == nullptr) {
        fprintf(stderr, "unit_tests: cannot create a temporary directory under '%s'\n", base);
        return EXIT_FAILURE;
    }
    snprintf(g_tmpdir, sizeof g_tmpdir, "%s", tmpl);
    CHECK_MSG(log_get_level() == LOG_INFO, "log level is LOG_INFO by default");
    log_set_level(LOG_QUIET);
    g_solver = exact_solver_create(1, false);

    test_util();
    test_rng();
    test_path_pool();
    test_graph_parsing();
    test_graph_contraction();
    test_diamond_structure();
    test_cover_counters();
    test_cover_redundant_and_feasibility();
    test_exact_stand_in();
    test_cleaning();
    test_separation();
    test_algorithm_parameters();
    test_algorithms_on_diamond();
    test_algorithms_with_locked_edges();
    test_real_instance();

    exact_solver_free(g_solver);
    cleanup_tmpdir();

    if (g_failures > 0) {
        fprintf(stderr, "%d of %d checks FAILED\n", g_failures, g_checks);
        return EXIT_FAILURE;
    }
    printf("all %d checks passed\n", g_checks);
    return EXIT_SUCCESS;
}
