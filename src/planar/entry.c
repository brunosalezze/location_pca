/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
/**
 * @file entry.c
 * @brief The planar algorithm: planarity check, primal dynamic program and
 *        Klein-Marx fallback, driven by the parsed command line.
 *
 * planar_run_all() replaces the old standalone planar-cut main(). The
 * traffic-counting mode reads edges.txt and od_nodes.txt, tests planarity
 * first (the verdict goes to stderr), reports one TSV row per run on stdout
 * and writes the cut through --solution-dir. Without od_nodes.txt the input
 * is a weighted graph with a rotation system (stdin with "-", see
 * tests/planar/README.md), solved in the human-readable test format.
 */
#include "entry.h"

#include "cut.h"
#include "../report.h"
#include "../timer.h"
#include "../util.h"

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <sys/stat.h>

int g_crossings = 32;

static bool word(FILE *in, char *buf, size_t cap)
{
    int ch;
    do {
        ch = fgetc(in);
        if (ch == '#') {
            do ch = fgetc(in); while (ch != '\n' && ch != EOF);
        }
    } while (ch != EOF && isspace(ch));
    if (ch == EOF) {
        if (ferror(in)) die("input read failed");
        return false;
    }
    size_t n = 0;
    do {
        if (n+1 == cap) die("input token too long");
        buf[n++] = (char)ch;
        ch = fgetc(in);
    } while (ch != EOF && !isspace(ch));
    buf[n] = '\0';
    return true;
}

static uint64_t number(FILE *in, bool *neg)
{
    char buf[64];
    if (!word(in, buf, sizeof buf)) die("unexpected end of input");
    char *digits = buf;
    *neg = *digits == '-';
    if (*digits == '-' || *digits == '+') ++digits;
    if (!isdigit((unsigned char)*digits)) die("expected an integer");
    errno = 0;
    char *end;
    uintmax_t value = strtoumax(digits, &end, 10);
    if (errno || *end || value > UINT64_MAX) die("invalid or out-of-range integer");
    return (uint64_t)value;
}

static int index_in(FILE *in, int lo, int hi)
{
    bool neg;
    uint64_t v = number(in, &neg);
    if (neg || v < (uint64_t)lo || v > (uint64_t)hi) die("integer outside allowed range");
    return (int)v;
}

static Plane load(FILE *in, int **term, int *k)
{
    int n = index_in(in, 1, INT_MAX), m = index_in(in, 0, INT_MAX/2);
    *k = index_in(in, 0, n);
    *term = mem((size_t)*k, sizeof(int));
    bool *seen = mem((size_t)n, sizeof(bool));
    for (int i = 0; i < *k; ++i) {
        int v = index_in(in, 0, n-1);
        if (seen[v]) die("duplicate terminal");
        seen[v] = true;
        (*term)[i] = v;
    }
    free(seen);
    Plane g = plane(n, m);
    Cost total = 0;
    for (int e = 0; e < m; ++e) {
        g.at[2*e] = index_in(in, 0, n-1);
        g.at[2*e+1] = index_in(in, 0, n-1);
        bool neg;
        g.w[e] = number(in, &neg);
        if (neg || !g.w[e]) die("edge costs must be positive");
        total = plus(total, g.w[e]);
        if (total == INF) die("sum of edge costs must be below UINT64_MAX");
    }
    for (int v = 0; v < n; ++v) {
        int degree = index_in(in, 0, 2*m), prev = -1;
        for (int j = 0; j < degree; ++j) {
            bool neg;
            uint64_t id = number(in, &neg);
            if (!id || id > (uint64_t)m) die("invalid rotation edge ID");
            int d = 2*((int)id-1)+(neg ? 1 : 0);
            if (g.at[d] != v || g.rot[d] >= 0 || d == g.head[v])
                die("wrong vertex or repeated dart in rotation");
            if (prev >= 0) g.rot[prev] = d;
            else g.head[v] = d;
            prev = d;
        }
        if (prev >= 0) g.rot[prev] = g.head[v];
    }
    for (int d = 0; d < 2*m; ++d) if (g.rot[d] < 0) die("dart missing from rotation");
    char buf[64];
    if (word(in, buf, sizeof buf)) die("unexpected trailing input");
    return g;
}

/* Group vertices and edges in linear time; the paper assumes each
   instance is connected. No cut edges are needed between components. */
static Answer components(const Plane *g, const int *term, int k, bool paper_only)
{
    size_t n = (size_t)g->n, m = (size_t)g->m;
    int *comp = mem(n, sizeof(int)), *local = mem(n, sizeof(int));
    int *q = mem(n, sizeof(int)), *vstart = mem(n+1, sizeof(int));
    int *edges = mem(m, sizeof(int)), *estart = mem(n+1, sizeof(int));
    int *local_e = mem(m, sizeof(int));
    int count = 0, front = 0, back = 0;
    for (int u = 0; u < g->n; ++u) comp[u] = -1;
    for (int u = 0; u < g->n; ++u) if (comp[u] < 0) {
        vstart[count] = back; comp[u] = count; q[back++] = u;
        while (front < back) {
            int v = q[front];
            local[v] = front++ - vstart[count];
            int s = g->head[v];
            if (s < 0) continue;
            int a = s;
            do {
                int w = g->at[a ^ 1];
                if (comp[w] < 0) { comp[w] = count; q[back++] = w; }
                a = g->rot[a];
            } while (a != s);
        }
        ++count;
    }
    vstart[count] = back;
    for (int e = 0; e < g->m; ++e) ++estart[comp[g->at[2*e]]+1];
    for (int i = 0; i < count; ++i) estart[i+1] += estart[i];
    int *pos = mem((size_t)count, sizeof(int));
    memcpy(pos, estart, (size_t)count*sizeof(int));
    for (int e = 0; e < g->m; ++e) {
        int c = comp[g->at[2*e]], p = pos[c]++;
        local_e[e] = p-estart[c]; edges[p] = e;
    }
    int *reps = mem((size_t)k, sizeof(int));
    Answer ans = {g->m, 0, mem(m, sizeof(bool))};
    for (int c = 0; c < count; ++c) {
        int ne = estart[c+1]-estart[c], nv = vstart[c+1]-vstart[c];
        if (!ne) continue;
        Plane h = plane(nv, ne);
        for (int j = 0; j < ne; ++j) {
            int e = edges[estart[c]+j];
            h.w[j] = g->w[e];
            for (int b = 0; b < 2; ++b) {
                int a = 2*e+b, r = g->rot[a], v = local[g->at[a]];
                h.at[2*j+b] = v; h.head[v] = 2*j+b;
                h.rot[2*j+b] = 2*local_e[r/2]+r%2;
            }
        }
        int nf;
        free(faces(&h, &nf));
        if ((int64_t)nv-ne+nf != 2) die("rotation system is not planar");
        int nr = 0;
        for (int i = 0; i < k; ++i) if (comp[term[i]] == c) reps[nr++] = local[term[i]];
        if (nr < 2) { plane_free(&h); continue; }
        Answer result;
        if (paper_only || nr == 2 || !solve_narrow(&h, reps, nr, &result)) {
            if (nr > MAX_TERMINALS) {
                fprintf(stderr, "location_pca: %d terminals exceed the limit of %d for Klein–Marx%s\n",
                        nr, MAX_TERMINALS, paper_only ? "" : " (primal DP exceeded its budget)");
                exit(EXIT_FAILURE);
            }
            result = solve(&h, reps, nr);
        }
        ans.cost += result.cost;
        for (int e = 0; e < ne; ++e) ans.cut[edges[estart[c]+e]] = result.cut[e];
        free(result.cut); plane_free(&h);
    }
    free(comp); free(local); free(q); free(vstart); free(estart);
    free(pos); free(edges); free(local_e); free(reps);
    return ans;
}

/* Independently verify the returned edge set on the original graph. */
static int root(int *p, int i)
{
    while (p[i] != i) { p[i] = p[p[i]]; i = p[i]; }
    return i;
}

static void verify(const Plane *g, const int *term, int k, const Answer *a)
{
    int *p = mem((size_t)g->n, sizeof(int));
    for (int v = 0; v < g->n; ++v) p[v] = v;
    Cost cost = 0;
    for (int e = 0; e < g->m; ++e) {
        if (a->cut[e]) cost += g->w[e];
        else p[root(p, g->at[2*e])] = root(p, g->at[2*e+1]);
    }
    if (cost != a->cost) die("internal error: cut cost mismatch");
    bool *seen = mem((size_t)g->n, sizeof(bool));
    for (int i = 0; i < k; ++i) {
        int r = root(p, term[i]);
        if (seen[r]) die("internal error: returned cut does not separate the terminals");
        seen[r] = true;
    }
    free(seen);
    free(p);
}

/* Writes the chosen stations in the solution file format of report.c. */
static void write_solution(const Options *options, uint64_t seed, const Inst *ins,
                           const Answer *a)
{
    if (mkdir(options->solution_dir, 0777) != 0 && errno != EEXIST)
        die("cannot create the solution directory");
    size_t size = strlen(options->solution_dir) + strlen(options->instance_name) + 48;
    char *path = mem(size, 1);
    snprintf(path, size, "%s/%s_planar_seed%" PRIu64 ".txt", options->solution_dir,
             options->instance_name, seed);
    FILE *file = fopen(path, "w");
    if (!file) die("cannot write the solution file");
    int count = 0;
    for (int e = 0; e < ins->g.m; ++e) count += a->cut[e];
    fprintf(file, "# location_pca %s solution\n", LOCATION_PCA_VERSION);
    fprintf(file, "# instance: %s\n# algorithm: planar\n# seed: %" PRIu64 "\n",
            options->instance_name, seed);
    fprintf(file, "# counting stations: %d\n# feasible: yes\n# lower bound: %" PRIu64 "\n",
            count, a->cost);
    fputs("# One edge identifier (first column of the edge file) per line.\n", file);
    for (int e = 0; e < ins->g.m; ++e)
        if (a->cut[e]) fprintf(file, "%s\n", ins->edge_id[e]);
    if (fclose(file) != 0) die("error while writing the solution file");
    free(path);
}

int planar_run_all(const Options *options)
{
    bool weighted = options->od_path == NULL;
    if (!weighted && options->print_header && !options->planarity_only) {
        report_table_header(stdout);
    }
    for (int32_t run = 0; run < options->runs; ++run) {
        uint64_t seed = options->seed + (uint64_t)run;
        if (!weighted) {
            Inst ins = tclp_read(NULL, options->edges_path, options->od_path);
            log_message(LOG_INFO, "instance %s: %d nodes, %d edges, %d O-D nodes",
                        options->instance_name, ins.g.n, ins.g.m, ins.k);
            /* The verdict belongs to the logs: stdout keeps the TSV contract. */
            Planarity verdict = embed(&ins.g);
            if (verdict != PLANAR) {
                fputs(verdict == NONPLANAR ? "planar no\n" : "planar unknown\n", stderr);
                if (verdict == PLANARITY_UNKNOWN)
                    fputs("location_pca: planarity search budget exhausted\n", stderr);
                tclp_free(&ins);
                return EXIT_FAILURE;
            }
            fputs("planar yes\n", stderr);
            if (options->planarity_only) { tclp_free(&ins); continue; }
            double begin = timer_now();
            Answer result = components(&ins.g, ins.term, ins.k, options->planar_paper_only);
            verify(&ins.g, ins.term, ins.k, &result);
            double elapsed = timer_now() - begin;
            if (result.cost > INT32_MAX) die("cut cost exceeds the table format");
            /* A minimum cut in a planar graph: the bound is the cost itself. */
            RunResult rr = {
                .cost = (int32_t)result.cost,
                .connected_pairs = 0,
                .lower_bound = (double)result.cost,
                .num_paths = 0,
                .iterations = 1,
                .exact_solves = 0,
                .time_total = elapsed,
            };
            report_table_row(stdout, options, run, seed, &rr);
            log_message(LOG_INFO, "run %d (seed %" PRIu64 "): %" PRIu64 " counting stations, %.2fs",
                        run, seed, result.cost, elapsed);
            if (options->solution_dir) write_solution(options, seed, &ins, &result);
            free(result.cut);
            tclp_free(&ins);
        } else {
            FILE *in = stdin;
            if (strcmp(options->edges_path, "-") != 0) {
                in = fopen(options->edges_path, "r");
                if (!in) die("cannot open the input file");
            }
            int *term, k;
            Plane g = load(in, &term, &k);
            if (in != stdin) fclose(in);
            Answer result = components(&g, term, k, options->planar_paper_only);
            verify(&g, term, k, &result);
            int count = 0;
            for (int e = 0; e < g.m; ++e) count += result.cut[e];
            printf("cost %" PRIu64 "\ncut %d\n", result.cost, count);
            for (int e = 0; e < g.m; ++e)
                if (result.cut[e]) printf("%d\n", e + 1);
            fflush(stdout);
            free(result.cut); free(term); plane_free(&g);
        }
    }
    return EXIT_SUCCESS;
}
