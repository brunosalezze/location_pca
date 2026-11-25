/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "cut.h"
#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>

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
                fprintf(stderr, "planar-cut: %d terminals exceed the limit of %d for Klein–Marx%s\n",
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

int main(int argc, char **argv)
{
    const char *path = NULL, *dir = NULL, *edges = NULL, *ods = NULL;
    bool only_planarity = false, paper_only = false;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0) {
            puts("usage: planar-cut [--klein-marx-only] [--crossings N] [FILE|-]\n"
                 "       planar-cut --tclp DIR | --edges EDGES --od OD_NODES\n"
                 "Read a weighted planar graph with a rotation system; see README.md.\n"
                 "--tclp DIR     solve a traffic counting location instance (DIR holds\n"
                 "               edges.txt and od_nodes.txt); the planarity of the\n"
                 "               network is tested before the cut is computed\n"
                 "--planarity-only  only test and report planarity in TCLP mode\n"
                 "--klein-marx-only  disable the exact primal DP fast path\n"
                 "--crossings N  cap on interior crossings per terminal (default 32;\n"
                 "               the paper proves 24 suffices for optimality)");
            return 0;
        }
        if (strcmp(argv[i], "--crossings") == 0) {
            paper_only = true;
            if (++i == argc) die("--crossings needs a value");
            char *end;
            long v = strtol(argv[i], &end, 10);
            if (*end || v < 1 || v > 1000) die("--crossings must be in 1..1000");
            g_crossings = (int)v;
            continue;
        }
        if (strcmp(argv[i], "--klein-marx-only") == 0) {
            paper_only = true;
            continue;
        }
        if (strcmp(argv[i], "--planarity-only") == 0) {
            only_planarity = true;
            continue;
        }
        if (strcmp(argv[i], "--tclp") == 0) {
            if (++i == argc) die("--tclp needs a directory");
            dir = argv[i];
            continue;
        }
        if (strcmp(argv[i], "--edges") == 0) {
            if (++i == argc) die("--edges needs a file");
            edges = argv[i];
            continue;
        }
        if (strcmp(argv[i], "--od") == 0) {
            if (++i == argc) die("--od needs a file");
            ods = argv[i];
            continue;
        }
        if (path) die("usage: planar-cut [--crossings N] [FILE|-]");
        path = argv[i];
    }
    if ((dir && (edges || ods)) || (!dir && (!!edges != !!ods)))
        die("use either --tclp DIR or --edges EDGES --od OD_NODES");
    if (only_planarity && !dir && !edges)
        die("--planarity-only requires TCLP input");
    if (dir || edges) {
        if (path) die("no input file with --tclp or --edges");
        Inst ins = tclp_read(dir, edges, ods);
        fprintf(stderr, "instance: %d nodes, %d edges, %d terminals\n",
                ins.g.n, ins.g.m, ins.k);
        /* Planarity first: the algorithm needs a planar embedding. */
        Planarity verdict = embed(&ins.g);
        if (verdict != PLANAR) {
            puts(verdict == NONPLANAR ? "planar no" : "planar unknown");
            if (verdict == PLANARITY_UNKNOWN)
                fputs("planar-cut: planarity search budget exhausted\n", stderr);
            tclp_free(&ins);
            return EXIT_FAILURE;
        }
        puts("planar yes");
        fflush(stdout);
        if (only_planarity) { tclp_free(&ins); return 0; }
        Answer result = components(&ins.g, ins.term, ins.k, paper_only);
        verify(&ins.g, ins.term, ins.k, &result);
        int count = 0;
        for (int e = 0; e < ins.g.m; ++e) count += result.cut[e];
        printf("cost %" PRIu64 "\ncut %d\n", result.cost, count);
        for (int e = 0; e < ins.g.m; ++e) if (result.cut[e]) printf("%s\n", ins.edge_id[e]);
        free(result.cut);
        tclp_free(&ins);
        return 0;
    }
    FILE *in = stdin;
    if (path && strcmp(path, "-") != 0) {
        in = fopen(path, "r");
        if (!in) { perror(path); return EXIT_FAILURE; }
    }
    int *term, k;
    Plane g = load(in, &term, &k);
    if (in != stdin) fclose(in);
    Answer result = components(&g, term, k, paper_only);
    verify(&g, term, k, &result);
    int count = 0;
    for (int e = 0; e < g.m; ++e) count += result.cut[e];
    printf("cost %" PRIu64 "\ncut %d\n", result.cost, count);
    for (int e = 0; e < g.m; ++e) if (result.cut[e]) printf("%d\n", e+1);
    free(result.cut); free(term); plane_free(&g);
    return 0;
}
