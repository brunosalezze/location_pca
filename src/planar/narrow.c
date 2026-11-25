/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "cut.h"
#include <errno.h>
#include <inttypes.h>
#include <limits.h>

/* Exact variable elimination on the primal graph with the terminals
   fixed. This is a fast path, independent of Klein–Marx enumeration.
   Assign vertices one of k labels and require different labels at the
   terminals. Cutting all differently labeled endpoints is feasible.
   Conversely, a feasible cut lets us label each residual component by
   its terminal (or arbitrarily if none), paying no more than the cut.
   Free-to-free edges give pair costs, and fixed endpoints give unary
   costs or constants. Minimize and reconstruct the labels in reverse.
   Bound the dense graph, table space, and estimated work before solving. */
enum { FREE_LIMIT = 1536 };
typedef struct { size_t rows; uint64_t work; bool stats; } Budget;

static uint64_t setting(const char *name, uint64_t fallback, uint64_t limit)
{
    const char *s = getenv(name);
    if (!s) return fallback;
    char *end;
    errno = 0;
    uintmax_t value = strtoumax(s, &end, 10);
    if (*s < '0' || *s > '9' || *end || errno || !value || value > limit) {
        fprintf(stderr, "location_pca: %s must be in 1..%" PRIu64 "\n", name, limit);
        exit(EXIT_FAILURE);
    }
    return (uint64_t)value;
}

static Budget budget(void)
{
    /* Charge every row for both a value and a traceback entry even when
       the value can be freed early. Leave headroom for the graph, dense
       pair matrix and the rest of the process under the 64 GB ceiling. */
    uint64_t mib = setting("CUT_DP_MIB", 1145, 60000);
    uint64_t bytes = mib*UINT64_C(1048576);
    if (bytes > SIZE_MAX) bytes = SIZE_MAX;
    return (Budget){(size_t)(bytes/(sizeof(Cost)+sizeof(int))),
                    setting("CUT_DP_WORK", UINT64_C(100000000000), UINT64_MAX),
                    getenv("CUT_DP_STATS") != NULL};
}

typedef struct {
    int v, degree, *scope;
    int bucket, next;
    size_t rows, factor_count;
    Cost *value;
    int *choice;
} Step;

static size_t cell(int n, int u, int v)
{
    return (size_t)u*(size_t)n+(size_t)v;
}

static bool plan(int n, int k, bool *adj, Step *step, int *rank, Budget limits)
{
    /* One unary table row per label and free vertex. */
    size_t tables = (size_t)n*(size_t)k;
    uint64_t work = 0;
    int width = 0;
    bool fits = true;
    for (int u = 0; u < n; ++u) rank[u] = -1;
    for (int i = 0; i < n; ++i) {
        int v = -1, degree = n, best_fill = INT_MAX;
        int neighbors[FREE_LIMIT];
        for (int u = 0; u < n; ++u) if (rank[u] < 0) {
            int count = 0, fill = 0;
            for (int p = 0; p < n; ++p)
                if (rank[p] < 0 && adj[cell(n, u, p)]) neighbors[count++] = p;
            for (int a = 0; a < count && fill <= best_fill; ++a)
                for (int b = 0; b < a; ++b) fill += !adj[cell(n, neighbors[a], neighbors[b])];
            if (fill < best_fill || (fill == best_fill && count < degree)) {
                v = u; degree = count; best_fill = fill;
                if (!fill && degree <= 1) break;
            }
        }
        if (degree > width) width = degree;
        size_t rows = 1;
        for (int j = 0; fits && j < degree; ++j) {
            if (rows > limits.rows/(size_t)k) {
                if (limits.stats) fprintf(stderr, "primal DP: %d^%d rows at step %d exceed the %zu-row budget\n",
                                          k, degree, i+1, limits.rows);
                fits = false;
                break;
            }
            rows *= (size_t)k;
        }
        if (fits && rows > limits.rows-tables) {
            if (limits.stats) fprintf(stderr, "primal DP: cumulative tables exceed %zu rows at step %d\n",
                                      limits.rows, i+1);
            fits = false;
        }
        if (fits) tables += rows;
        else if (!limits.stats) return false;
        step[i] = (Step){.v = v, .degree = degree, .rows = rows,
                         .scope = mem((size_t)degree, sizeof(int)), .bucket = -1};
        rank[v] = i;
        int count = 0;
        for (int p = 0; p < n; ++p)
            if (rank[p] < 0 && adj[cell(n, v, p)]) step[i].scope[count++] = p;
        for (int a = 0; a < degree; ++a)
            for (int b = 0; b < a; ++b) {
                int p = step[i].scope[a], q = step[i].scope[b];
                adj[cell(n, p, q)] = adj[cell(n, q, p)] = true;
            }
    }
    if (!fits) {
        fprintf(stderr, "primal DP: rejected elimination order has width %d (largest table: %d^%d rows)\n",
                width, k, width);
        return false;
    }
    for (int i = 0; i < n; ++i) {
        Step *s = &step[i];
        size_t unit = (size_t)s->degree+1+s->factor_count;
        if (s->rows > (limits.work-work)/unit/(size_t)k) {
            if (limits.stats) fprintf(stderr, "primal DP: work estimate exceeds %" PRIu64
                                      " at step %d\n", limits.work, i+1);
            return false;
        }
        work += (uint64_t)s->rows*unit*(size_t)k;
        if (s->degree) {
            int first = n;
            for (int j = 0; j < s->degree; ++j)
                if (rank[s->scope[j]] < first) first = rank[s->scope[j]];
            s->next = step[first].bucket;
            step[first].bucket = i;
            step[first].factor_count += (size_t)s->degree+1;
        }
    }
    if (limits.stats) fprintf(stderr, "primal DP: width %d, %zu rows, %" PRIu64 " work units\n",
                              width, tables, work);
    return true;
}

static size_t row_index(const Step *s, const int *label, int k)
{
    size_t row = 0;
    for (int j = s->degree-1; j >= 0; --j)
        row = row*(size_t)k+(size_t)label[s->scope[j]];
    return row;
}

static void minimize(int n, int k, const Cost *pair, const Cost *attached,
                     const Cost *total, Step *step, int *label)
{
    for (int i = 0; i < n; ++i) {
        Step *s = &step[i];
        int v = s->v;
        s->choice = mem(s->rows, sizeof(int));
        s->value = mem(s->rows, sizeof(Cost));
        for (size_t row = 0; row < s->rows; ++row) {
            size_t x = row;
            for (int j = 0; j < s->degree; ++j) {
                label[s->scope[j]] = (int)(x%(size_t)k);
                x /= (size_t)k;
            }
            Cost best = INF;
            for (int c = 0; c < k; ++c) {
                label[v] = c;
                Cost cost = total[v]-attached[cell(k, v, c)];
                for (int j = 0; j < s->degree; ++j) {
                    int u = s->scope[j];
                    if (label[u] != c) cost += pair[cell(n, v, u)];
                }
                for (int f = s->bucket; f >= 0; f = step[f].next)
                    cost += step[f].value[row_index(&step[f], label, k)];
                if (cost < best) { best = cost; s->choice[row] = c; }
            }
            s->value[row] = best;
        }
        for (int f = s->bucket; f >= 0; f = step[f].next) {
            free(step[f].value);
            step[f].value = NULL;
        }
    }
    for (int i = n-1; i >= 0; --i)
        label[step[i].v] = step[i].choice[row_index(&step[i], label, k)];
}

static bool solve_piece(const Plane *g, const int *term, int k, Answer *out, Budget budget)
{
    int n = g->n-k;
    if (n > FREE_LIMIT || (n && (size_t)k > budget.rows/(size_t)n)) return false;
    int *fixed = mem((size_t)g->n, sizeof(int)), *id = mem((size_t)g->n, sizeof(int));
    for (int v = 0; v < g->n; ++v) fixed[v] = -1;
    for (int c = 0; c < k; ++c) fixed[term[c]] = c;
    for (int v = 0, u = 0; v < g->n; ++v) if (fixed[v] < 0) id[v] = u++;
    bool *adj = mem(cell(n, n, 0), sizeof(bool));
    for (int e = 0; e < g->m; ++e) {
        int u = g->at[2*e], v = g->at[2*e+1];
        if (u != v && fixed[u] < 0 && fixed[v] < 0)
            adj[cell(n, id[u], id[v])] = adj[cell(n, id[v], id[u])] = true;
    }
    Step *step = mem((size_t)n, sizeof(Step));
    int *rank = mem((size_t)n, sizeof(int));
    bool ok = plan(n, k, adj, step, rank, budget);
    free(rank); free(adj);
    if (ok) {
        Cost *pair = mem(cell(n, n, 0), sizeof(Cost));
        Cost *attached = mem(cell(k, n, 0), sizeof(Cost));
        Cost *total = mem((size_t)n, sizeof(Cost));
        for (int e = 0; e < g->m; ++e) {
            int u = g->at[2*e], v = g->at[2*e+1];
            if (u == v) continue;
            if (fixed[u] >= 0) { int t = u; u = v; v = t; }
            if (fixed[u] >= 0) continue;
            int a = id[u];
            if (fixed[v] >= 0) {
                attached[cell(k, a, fixed[v])] += g->w[e];
                total[a] += g->w[e];
            } else {
                int b = id[v];
                pair[cell(n, a, b)] += g->w[e];
                pair[cell(n, b, a)] += g->w[e];
            }
        }
        int *label = mem((size_t)n, sizeof(int));
        minimize(n, k, pair, attached, total, step, label);
        for (int v = 0; v < g->n; ++v) if (fixed[v] < 0) fixed[v] = label[id[v]];
        *out = (Answer){g->m, 0, mem((size_t)g->m, sizeof(bool))};
        for (int e = 0; e < g->m; ++e)
            if (fixed[g->at[2*e]] != fixed[g->at[2*e+1]]) {
                out->cut[e] = true;
                out->cost += g->w[e];
            }
        free(label); free(pair); free(attached); free(total);
    }
    for (int i = 0; i < n; ++i) {
        free(step[i].scope); free(step[i].value); free(step[i].choice);
    }
    free(step); free(id); free(fixed);
    return ok;
}

/* Persistency of an isolating cut. For any optimum labeling L with
   color-i class A, let S be a minimum source side separating i from all
   other fixed colors. A intersect S is also a feasible source side,
   so w(S minus A, V minus S) <= w(S minus A, A intersect S).
   Recoloring S to i thus cannot increase the objective. Process the
   colors sequentially, respecting every previously fixed vertex. */
static void persist(const Plane *g, const int *term, int k, int *fixed)
{
    Plane h = *g;
    h.at = mem(2*(size_t)g->m, sizeof(int));
    int *map = mem((size_t)g->n, sizeof(int));
    bool *side = mem((size_t)g->n, sizeof(bool));
    for (int i = 0; i < k; ++i) {
        int a = term[i], b = term[i == 0 ? 1 : 0];
        for (int v = 0; v < g->n; ++v)
            map[v] = fixed[v] < 0 ? v : fixed[v] == a ? a : b;
        for (int d = 0; d < 2*g->m; ++d) h.at[d] = map[g->at[d]];
        mincut_max_side(&h, a, b, side);
        for (int v = 0; v < g->n; ++v)
            if (fixed[v] < 0 && side[v]) fixed[v] = a;
    }
    free(side); free(map); free(h.at);
}

/* Removing all fixed vertices splits the objective as well as the graph:
   every free/free or free/terminal edge belongs to one piece, and edges
   between different fixed vertices are mandatory cuts. Within a piece,
   relabel all vertices with nonincident labels to any one incident label:
   this cannot increase any edge cost. A piece touching at most one
   terminal needs no cut. Solve the pieces sequentially to reuse space. */
bool solve_narrow(const Plane *g, const int *term, int k, Answer *out)
{
    size_t n = (size_t)g->n, m = (size_t)g->m;
    Budget limits = budget();
    int *fixed = mem(n, sizeof(int));
    bool *seen = mem(n, sizeof(bool));
    int *vertex = mem(n, sizeof(int)), *boundary = mem(n, sizeof(int));
    int *local = mem(n, sizeof(int)), *edge = mem(m, sizeof(int));
    Answer ans = {g->m, 0, mem(m, sizeof(bool))};
    for (int u = 0; u < g->n; ++u) fixed[u] = local[u] = -1;
    for (int i = 0; i < k; ++i) fixed[term[i]] = term[i];
    persist(g, term, k, fixed);
    for (int e = 0; e < g->m; ++e) {
        int u = g->at[2*e], v = g->at[2*e+1];
        if (fixed[u] >= 0 && fixed[v] >= 0 && fixed[u] != fixed[v]) {
            ans.cut[e] = true;
            ans.cost += g->w[e];
        }
    }

    bool ok = true;
    for (int u = 0; u < g->n && ok; ++u) if (fixed[u] < 0 && !seen[u]) {
        int nv = 1, nt = 0, ne = 0;
        vertex[0] = u; seen[u] = true;
        for (int i = 0; i < nv; ++i) {
            int v = vertex[i], first = g->head[v];
            local[v] = i;
            if (first < 0) continue;
            int d = first;
            do {
                int w = g->at[d ^ 1];
                if (fixed[w] >= 0) {
                    int t = fixed[w];
                    if (local[t] < 0) { local[t] = nt; boundary[nt++] = t; }
                    edge[ne++] = d/2;
                } else if (w != v) {
                    if (!seen[w]) { seen[w] = true; vertex[nv++] = w; }
                    if (!(d & 1)) edge[ne++] = d/2;
                }
                d = g->rot[d];
            } while (d != first);
        }
        if (nt > 1) {
            if (limits.stats)
                fprintf(stderr, "primal DP piece: %d free vertices, %d labels\n", nv, nt);
            Plane h = plane(nv+nt, ne);
            for (int i = 0; i < nt; ++i) local[boundary[i]] += nv;
            for (int i = 0; i < ne; ++i) {
                int e = edge[i];
                for (int j = 0; j < 2; ++j) {
                    int v = g->at[2*e+j];
                    h.at[2*i+j] = local[fixed[v] < 0 ? v : fixed[v]];
                }
                h.w[i] = g->w[e];
            }
            int *reps = mem((size_t)nt, sizeof(int));
            for (int i = 0; i < nt; ++i) reps[i] = nv+i;
            Answer part = {0};
            ok = solve_piece(&h, reps, nt, &part, limits);
            if (ok) {
                ans.cost += part.cost;
                for (int i = 0; i < ne; ++i) ans.cut[edge[i]] = part.cut[i];
                free(part.cut);
            }
            free(reps);
            plane_free(&h);
        }
        for (int i = 0; i < nt; ++i) local[boundary[i]] = -1;
    }
    free(fixed); free(seen); free(vertex); free(boundary); free(local); free(edge);
    if (ok) *out = ans;
    else free(ans.cut);
    return ok;
}
