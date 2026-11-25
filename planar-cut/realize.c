/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#define _POSIX_C_SOURCE 200809L
#include "cut.h"
#include <time.h>
#include <stdatomic.h>
static double rtime[4];
static unsigned long long rcount[8];
static int prof = -1;
static bool profiling(void)
{ if (prof < 0) prof = getenv("CUT_PROFILE") != NULL; return prof; }
void profile_init(void) { (void)profiling(); }
static double elapsed(const struct timespec *a, const struct timespec *b)
{ return (double)(b->tv_sec-a->tv_sec) + 1e-9*(double)(b->tv_nsec-a->tv_nsec); }
#define PROF(i, code) do { struct timespec t0_, t1_; \
    clock_gettime(CLOCK_MONOTONIC, &t0_); code; clock_gettime(CLOCK_MONOTONIC, &t1_); \
    rtime[i] += elapsed(&t0_, &t1_); } while (0)

typedef struct {
    int size, *values;
} Domain;

typedef struct {
    int n, *scope, *witness;
    Cost *value;
    const Link *link;
    bool active;
} Factor;

static Cost link_cost(const Open *h, const Link *l, int u, int v)
{
    switch (l->kind) {
        case CHORD: return u != v && (h->center[u] == v || h->center[v] == u) ? 0 : INF;
        case ORDER: return u < v ? 0 : INF;
        case PATH: return h->dist[(size_t)u*(size_t)h->n+(size_t)v];
    }
    die("unknown constraint");
    return INF;
}

/* Test all components of the surviving induced graph after a proposed
   separator; found[] is reused to return their component numbers. */
static bool balanced(int n, const bool *a, const int *vertices, int count,
                     const bool *removed, int *found)
{
    int q[count];
    for (int i = 0; i < count; ++i) found[i] = -1;
    for (int i = 0; i < count; ++i) if (!removed[i] && found[i] < 0) {
        int size = 1;
        q[0] = i; found[i] = i;
        for (int j = 0; j < size; ++j) {
            const bool *row = a + (size_t)vertices[q[j]]*(size_t)n;
            for (int u = 0; u < count; ++u)
                if (!removed[u] && found[u] < 0 && row[vertices[u]]) {
                    found[u] = i; q[size++] = u;
                }
        }
        if (3*size > 2*count) return false;
    }
    return true;
}

static bool find_separator(int n, const bool *a, const int *v, int count,
                           bool *removed, int *found, int start, int remaining)
{
    if (!remaining) return balanced(n, a, v, count, removed, found);
    for (int i = start; i <= count-remaining; ++i) {
        removed[i] = true;
        if (find_separator(n, a, v, count, removed, found, i+1, remaining-1)) return true;
        removed[i] = false;
    }
    return false;
}

static int ceil_sqrt(int n)
{
    int r = 0;
    while (r*r < n) ++r;
    return r;
}

/* Planar separator fallback certifies the O(sqrt(n)) width; in the
   common case min-fill below finds a smaller order far more cheaply. */
static void separator_order(int n, const bool *a, const int *v, int count, int *order)
{
    if (count <= 2) {
        memcpy(order, v, (size_t)count*sizeof(int));
        return;
    }
    int found[count], part[count], s = 0, at = 0;
    bool removed[count];
    memset(removed, 0, sizeof removed);
    while (!find_separator(n, a, v, count, removed, found, 0, s)) {
        if (++s > ceil_sqrt(8*count)) die("topology constraint graph is not planar");
    }
    for (int i = 0; i < count; ++i) if (found[i] == i) {
        int size = 0;
        for (int j = i; j < count; ++j) if (found[j] == i) part[size++] = v[j];
        separator_order(n, a, part, size, order+at);
        at += size;
    }
    for (int i = 0; i < count; ++i) if (removed[i]) order[at++] = v[i];
}

static void elimination(const Topology *t, int *order, const Domain *d)
{
    int n = t->n, width = 0;
    size_t cells = (size_t)n*(size_t)n;
    bool *original = mem(cells, sizeof(bool)), *a = mem(cells, sizeof(bool));
    bool used[n];
    int neighbors[n];
    memset(used, 0, sizeof used);
    for (int i = 0; i < t->l; ++i) {
        int u = t->link[i].u, v = t->link[i].v;
        original[(size_t)u*(size_t)n+(size_t)v] = true;
        original[(size_t)v*(size_t)n+(size_t)u] = true;
    }
    memcpy(a, original, cells*sizeof(bool));
    for (int i = 0; i < n; ++i) {
        int best = -1, lo = n*n, degree = n;
        for (int u = 0; u < n; ++u) if (!used[u]) {
            const bool *row = a + (size_t)u*(size_t)n;
            int fill = 0, size = 0;
            for (int v = 0; v < n; ++v) if (!used[v] && row[v]) neighbors[size++] = v;
            for (int j = 0; j < size; ++j)
                for (int k = 0; k < j; ++k)
                    fill += !a[(size_t)neighbors[j]*(size_t)n+(size_t)neighbors[k]];
            if (fill < lo || (fill == lo && (size < degree ||
                    (size == degree && d[u].size > d[best].size)))) {
                best = u; lo = fill; degree = size;
            }
        }
        if (degree > width) width = degree;
        order[i] = best; used[best] = true;
        const bool *row = a+(size_t)best*(size_t)n;
        for (int v = 0; v < n; ++v) if (!used[v] && row[v])
            for (int u = v+1; u < n; ++u) if (!used[u] && row[u]) {
                a[(size_t)v*(size_t)n+(size_t)u] = true;
                a[(size_t)u*(size_t)n+(size_t)v] = true;
            }
    }
    if (width > 8*ceil_sqrt(n)) {
        for (int i = 0; i < n; ++i) neighbors[i] = i;
        separator_order(n, original, neighbors, n, order);
    }
    free(a); free(original);
}

static bool contains(const Factor *f, int v)
{
    for (int j = 0; j < f->n; ++j) if (f->scope[j] == v) return true;
    return false;
}

static Cost lookup(const Open *h, const Factor *f, const Domain *d, const int *x)
{
    if (f->link) {
        int u = f->link->u, v = f->link->v;
        return link_cost(h, f->link, d[u].values[x[u]], d[v].values[x[v]]);
    }
    size_t i = 0;
    for (int j = 0; j < f->n; ++j) {
        int u = f->scope[j];
        i = i*(size_t)d[u].size+(size_t)x[u];
    }
    return f->value[i];
}

static void decode(const Factor *f, const Domain *d, int *x, size_t i)
{
    for (int j = f->n; j > 0; --j) {
        int u = f->scope[j-1];
        x[u] = (int)(i % (size_t)d[u].size);
        i /= (size_t)d[u].size;
    }
}

/* Fast feasibility of the chord, order and port-path links over the
   boundary positions. Almost every topology is infeasible here;
   rejecting it before the elimination DP keeps the enumeration
   affordable. If the node budget runs out, fall through to the exact
   DP instead of guessing. */
static bool ports_rec(const Open *h, const Topology *t, int *x, int u, long *budget);

static bool ports_ok(const Open *h, const Topology *t, int *x, long budget)
{
    return ports_rec(h, t, x, 0, &budget);
}

/* Cheapest cost over port assignments of the exterior trees: any tree
   connecting a group's ports pays at least the largest pairwise
   distance between them (exact for two ports). Reject the topology
   before the elimination DP when even that reaches the limit. */
static bool cheap_rec(const Open *h, const Topology *t, int *x, int u,
                      Cost *high, Cost partial, Cost limit, long *budget)
{
    if (u == t->b) {
        /* Full assignment: each exterior tree is at least half its
           metric MST, and at least its widest pair. */
        Cost total = 0;
        for (int g = 0; g <= t->b; ++g) {
            Cost bound = high[g];
            if (!bound) continue;
            int members[t->b], m = 0;
            for (int w = 0; w < t->b; ++w) if (t->eg[w] == g) members[m++] = x[w];
            Cost span[t->b];
            for (int a = 1; a < m; ++a)
                span[a] = h->dist[(size_t)members[0]*(size_t)h->n+(size_t)members[a]];
            span[0] = 0;
            Cost mst = 0;
            for (int a = 1; a < m; ++a) {
                int near = -1;
                for (int b2 = 0; b2 < m; ++b2) if (span[b2] && (near < 0 || span[b2] < span[near])) near = b2;
                if (near < 0) break;
                mst += span[near];
                span[near] = 0;
                for (int b2 = 0; b2 < m; ++b2) if (span[b2]) {
                    Cost d = h->dist[(size_t)members[near]*(size_t)h->n+(size_t)members[b2]];
                    if (d < span[b2]) span[b2] = d;
                }
            }
            if (mst/2 > bound) bound = mst/2;
            total = plus(total, bound);
            if (total >= limit) return false;
        }
        return true;
    }
    if (--*budget <= 0) return true;
    int g = t->gap[u], at = 0;
    for (int i = 0; i < g; ++i) at += h->dom_n[i];
    for (int i = 0; i < h->dom_n[g]; ++i) {
        int v = h->dom[at+i], ok = 1;
        Cost run = partial, old = high[t->eg[u]];
        for (int j = 0; j < t->l && ok; ++j) {
            const Link *l = &t->link[j];
            int other = l->u == u ? l->v : (l->v == u ? l->u : -1);
            if (other < 0 || other >= u) continue;
            int a = l->u == u ? v : x[other], b = l->u == u ? x[other] : v;
            switch (l->kind) {
                case CHORD: ok = a != b && (h->center[a] == b || h->center[b] == a); break;
                case ORDER: ok = a < b; break;
                case PATH: ok = h->dist[(size_t)a*(size_t)h->n+(size_t)b] < INF; break;
            }
        }
        if (!ok) continue;
        int grp = t->eg[u];
        for (int w = 0; w < u; ++w) if (t->eg[w] == grp) {
            Cost d = h->dist[(size_t)v*(size_t)h->n+(size_t)x[w]];
            if (d < INF && d > high[grp]) high[grp] = d;
        }
        run += high[grp]-old;
        if (run >= limit) { high[grp] = old; continue; }
        x[u] = v;
        if (cheap_rec(h, t, x, u+1, high, run, limit, budget)) return true;
        high[grp] = old;
    }
    x[u] = -1;
    return false;
}

static bool cheap_cost(const Open *h, const Topology *t, int *x, Cost limit)
{
    Cost high[t->b+1];
    for (int i = 0; i <= t->b; ++i) high[i] = 0;
    long budget = 20000;
    return cheap_rec(h, t, x, 0, high, 0, limit, &budget);
}

static bool ports_rec(const Open *h, const Topology *t, int *x, int u, long *budget)
{
    if (u == t->b) return true;
    if (--*budget <= 0) return true;
    int g = t->gap[u], at = 0;
    for (int i = 0; i < g; ++i) at += h->dom_n[i];
    for (int i = 0; i < h->dom_n[g]; ++i) {
        int v = h->dom[at+i], ok = 1;
        for (int j = 0; j < t->l && ok; ++j) {
            const Link *l = &t->link[j];
            int other = l->u == u ? l->v : (l->v == u ? l->u : -1);
            if (other < 0 || other >= u) continue;
            int a = l->u == u ? v : x[other], b = l->u == u ? x[other] : v;
            switch (l->kind) {
                case CHORD: ok = a != b && (h->center[a] == b || h->center[b] == a); break;
                case ORDER: ok = a < b; break;
                case PATH: ok = h->dist[(size_t)a*(size_t)h->n+(size_t)b] < INF; break;
            }
        }
        if (!ok) continue;
        x[u] = v;
        if (ports_rec(h, t, x, u+1, budget)) return true;
    }
    x[u] = -1;
    return false;
}

/* Arc consistency is only a pruning step, not a feasibility test. */
static bool prune(const Open *h, const Topology *t, Domain *d, Cost limit)
{
    bool changed;
    do {
        changed = false;
        for (int j = 0; j < t->l; ++j) {
            const Link *l = &t->link[j];
            /* Unconstrained branch vertices can have large domains.
               Skip expensive PATH pruning there; the DP is exact anyway. */
            if (l->kind == PATH && (l->u >= t->b || l->v >= t->b)) continue;
            for (int rev = 0; rev < 2; ++rev) {
                Domain *a = &d[rev ? l->v : l->u];
                const Domain *b = &d[rev ? l->u : l->v];
                int count = 0;
                for (int v = 0; v < a->size; ++v) {
                    for (int u = 0; u < b->size; ++u) {
                        int x = a->values[v], y = b->values[u];
                        if (link_cost(h, l, rev ? y : x, rev ? x : y) < limit) {
                            a->values[count++] = x;
                            break;
                        }
                    }
                }
                if (!count) return false;
                if (count != a->size) { changed = true; a->size = count; }
            }
        }
    } while (changed);
    return true;
}

/* Section 6: evaluate a vertex-elimination tree decomposition. Each
   message is indexed by the bag minus its eliminated vertex; keep
   argmins to reconstruct the topology embedding. */
Cost realize(const Open *h, const Topology *t, Cost limit, int *map)
{
    int n = t->n, l = t->l;
    int need = 0;
    for (int u = 0; u < t->b; ++u) need += h->dom_n[t->gap[u]];
    int *values = mem((size_t)(need ? need : 1), sizeof(int));
    int *x = mem((size_t)n, sizeof(int));
    Domain d[n];
    int at = 0;
    for (int u = 0; u < n; ++u) {
        if (u >= t->b) {
            d[u] = (Domain){h->n, h->all}; /* never rewritten by prune */
        } else {
            int g = t->gap[u], sz = h->dom_n[g], base = 0;
            for (int i = 0; i < g; ++i) base += h->dom_n[i];
            d[u] = (Domain){sz, values+at};
            memcpy(d[u].values, h->dom+base, (size_t)sz*sizeof(int));
            at += sz;
        }
    }
    for (int u = 0; u < n; ++u) x[u] = -1;
    bool on = profiling();
    struct timespec s0 = {0}, s1 = {0}, s2 = {0}, s3 = {0};
    if (on) clock_gettime(CLOCK_MONOTONIC, &s0);
    bool ok1 = ports_ok(h, t, x, 2000);
    bool ok1b = ok1 && cheap_cost(h, t, x, limit);
    if (on) clock_gettime(CLOCK_MONOTONIC, &s1);
    bool ok2 = ok1b && prune(h, t, d, limit);
    if (on) clock_gettime(CLOCK_MONOTONIC, &s2);
    if (on) { rtime[0] += elapsed(&s0, &s1); rtime[1] += elapsed(&s1, &s2); }
    if (!ok2) { if (on) ++rcount[ok1b ? 1 : 0]; free(values); free(x); return INF; }
    Factor *f = mem((size_t)(l+n), sizeof(Factor));
    int order[n], pulled[l+n], scope[n], ends[l][2];
    bool mark[n];
    elimination(t, order, d);
    for (int i = 0; i < l; ++i) {
        ends[i][0] = t->link[i].u; ends[i][1] = t->link[i].v;
        f[i] = (Factor){.n = 2, .scope = ends[i], .link = &t->link[i], .active = true};
    }
    Cost total = 0;
    for (int step = 0; step < n; ++step) {
        int v = order[step], count = 0, size = 0, idx = l+step;
        Factor *g = &f[idx];
        memset(mark, 0, sizeof mark);
        mark[v] = true;
        for (int i = 0; i < idx; ++i) if (f[i].active && contains(&f[i], v)) {
            pulled[count++] = i; f[i].active = false;
            for (int j = 0; j < f[i].n; ++j) {
                int u = f[i].scope[j];
                if (!mark[u]) { mark[u] = true; scope[size++] = u; }
            }
        }
        g->n = size; g->scope = mem((size_t)size, sizeof(int));
        memcpy(g->scope, scope, (size_t)size*sizeof(int));
        size_t cells = 1;
        for (int j = 0; j < size; ++j) {
            size_t s = (size_t)d[scope[j]].size;
            if (cells > SIZE_MAX/s) die("realization table too large");
            cells *= s;
        }
        g->value = mem(cells, sizeof(Cost));
        g->witness = mem(cells, sizeof(int));
        g->active = size > 0;
        bool finite = false;
        for (size_t i = 0; i < cells; ++i) {
            decode(g, d, x, i);
            Cost best = limit;
            for (int j = 0; j < d[v].size; ++j) {
                x[v] = j;
                Cost sum = 0;
                for (int k = 0; k < count && sum < best; ++k)
                    sum = plus(sum, lookup(h, &f[pulled[k]], d, x));
                if (sum < best) { best = sum; g->witness[i] = j; }
            }
            g->value[i] = best;
            finite |= best < limit;
        }
        for (int i = 0; i < count; ++i) {
            Factor *old = &f[pulled[i]];
            free(old->value); old->value = NULL;
        }
        if (!finite) { total = INF; break; }
        if (!size) total = plus(total, g->value[0]);
        if (total >= limit) break;
    }
    if (on) { clock_gettime(CLOCK_MONOTONIC, &s3); rtime[2] += elapsed(&s2, &s3); }
    if (on) ++rcount[total == INF ? 3 : 4];
    if (total < limit) {
        for (int step = n; step > 0; --step) {
            int u = order[step-1];
            Factor *g = &f[l+step-1];
            size_t i = 0;
            for (int j = 0; j < g->n; ++j) {
                int v = g->scope[j];
                i = i*(size_t)d[v].size+(size_t)x[v];
            }
            x[u] = g->witness[i];
            map[u] = d[u].values[x[u]];
        }
    }
    for (int i = l; i < l+n; ++i) { free(f[i].scope); free(f[i].witness); free(f[i].value); }
    free(values); free(x); free(f);
    if (prof) fprintf(stderr, "  realize: ok_ports=%.2fs ok_prune=%.2fs dp=%.2fs | "
        "rej_ports=%llu rej_prune=%llu pass=%llu dp_inf=%llu dp_fin=%llu\n",
        rtime[0], rtime[1], rtime[2], rcount[0], rcount[1], rcount[2], rcount[3], rcount[4]);
    return total;
}
