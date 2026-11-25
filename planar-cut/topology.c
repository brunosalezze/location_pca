/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#define _POSIX_C_SOURCE 200809L
#include "cut.h"
#include <time.h>
#include <stdatomic.h>

/* Theorem 4.7 gives a linear bound. Use a conservative constant
   that also covers the sum of Section 5's individual charges. */
enum { CROSSINGS_PER_TERMINAL = 32 };

static unsigned long long stat[9];
static double stat_time;
static int prof = -1;
static bool profiling(void)
{ if (prof < 0) prof = getenv("CUT_PROFILE") != NULL; return prof; }
#define TICK(i) (prof ? (++stat[i]) : 0)
#define PROFILE() (profiling())

typedef struct {
    const Open *h;
    Answer *best;
    Topology t;
    int *ig, *eg, *icap, *labels;
    int *parent, *last, *emap, *eleaf;
    Mask *parity;
    bool *chord_gap, *path_gap;
    int *mapping;
    bool *cut;
    int chords, max_star;
    Mask isolated;
    Cost limit;
    const Cost *other;
    Cost lb;
    atomic_bool *stop;
} Enum;

static int leader(Enum *e, int v, Mask *p)
{
    *p = 0;
    while (e->parent[v] != v) {
        *p ^= e->parity[v];
        v = e->parent[v];
    }
    return v;
}

static void evaluate(Enum *e)
{
    const Open *h = e->h;
    TICK(0);
    struct timespec t0, t1;
    bool on = profiling();
    if (on) clock_gettime(CLOCK_MONOTONIC, &t0);
    Cost c = realize(h, &e->t, e->limit, e->mapping);
    if (on) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        stat_time += (double)(t1.tv_sec-t0.tv_sec) + 1e-9*(double)(t1.tv_nsec-t0.tv_nsec);
    }
    if (c == INF) TICK(8);
    if (c >= e->limit) return;
    memset(e->cut, 0, (size_t)h->input_m * sizeof(bool));
    Cost cost = 0;
    for (int i = e->chords; i < e->t.l; ++i) {
        Link a = e->t.link[i];
        if (a.kind != PATH) continue;
        int u = e->mapping[a.u], v = e->mapping[a.v];
        const int *p = h->first + (size_t)u*(size_t)h->n;
        while (u != v) {
            int d = p[v];
            if (d < 0) die("invalid realization traceback");
            int k = d/2, id = h->id[k];
            if (!e->cut[id]) { e->cut[id] = true; cost += h->w[k]; }
            v = d % 2 ? h->v[k] : h->u[k];
        }
    }
    TICK(7);
    #pragma omp critical (multicut_best)
    {
        Cost ub = INF;
        for (Mask x = e->isolated; x; x = (x-1)&e->isolated)
            if (cost < e->best[x].cost) {
                e->best[x].cost = cost;
                memcpy(e->best[x].cut, e->cut, (size_t)h->input_m*sizeof(bool));
            }
        for (Mask x = 1; x < (Mask)(1u << h->k); ++x) {
            Cost total = plus(e->best[x].cost, e->other[x]);
            if (total < ub) ub = total;
        }
        if (ub <= e->lb) atomic_store_explicit(e->stop, true, memory_order_relaxed);
    }
}

typedef struct Task {
    int from, lo, hi;
    const struct Task *next;
} Task;

/* Root each exterior tree at its first leaf. Enumerate plane binary
   refinements by splitting the remaining leaves into two intervals. */
static void shapes(Enum *e, const Task *task)
{
    TICK(2);
    if (atomic_load_explicit(e->stop, memory_order_relaxed)) return;
    if (!task) { evaluate(e); return; }
    int j = e->t.l++;
    if (task->lo+1 == task->hi) {
        e->t.link[j] = (Link){task->from, e->eleaf[task->lo], PATH};
        shapes(e, task->next);
    } else {
        int v = e->t.n++;
        e->t.link[j] = (Link){task->from, v, PATH};
        for (int mid = task->lo+1; mid < task->hi; ++mid) {
            Task right = {v, mid, task->hi, task->next};
            Task left = {v, task->lo, mid, &right};
            shapes(e, &left);
        }
        --e->t.n;
    }
    --e->t.l;
}

/* A separating cycle is odd on the chords crossing the terminal arc.
   Contract the exterior trees and collect the signatures of the
   fundamental cycles; parity on any arc is linear over GF(2). */
static void valid(Enum *e, int groups)
{
    TICK(1);
    if (atomic_load_explicit(e->stop, memory_order_relaxed)) return;
    for (int i = 0; i < groups; ++i) { e->parent[i] = i; e->parity[i] = 0; }
    Mask separated[MAX_TERMINALS] = {0};
    for (int a = 0; a < e->chords; ++a) {
        Link j = e->t.link[a];
        int u = j.u, v = j.v;
        Mask p, q;
        int i = leader(e, e->eg[u], &p), l = leader(e, e->eg[v], &q);
        Mask weight = p ^ q ^ (1u << e->t.gap[u]) ^ (1u << e->t.gap[v]);
        if (i != l) {
            e->parent[i] = l; e->parity[i] = weight;
        } else if (weight) {
            for (int f = 0; f < e->h->k; ++f) {
                bool odd = false;
                for (int g = f+1; g < e->h->k; ++g) {
                    odd ^= (weight & (1u << (g-1))) != 0;
                    if (odd) {
                        separated[f] |= 1u << g;
                        separated[g] |= 1u << f;
                    }
                }
            }
        }
    }
    Mask iso = 0, all = (1u << e->h->k)-1;
    for (int i = 0; i < e->h->k; ++i)
        if ((separated[i] | (1u << i)) == all) iso |= 1u << e->labels[i];
    if (!iso) return;
    e->isolated = iso;
    e->limit = 0;
    #pragma omp critical (multicut_best)
    {
        for (Mask x = iso; x; x = (x-1)&iso)
            if (e->best[x].cost > e->limit) e->limit = e->best[x].cost;
    }
    if (!e->limit) return;

    Task tasks[groups];
    int j = 0, nt = 0;
    for (int g = 0; g < groups; ++g) {
        int lo = j;
        for (int v = 0; v < e->t.b; ++v)
            if (e->eg[v] == g) e->eleaf[j++] = v;
        if (lo+1 < j) {
            tasks[nt] = (Task){e->eleaf[lo], lo+1, j, nt ? tasks+nt-1 : NULL};
            ++nt;
        }
    }
    shapes(e, nt ? tasks+nt-1 : NULL);
}

/* Appending x to a block can only introduce a crossing involving x
   and the previous occurrence p of that block. */
static bool noncrossing(const int *part, int p, int x)
{
    for (int i = p+1; i < x; ++i)
        for (int j = 0; j < p; ++j)
            if (part[i] == part[j]) return false;
    return true;
}

static void exterior(Enum *e, int x, int groups)
{
    TICK(3);
    if (atomic_load_explicit(e->stop, memory_order_relaxed)) return;
    if (x == e->t.b) { valid(e, groups); return; }
    for (int g = 0; g < groups; ++g) {
        int p = e->last[g];
        int a = e->t.gap[p], b = e->t.gap[x];
        if (!e->path_gap[a*e->h->k+b] || !noncrossing(e->eg, p, x)) continue;
        e->eg[x] = g; e->last[g] = x;
        exterior(e, x+1, groups);
        e->last[g] = p;
    }
    e->eg[x] = groups; e->last[groups] = x;
    exterior(e, x+1, groups+1);
}

static void stars(Enum *e, int g, int groups)
{
    TICK(4);
    if (atomic_load_explicit(e->stop, memory_order_relaxed)) return;
    if (g == groups) {
        e->chords = e->t.l;
        /* Sec. 6: consecutive constrained positions within each gap. */
        for (int i = 1; i < e->t.b; ++i)
            if (e->t.gap[i-1] == e->t.gap[i])
                e->t.link[e->t.l++] = (Link){i-1, i, ORDER};
        exterior(e, 0, 0);
        e->t.l = e->chords;
        return;
    }
    int count = 0, first = -1;
    for (int i = 0; i < e->t.b; ++i) if (e->ig[i] == g) {
        if (count == 0) first = i;
        ++count;
    }
    if (count < 2) return;
    for (int c = first; c < e->t.b; ++c) if (e->ig[c] == g) {
        int old = e->t.l;
        bool possible = true;
        for (int i = first; i < e->t.b; ++i) if (i != c && e->ig[i] == g) {
            int a = e->t.gap[c], b = e->t.gap[i];
            if (!e->chord_gap[a*e->h->k+b]) { possible = false; break; }
            e->t.link[e->t.l++] = (Link){c, i, CHORD};
        }
        if (possible) stars(e, g+1, groups);
        e->t.l = old;
        if (count == 2) break; /* a single edge has no distinguished end */
    }
}

static void interior(Enum *e, int x, int groups, int singles)
{
    TICK(5);
    if (atomic_load_explicit(e->stop, memory_order_relaxed)) return;
    if (x-groups > g_crossings*e->h->k) return;
    if (singles > e->t.b-x) return;
    if (x == e->t.b) { stars(e, 0, groups); return; }
    for (int g = 0; g < groups; ++g) {
        int p = e->icap[g];
        if (!noncrossing(e->ig, p, x)) continue;
        int size = 0;
        for (int i = 0; i <= p; ++i) size += e->ig[i] == g;
        if (size == e->max_star) continue;
        e->ig[x] = g; e->icap[g] = x;
        interior(e, x+1, groups, singles-(size == 1));
        e->icap[g] = p;
    }
    if (groups < e->t.b/2) {
        e->ig[x] = groups; e->icap[groups] = x;
        interior(e, x+1, groups+1, singles+1);
    }
}

static void gaps(Enum *e, const int *count, int i, int left, int pos)
{
    TICK(6);
    if (atomic_load_explicit(e->stop, memory_order_relaxed)) return;
    if (i == e->h->k) {
        if (!left) interior(e, 0, 0, 0);
        return;
    }
    if (left > count[i]) return;
    for (int j = 0; j <= left && j <= count[i]-count[i+1]; ++j) {
        if (j) e->t.gap[pos+j-1] = i;
        gaps(e, count, i+1, left-j, pos+j);
    }
}

void enumerate(const Open *h, Answer *best, const Cost *other, Cost lb)
{
    int n = 0, spokes = 0, count[MAX_TERMINALS+1] = {0}, labels[MAX_TERMINALS];
    int next = 0, max_star = 0;
    int *size = mem((size_t)h->b, sizeof(int));
    for (int i = 0; i < h->b; ++i)
        if (++size[h->center[i]] > max_star) max_star = size[h->center[i]];
    for (int i = 0; i < h->b; ++i) {
        if (h->label[i] >= 0) labels[next++] = h->label[i];
        else {
            spokes += h->center[i] != i;
            if (size[h->center[i]] > 1) { ++n; ++count[h->gap[i]]; }
        }
    }
    free(size);
    if (n > 2*g_crossings*h->k) n = 2*g_crossings*h->k;
    if (n > 2*spokes) n = 2*spokes;
    if (n < 2) return;
    for (int i = h->k-1; i >= 0; --i) count[i] += count[i+1];

    Cost ub = INF;
    for (Mask x = 1; x < (Mask)(1u << h->k); ++x) {
        Cost c = plus(best[x].cost, other[x]);
        if (c < ub) ub = c;
    }
    if (ub <= lb) return;

    size_t s = (size_t)n, kk = (size_t)h->k*(size_t)h->k;
    bool *chord_gap = mem(kk, sizeof(bool)), *path_gap = mem(kk, sizeof(bool));
    for (int u = 0; u < h->b; ++u) if (h->gap[u] >= 0) {
        for (int v = 0; v < h->b; ++v) if (h->gap[v] >= 0) {
            int a = h->gap[u]*h->k+h->gap[v];
            if (u != v && (h->center[u] == v || h->center[v] == u)) chord_gap[a] = true;
            if (h->dist[(size_t)u*(size_t)h->n+(size_t)v] < INF) path_gap[a] = true;
        }
    }

    /* One task per (port count, ports in gap 0, ports in gap 1); the
       rest of the gap distribution recurses inside the task. Tasks run
       in parallel and share best[] through a critical section. */
    int cap0 = count[0]-count[1], cap1 = count[1]-count[2], total = 0;
    for (int i = 2; i <= n; ++i)
        for (int j0 = 0; j0 <= i && j0 <= cap0; ++j0)
            total += (i-j0 < cap1 ? i-j0 : cap1)+1;
    int *tasks = mem((size_t)total*3, sizeof(int));
    for (int i = 2, tk = 0; i <= n; ++i)
        for (int j0 = 0; j0 <= i && j0 <= cap0; ++j0)
            for (int j1 = 0; j1 <= i-j0 && j1 <= cap1; ++j1) {
                tasks[3*tk] = i; tasks[3*tk+1] = j0; tasks[3*tk+2] = j1; ++tk;
            }

    atomic_bool stop = false;
    (void)profiling();
    profile_init();
    #pragma omp parallel for schedule(dynamic)
    for (int tk = 0; tk < total; ++tk) {
        if (stop) continue;
        Enum e = {
            .h = h, .best = best, .labels = labels,
            .other = other, .lb = lb, .max_star = max_star, .stop = &stop,
            .t = {.gap = mem(s, sizeof(int)), .link = mem(4*s, sizeof(Link))},
            .ig = mem(s, sizeof(int)), .eg = mem(s, sizeof(int)),
            .icap = mem(s, sizeof(int)), .last = mem(s, sizeof(int)),
            .parent = mem(s, sizeof(int)), .parity = mem(s, sizeof(Mask)),
            .eleaf = mem(s, sizeof(int)), .mapping = mem(2*s, sizeof(int)),
            .cut = mem((size_t)h->input_m, sizeof(bool)),
            .chord_gap = chord_gap, .path_gap = path_gap
        };
        e.t.eg = e.eg;
        int i = tasks[3*tk], j0 = tasks[3*tk+1], j1 = tasks[3*tk+2];
        e.t.b = e.t.n = i;
        for (int p = 0; p < j0; ++p) e.t.gap[p] = 0;
        for (int p = j0; p < j0+j1; ++p) e.t.gap[p] = 1;
        gaps(&e, count, 2, i-j0-j1, j0+j1);
        free(e.t.gap); free(e.t.link); free(e.ig); free(e.eg); free(e.icap);
        free(e.last); free(e.parent); free(e.parity); free(e.eleaf); free(e.mapping);
        free(e.cut);
    }
    if (PROFILE()) fprintf(stderr, "gaps=%llu interior=%llu stars=%llu exterior=%llu "
                            "valid=%llu shapes=%llu evaluate=%llu realize_time=%.2f\n",
                            stat[6], stat[5], stat[4], stat[3], stat[1], stat[2], stat[0], stat_time);
    if (PROFILE()) fprintf(stderr, "improved=%llu infeasible=%llu stop=%d ub=%llu lb=%llu\n",
                           stat[7], stat[8], stop, (unsigned long long)ub, (unsigned long long)lb);
    free(tasks); free(chord_gap); free(path_gap);
}
