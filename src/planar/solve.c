/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "cut.h"

/* Contract all other terminals into the sink of an isolating min cut. */
static Cost isolate(const Plane *g, const int *term, int k, int i, bool *cut)
{
    int *vmap = mem((size_t)g->n, sizeof(int));
    for (int v = 0; v < g->n; ++v) vmap[v] = v;
    int sink = term[i == 0 ? 1 : 0];
    for (int j = 0; j < k; ++j) if (j != i) vmap[term[j]] = sink;
    Plane h = *g;
    h.at = mem(2*(size_t)g->m, sizeof(int));
    for (int d = 0; d < 2*g->m; ++d) h.at[d] = vmap[g->at[d]];
    Cost cost = mincut(&h, term[i], sink, cut);
    free(h.at); free(vmap);
    return cost;
}

Answer solve(const Plane *g, const int *term, int k)
{
    size_t m = (size_t)g->m;
    Mask end = 1u << k;
    if (k < 2) return (Answer){g->m, 0, mem(m, sizeof(bool))};
    Answer *dp = mem(end, sizeof(Answer));
    int *tag = mem((size_t)g->n, sizeof(int));
    for (int v = 0; v < g->n; ++v) tag[v] = -1;

    /* Section 3 lower bounds: every (X,Y)-mcut pays the min cut of each
       pair it separates, so a(T) >= min_x lb(x,T-x) + lb(T-x). The
       topology enumeration stops early once this bound is met. */
    Cost *pair = mem((size_t)k*(size_t)k, sizeof(Cost));
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < i; ++j) {
            Cost c = mincut(g, term[i], term[j], NULL);
            pair[(size_t)i*(size_t)k+(size_t)j] = c;
            pair[(size_t)j*(size_t)k+(size_t)i] = c;
        }
    Cost *lb = mem(end, sizeof(Cost));
    for (Mask s = 1; s < end; ++s) {
        if (!(s & (s-1))) continue; /* singletons cost nothing */
        lb[s] = INF;
        for (Mask x = s; x; x = (x-1) & s) {
            Mask rest = s ^ x;
            Cost bound = 0;
            for (int i = 0; i < k; ++i) if (x & (1u << i)) {
                for (int j = 0; j < k; ++j) if ((i != j && (x & (1u << j))) ||
                                                (rest & (1u << j))) {
                    Cost c = pair[(size_t)i*(size_t)k+(size_t)j];
                    if (c > bound) bound = c;
                }
            }
            Cost c = plus(bound, lb[rest]);
            if (c < lb[s]) lb[s] = c;
        }
    }
    free(pair);

    dp[0].cut = mem(m, sizeof(bool));
    for (Mask s = 1; s < end; ++s) {
        dp[s] = (Answer){g->m, 0, mem(m, sizeof(bool))};
        if (!(s & (s-1))) continue;
        int reps[MAX_TERMINALS], q = 0;
        for (int i = 0; i < k; ++i) if (s & (1u << i)) {
            tag[term[i]] = q; reps[q++] = term[i];
        }
        if (q == 2) {
            dp[s].cost = mincut(g, reps[0], reps[1], dp[s].cut);
            tag[reps[0]] = tag[reps[1]] = -1;
            continue;
        }
        Mask count = 1u << q;
        Mask *original = mem(count, sizeof(Mask));
        Answer *best = mem(count, sizeof(Answer));
        original[0] = 0;
        for (int i = 0; i < k; ++i) if (s & (1u << i)) {
            Mask bit = 1u << tag[term[i]];
            for (Mask x = 0; x < bit; ++x) original[x | bit] = original[x] | (1u << i);
        }
        Mask *incident = mem(m, sizeof(Mask));
        for (int j = 0; j < g->m; ++j) {
            int u = g->at[2*j], v = g->at[2*j+1];
            if (u == v) continue;
            if (tag[u] >= 0) incident[j] |= 1u << tag[u];
            if (tag[v] >= 0) incident[j] |= 1u << tag[v];
        }
        for (Mask x = 1; x < count; ++x) {
            best[x] = (Answer){g->m, 0, mem(m, sizeof(bool))};
            for (int j = 0; j < g->m; ++j) if (incident[j] & x) {
                best[x].cut[j] = true; best[x].cost += g->w[j];
            }
        }
        free(incident);
        for (int i = 0; i < q; ++i) tag[reps[i]] = -1;

        /* Each terminal component boundary costs at least its isolating
           cut, and each edge of a multiway cut pays for at most two such
           boundaries. Accumulate the rounded half-sum without overflow. */
        Cost half = 0;
        unsigned odd = 0;
        for (int i = 0; i < q; ++i) {
            Answer *a = &best[1u << i];
            a->cost = isolate(g, reps, q, i, a->cut);
            half += a->cost / 2;
            odd += (unsigned)(a->cost % 2);
        }
        half += (odd + 1) / 2;
        if (half > lb[s]) lb[s] = half;

        /* One cut-open graph per terminal subset. Every realized
           topology can serve all X contained in its isolated labels. */
        Open h;
        open_build(g, reps, q, &h);
        shortest(&h);
        Cost *other = mem(count, sizeof(Cost));
        for (Mask x = 1; x < count; ++x) other[x] = dp[s ^ original[x]].cost;
        enumerate(&h, best, other, lb[s]);
        free(other);
        open_free(&h);

        /* Section 3, Lemma 3.1: choose by sum, then take the union.
           The union is feasible; any overlap is charged only once. */
        Cost min = INF;
        Mask choice = 0;
        for (Mask x = 1; x < count; ++x) {
            Cost c = plus(best[x].cost, dp[s ^ original[x]].cost);
            if (c < min) { min = c; choice = x; }
        }
        if (!choice) die("no feasible terminal partition");
        const bool *a = best[choice].cut, *b = dp[s ^ original[choice]].cut;
        for (int j = 0; j < g->m; ++j) if (a[j] || b[j]) {
            dp[s].cut[j] = true; dp[s].cost += g->w[j];
        }
        for (Mask x = 1; x < count; ++x) free(best[x].cut);
        free(best); free(original);
    }
    Answer ans = dp[end-1];
    for (Mask s = 0; s+1 < end; ++s) free(dp[s].cut);
    free(dp); free(tag); free(lb);
    return ans;
}
