/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "cut.h"
#include <limits.h>

void die(const char *msg)
{
    fprintf(stderr, "location_pca: %s\n", msg);
    exit(EXIT_FAILURE);
}

/* Undirected s–t min cut (Dinic, uint64 capacities). */
static Cost flow_cut(const Plane *g, int s, int t, bool *cut, bool *side)
{
    int n = g->n;
    size_t arcs = 2 * (size_t)g->m;
    int *to = mem(2*arcs, sizeof(int)), *next = mem(2*arcs, sizeof(int));
    int *head = mem((size_t)n, sizeof(int)), *dist = mem((size_t)n, sizeof(int));
    int *it = mem((size_t)n, sizeof(int)), *queue = mem((size_t)n, sizeof(int));
    Cost *cap = mem(2*arcs, sizeof(Cost));
    for (int v = 0; v < n; ++v) head[v] = -1;
    for (int e = 0; e < g->m; ++e) {
        int u = g->at[2*e], v = g->at[2*e+1], a = 2*e, b = 2*e+1;
        to[2*a] = v; cap[2*a] = g->w[e]; next[2*a] = head[u]; head[u] = 2*a;
        to[2*a+1] = u; cap[2*a+1] = 0; next[2*a+1] = head[v]; head[v] = 2*a+1;
        to[2*b] = u; cap[2*b] = g->w[e]; next[2*b] = head[v]; head[v] = 2*b;
        to[2*b+1] = v; cap[2*b+1] = 0; next[2*b+1] = head[u]; head[u] = 2*b+1;
    }
    Cost flow = 0;
    int *path = mem((size_t)n, sizeof(int)), *pv = mem((size_t)n, sizeof(int));
    for (;;) {
        for (int v = 0; v < n; ++v) dist[v] = -1;
        dist[s] = 0;
        int front = 0, back = 0;
        queue[back++] = s;
        while (front < back) {
            int u = queue[front++];
            for (int d = head[u]; d >= 0; d = next[d])
                if (cap[d] && dist[to[d]] < 0) { dist[to[d]] = dist[u]+1; queue[back++] = to[d]; }
        }
        if (dist[t] < 0) break;
        for (int v = 0; v < n; ++v) it[v] = head[v];
        int u = s, len = 0;
        for (;;) {
            if (u == t) {
                Cost push = INF;
                for (int i = 0; i < len; ++i)
                    if (cap[path[i]] < push) push = cap[path[i]];
                for (int i = 0; i < len; ++i) {
                    cap[path[i]] -= push;
                    cap[path[i] ^ 1] += push;
                }
                flow = plus(flow, push);
                if (flow == INF) die("min-cut value overflow");
                u = s; len = 0;
                continue;
            }
            int d = it[u];
            while (d >= 0 && !(cap[d] && dist[to[d]] == dist[u]+1)) d = next[d];
            it[u] = d;
            if (d >= 0) {
                path[len] = d; pv[len++] = u;
                u = to[d];
            } else {
                dist[u] = -1;
                if (!len) break;
                u = pv[--len];
            }
        }
    }
    if (cut) for (int e = 0; e < g->m; ++e)
        cut[e] = (dist[g->at[2*e]] >= 0) != (dist[g->at[2*e+1]] >= 0);
    if (side) {
        /* The complement of vertices that can reach the sink in the
           residual graph is the largest minimum source side. */
        for (int v = 0; v < n; ++v) side[v] = true;
        int front = 0, back = 0;
        queue[back++] = t; side[t] = false;
        while (front < back) {
            int v = queue[front++];
            for (int d = head[v]; d >= 0; d = next[d]) if (cap[d ^ 1] && side[to[d]]) {
                side[to[d]] = false; queue[back++] = to[d];
            }
        }
    }
    free(path); free(pv);
    free(to); free(next); free(head); free(dist); free(it); free(queue); free(cap);
    return flow;
}

Cost mincut(const Plane *g, int s, int t, bool *cut)
{
    return flow_cut(g, s, t, cut, NULL);
}

void mincut_max_side(const Plane *g, int s, int t, bool *side)
{
    (void)flow_cut(g, s, t, NULL, side);
}

void *mem(size_t n, size_t s)
{
    if (s && n > SIZE_MAX / s) die("allocation size overflow");
    void *p = calloc(n ? n : 1, s);
    if (!p) die("out of memory");
    return p;
}

Cost plus(Cost a, Cost b)
{
    return b >= INF - a ? INF : a + b;
}

Plane plane(int n, int m)
{
    if (m > INT_MAX / 2) die("too many edges");
    size_t nn = (size_t)n, mm = (size_t)m;
    Plane g = {n, m, mem(2 * mm, sizeof(int)), mem(2 * mm, sizeof(int)),
               mem(nn, sizeof(int)), mem(mm, sizeof(Cost))};
    for (int i = 0; i < n; ++i) g.head[i] = -1;
    for (int i = 0; i < 2*m; ++i) g.at[i] = g.rot[i] = -1;
    return g;
}

void plane_free(Plane *g)
{
    free(g->at); free(g->rot); free(g->head); free(g->w);
}

int *faces(const Plane *g, int *nf)
{
    int *f = mem(2*(size_t)g->m, sizeof(int));
    for (int d = 0; d < 2*g->m; ++d) f[d] = -1;
    *nf = g->m == 0;
    for (int d = 0; d < 2*g->m; ++d) if (f[d] < 0) {
        int a = d;
        do {
            f[a] = *nf;
            a = g->rot[a ^ 1];
        } while (a != d);
        ++*nf;
    }
    return f;
}

/* B_SOLVE lines 0–1. Port i gets old dart i and both cycle darts.
   This also implements the one-port loop and two-port digon. */
static Plane expand(const Plane *g, const int *term, int k, int *spoke)
{
    int *port = mem(2*(size_t)g->m, sizeof(int));
    int *rim = mem(2*(size_t)g->m, sizeof(int));
    int nn = g->n, mm = g->m;
    for (int d = 0; d < 2*g->m; ++d) {
        port[d] = g->at[d]; rim[d] = -1;
    }
    for (int i = 0; i < k; ++i) {
        int s = g->head[term[i]], d = s;
        if (s < 0) die("isolated terminal in connected component");
        spoke[i] = mm;
        do {
            if (d != s) {
                if (nn == INT_MAX) die("too many expanded vertices");
                port[d] = nn++;
            }
            if (mm == INT_MAX / 2) die("too many expanded edges");
            rim[d] = mm++;
            d = g->rot[d];
        } while (d != s);
    }
    Plane p = plane(nn, mm);
    memcpy(p.w, g->w, (size_t)g->m * sizeof(Cost));
    for (int d = 0; d < 2*g->m; ++d) {
        p.at[d] = port[d]; p.head[port[d]] = d;
        if (rim[d] < 0) {
            p.rot[d] = g->rot[d];
        } else {
            int next = g->rot[d], a = 2*rim[d], rev = a+1;
            p.at[a] = port[d]; p.at[rev] = port[next];
            p.rot[d] = a;
            p.rot[rev] = next;
            p.w[rim[d]] = INF;
        }
    }
    for (int d = 0; d < 2*g->m; ++d)
        if (rim[d] >= 0) p.rot[2*rim[g->rot[d]]] = 2*rim[d]+1;
    free(port); free(rim);
    return p;
}

/* Implemented in steiner.c. */
void steiner(const Plane *dual, int old_m, const int *rep, int k, bool *tree);

void open_build(const Plane *g, const int *term, int k, Open *h)
{
    int *rep = mem((size_t)k, sizeof(int));
    Plane p = expand(g, term, k, rep);
    int nf;
    int *f = faces(&p, &nf);
    if ((int64_t)p.n-p.m+nf != 2) die("terminal expansion lost the planar embedding");
    Plane d = plane(nf, p.m);
    for (int a = 0; a < 2*p.m; ++a) {
        d.at[a] = f[a]; d.rot[a] = p.rot[a ^ 1]; d.head[f[a]] = a;
    }
    memcpy(d.w, p.w, (size_t)p.m * sizeof(Cost));
    for (int i = 0; i < k; ++i) rep[i] = f[2*rep[i]+1];
    free(f); plane_free(&p);

    bool *tree = mem((size_t)d.m, sizeof(bool));
    steiner(&d, g->m, rep, k, tree);
    int *corner = mem(2*(size_t)d.m, sizeof(int));
    int *out = mem(2*(size_t)d.m, sizeof(int));
    int *root = mem((size_t)d.n, sizeof(int));
    int *tag = mem((size_t)d.n, sizeof(int));
    for (int v = 0; v < d.n; ++v) root[v] = tag[v] = -1;
    for (int i = 0; i < k; ++i) tag[rep[i]] = i;

    int start = d.head[rep[0]];
    while (!tree[start / 2]) start = d.rot[start];
    int b = 0, a = start;
    do {
        if (b == 2*d.m) die("invalid Steiner contour");
        out[b] = a;
        corner[a] = b++;
        int next = d.rot[a ^ 1];
        while (!tree[next / 2]) next = d.rot[next];
        for (int x = d.rot[a ^ 1]; x != next; x = d.rot[x])
            corner[x] = b; /* next corner, with wraparound corrected below */
        a = next;
    } while (a != start);
    int ntree = 0;
    for (int e = 0; e < d.m; ++e) ntree += tree[e];
    if (b != 2*ntree) die("Steiner contour does not cover the tree");
    for (int x = d.rot[out[b-1] ^ 1]; x != start; x = d.rot[x]) corner[x] = 0;

    *h = (Open){.n = b, .b = b, .k = k, .input_m = g->m};
    h->label = mem((size_t)b, sizeof(int));
    h->gap = mem((size_t)b, sizeof(int));
    h->center = mem((size_t)b, sizeof(int));
    int gap = -1;
    for (int i = 0; i < b; ++i) {
        int v = d.at[out[i]];
        h->label[i] = tag[v];
        if (tag[v] >= 0) ++gap;
        h->gap[i] = tag[v] >= 0 ? -1 : gap;
        if (root[v] < 0) root[v] = i;
        h->center[i] = root[v];
    }
    /* Non-tree vertices have a single unsplit copy. */
    for (int v = 0; v < d.n; ++v) if (root[v] < 0) {
        if (h->n == INT_MAX) die("too many cut-open vertices");
        root[v] = h->n++;
        if (d.head[v] >= 0) {
            int x = d.head[v];
            do { corner[x] = root[v]; x = d.rot[x]; } while (x != d.head[v]);
        }
    }
    size_t cap = (size_t)g->m + (size_t)b;
    if (cap > INT_MAX/2) die("too many cut-open edges");
    h->u = mem(cap, sizeof(int)); h->v = mem(cap, sizeof(int));
    h->id = mem(cap, sizeof(int)); h->w = mem(cap, sizeof(Cost));
    for (int e = 0; e < g->m; ++e) if (!tree[e]) {
        int j = h->m++;
        h->u[j] = corner[2*e]; h->v[j] = corner[2*e+1];
        h->id[j] = e; h->w[j] = g->w[e];
    }
    for (int i = 0; i < b; ++i) {
        int e = out[i] / 2;
        if (e >= g->m) continue; /* forbidden terminal spokes */
        int j = h->m++;
        h->u[j] = i; h->v[j] = (i+1) % b; h->id[j] = e; h->w[j] = g->w[e];
    }
    h->dom_n = mem((size_t)k, sizeof(int));
    /* Every enumerated port carries a chord, so its image must touch a
       star edge: restrict candidate lists to usable boundary copies. */
    int *useful = mem((size_t)b, sizeof(int));
    int *size = mem((size_t)b, sizeof(int));
    for (int i = 0; i < b; ++i) ++size[h->center[i]];
    for (int i = 0; i < b; ++i) if (h->gap[i] >= 0) {
        useful[i] = size[h->center[i]] > 1;
        if (useful[i]) ++h->dom_n[h->gap[i]];
    }
    free(size);
    h->dom = mem((size_t)(b ? b : 1), sizeof(int));
    int *cursor = mem((size_t)k, sizeof(int));
    for (int i = 0; i < k; ++i) cursor[i] = i ? cursor[i-1] + h->dom_n[i-1] : 0;
    for (int i = 0; i < b; ++i) if (useful[i]) h->dom[cursor[h->gap[i]]++] = i;
    free(cursor); free(useful);
    h->all = mem((size_t)h->n, sizeof(int));
    for (int i = 0; i < h->n; ++i) h->all[i] = i;
    free(tree); free(rep); free(root); free(tag); free(out); free(corner);
    plane_free(&d);
}

void open_free(Open *h)
{
    free(h->label); free(h->gap); free(h->center);
    free(h->u); free(h->v); free(h->id); free(h->w); free(h->dist); free(h->first);
    free(h->dom); free(h->dom_n); free(h->all);
}
