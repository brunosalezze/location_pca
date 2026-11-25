/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "cut.h"
#include <limits.h>

typedef struct {
    int n, *head, *to, *next;
    Cost *w;
} Adj;

static Adj adj(int n, int m)
{
    Adj a = {n, mem((size_t)n, sizeof(int)), mem(2*(size_t)m, sizeof(int)),
             mem(2*(size_t)m, sizeof(int)), mem((size_t)m, sizeof(Cost))};
    for (int v = 0; v < n; ++v) a.head[v] = -1;
    return a;
}

static void arc(Adj *a, int d, int u, int v)
{
    a->to[d] = v; a->next[d] = a->head[u]; a->head[u] = d;
}

static void adj_free(Adj *a)
{
    free(a->head); free(a->to); free(a->next); free(a->w);
}

typedef struct {
    int count, *heap, *pos;
    Cost *dist;
} Queue;

static void move(Queue *q, int i, int v)
{
    q->heap[i] = v; q->pos[v] = i;
}

static void raise(Queue *q, int v)
{
    int i = q->pos[v];
    if (i < 0) i = q->count++;
    while (i > 0) {
        int p = (i-1)/2, u = q->heap[p];
        if (q->dist[u] <= q->dist[v]) break;
        move(q, i, u); i = p;
    }
    move(q, i, v);
}

static int pop(Queue *q)
{
    int result = q->heap[0], v = q->heap[--q->count], i = 0;
    if (q->count) {
        while (i < q->count/2) {
            int c = 2*i+1;
            if (c+1 < q->count && q->dist[q->heap[c+1]] < q->dist[q->heap[c]]) ++c;
            if (q->dist[v] <= q->dist[q->heap[c]]) break;
            move(q, i, q->heap[c]); i = c;
        }
        move(q, i, v);
    }
    q->pos[result] = -1;
    return result;
}

static void close_dist(const Adj *a, Cost *dist, int *parent, Queue *q)
{
    q->count = 0; q->dist = dist;
    for (int v = 0; v < a->n; ++v) q->pos[v] = -1;
    for (int v = 0; v < a->n; ++v) if (dist[v] < INF) raise(q, v);
    while (q->count) {
        int u = pop(q);
        for (int d = a->head[u]; d >= 0; d = a->next[d]) {
            int v = a->to[d];
            Cost c = plus(dist[u], a->w[d/2]);
            if (c < dist[v]) {
                dist[v] = c;
                parent[v] = d;
                raise(q, v);
            }
        }
    }
}

/* Dreyfus–Wagner on the old dual. Each rep is a leaf: a Steiner
   tree with an extra M-spoke loses to one spoke per rep plus the
   entire old dual. So drop the common k*M and use multi-source seeds. */
void steiner(const Plane *dual, int old_m, const int *rep, int k, bool *tree)
{
    int n = dual->n;
    size_t stride = (size_t)n;
    Mask end = 1u << k;
    if (stride > SIZE_MAX/end) die("Steiner table too large");
    size_t cells = stride*end;
    Cost *value = mem(cells, sizeof(Cost));
    int *parent = mem(cells, sizeof(int));
    Queue q = {0, mem(stride, sizeof(int)), mem(stride, sizeof(int)), NULL};
    Adj a = adj(n, old_m);
    for (int d = 0; d < 2*old_m; ++d)
        arc(&a, d, dual->at[d], dual->at[d ^ 1]);
    memcpy(a.w, dual->w, (size_t)old_m * sizeof(Cost));
    for (Mask set = 1; set < end; ++set) {
        Cost *row = value + set*stride;
        int *p = parent + set*stride;
        for (int v = 0; v < n; ++v) { row[v] = INF; p[v] = -1; }
        if (!(set & (set-1))) {
            int i = 0;
            while ((1u << i) != set) ++i;
            int d = dual->head[rep[i]];
            do {
                int v = dual->at[d ^ 1];
                row[v] = 0; p[v] = d; /* high-cost spoke marks a leaf */
                d = dual->rot[d];
            } while (d != dual->head[rep[i]]);
        } else {
            for (Mask s = (set-1)&set; s; s = (s-1)&set) {
                Mask t = set ^ s;
                if (s >= t) continue;
                for (int v = 0; v < n; ++v) {
                    size_t j = (size_t)v;
                    Cost c = plus(value[s*stride+j], value[t*stride+j]);
                    if (c < row[v]) {
                        row[v] = c; p[v] = -(int)s-1;
                    }
                }
            }
        }
        close_dist(&a, row, p, &q);
    }
    Mask full = end-1;
    int v = 0;
    for (int u = 1; u < n; ++u)
        if (value[full*stride+(size_t)u] < value[full*stride+(size_t)v]) v = u;
    if (value[full*stride+(size_t)v] == INF) die("no dual Steiner tree");
    Mask *sets = mem((size_t)k, sizeof(Mask));
    int *vertices = mem((size_t)k, sizeof(int));
    int top = 0;
    for (;;) {
        int d = parent[full*stride+(size_t)v];
        if (d >= 0) {
            tree[d/2] = true;
            if (d/2 < old_m) { v = a.to[d ^ 1]; continue; }
            if (!top) break;
            --top; full = sets[top]; v = vertices[top];
        } else {
            Mask s = (Mask)(-1-d);
            if (!s || top == k) die("invalid Steiner traceback");
            sets[top] = full ^ s; vertices[top++] = v; full = s;
        }
    }
    free(sets); free(vertices); free(parent); free(value);
    free(q.heap); free(q.pos); adj_free(&a);
}

void shortest(Open *h)
{
    size_t n = (size_t)h->n;
    if (h->m > INT_MAX/2 || n > SIZE_MAX/n) die("cut-open graph too large");
    h->dist = mem(n*n, sizeof(Cost));
    h->first = mem(n*n, sizeof(int)); /* predecessor dart, indexed by source then dest */
    Queue q = {0, mem(n, sizeof(int)), mem(n, sizeof(int)), NULL};
    Adj a = adj(h->n, h->m);
    for (int e = 0; e < h->m; ++e) {
        arc(&a, 2*e, h->u[e], h->v[e]); arc(&a, 2*e+1, h->v[e], h->u[e]);
        a.w[e] = h->w[e];
    }
    for (size_t u = 0; u < n; ++u) {
        Cost *row = h->dist + u*n;
        int *p = h->first + u*n;
        for (size_t v = 0; v < n; ++v) { row[v] = INF; p[v] = -1; }
        row[u] = 0;
        close_dist(&a, row, p, &q);
    }
    free(q.heap); free(q.pos); adj_free(&a);
}
