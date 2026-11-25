/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "cut.h"

/* Embed each biconnected block by adding the paths of an open ear
   decomposition. Try every face that admits the next ear, backtracking
   if needed. A search limit is inconclusive, not evidence of nonplanarity. */
enum { SEARCH_BUDGET = 20000000 };

typedef struct { int *a, n, cap; } List;

static void push(List *l, int x)
{
    if (l->n == l->cap) {
        l->cap = l->cap ? 2*l->cap : 8;
        l->a = realloc(l->a, (size_t)l->cap*sizeof *l->a);
        if (!l->a) die("out of memory");
    }
    l->a[l->n++] = x;
}

typedef struct {
    Plane *g;
    int *head, *next, *disc, *low, *prev, *block_head;
    int clock, budget;
    List stack;
    Planarity result;
} Search;

/* With the predecessor face rule, a face exiting v on dart d has its
   corner immediately after d in the vertex rotation. */
static void insert(Search *s, int d, int a)
{
    int b = s->g->rot[d];
    s->g->rot[d] = a; s->prev[a] = d;
    s->g->rot[a] = b; s->prev[b] = a;
}

static void remove_dart(Search *s, int a)
{
    int d = s->prev[a], b = s->g->rot[a];
    s->g->rot[d] = b; s->prev[b] = d;
}

static void add_at_vertex(Search *s, int a)
{
    int *h = &s->block_head[s->g->at[a]];
    if (*h < 0) {
        *h = a;
        s->g->rot[a] = s->prev[a] = a;
    } else insert(s, *h, a);
}

static Planarity place(Search *s, List *ears, int i, int n)
{
    if (i == n) return PLANAR;
    Plane *g = s->g;
    List ear = ears[i];
    int a = ear.a[0], b = ear.a[ear.n-1] ^ 1;
    int h = s->block_head[g->at[a]], v = g->at[b], d = h;
    do {
        int f = d;
        do {
            if (--s->budget == 0) return PLANARITY_UNKNOWN;
            if (g->at[f] == v) {
                insert(s, d, a); insert(s, f, b);
                for (int j = 1; j < ear.n; ++j) {
                    int p = ear.a[j-1] ^ 1, q = ear.a[j];
                    g->rot[p] = s->prev[p] = q;
                    g->rot[q] = s->prev[q] = p;
                    s->block_head[g->at[q]] = q;
                }
                Planarity r = place(s, ears, i+1, n);
                if (r != NONPLANAR) return r;
                remove_dart(s, a); remove_dart(s, b);
                for (int j = 1; j < ear.n; ++j) s->block_head[g->at[ear.a[j]]] = -1;
            }
            f = s->prev[f ^ 1];
        } while (f != d);
        d = g->rot[d];
    } while (d != h);
    return NONPLANAR;
}

/* Compute ears in a fixed order. After the initial cycle, run a
   multi-source BFS from its vertices (and those of previous ears). An
   unused edge joining different BFS trees supplies the next ear. */
static int decompose(const Plane *h, const int *hd, const int *nx, List *ears)
{
    int *src = mem((size_t)h->n, sizeof(int)), *par = mem((size_t)h->n, sizeof(int));
    int *q = mem((size_t)h->n, sizeof(int));
    bool *in = mem((size_t)h->n, sizeof(bool)), *used = mem((size_t)h->m, sizeof(bool));
    for (int v = 0; v < h->n; ++v) src[v] = -1;
    int a = h->at[0], b = h->at[1], front = 0, back = 0;
    q[back++] = a; src[a] = a;
    while (front < back && src[b] < 0) {
        int u = q[front++];
        for (int d = hd[u]; d >= 0; d = nx[d]) {
            int v = h->at[d ^ 1];
            if (d/2 == 0 || src[v] >= 0) continue;
            src[v] = a; par[v] = d; q[back++] = v;
        }
    }
    if (src[b] < 0) die("internal error: no initial cycle in block");
    push(&ears[0], 0);
    for (int v = b; v != a; v = h->at[par[v]]) push(&ears[0], par[v] ^ 1);

    int n = 0, left = h->m;
    while (true) {
        for (int i = 0; i < ears[n].n; ++i) {
            int d = ears[n].a[i];
            used[d/2] = true;
            in[h->at[d]] = in[h->at[d ^ 1]] = true;
            --left;
        }
        ++n;
        if (!left) break;
        front = back = 0;
        for (int v = 0; v < h->n; ++v) {
            src[v] = in[v] ? v : -1;
            if (in[v]) q[back++] = v;
        }
        int join = -1;
        while (front < back && join < 0) {
            int u = q[front++];
            for (int d = hd[u]; d >= 0; d = nx[d]) {
                if (used[d/2]) continue;
                int v = h->at[d ^ 1];
                if (src[v] < 0) {
                    src[v] = src[u]; par[v] = d; q[back++] = v;
                } else if (src[v] != src[u]) {
                    join = d;
                    break;
                }
            }
        }
        if (join < 0) die("internal error: no open ear in block");
        List *ear = &ears[n];
        for (int u = h->at[join]; !in[u]; u = h->at[par[u]]) push(ear, par[u]);
        for (int l = 0, r = ear->n-1; l < r; ++l, --r) {
            int d = ear->a[l]; ear->a[l] = ear->a[r]; ear->a[r] = d;
        }
        push(ear, join);
        for (int v = h->at[join ^ 1]; !in[v]; v = h->at[par[v]]) push(ear, par[v] ^ 1);
    }
    free(src); free(par); free(q); free(in); free(used);
    return n;
}

static void merge_rotation(Search *s, int v, int first)
{
    int old = s->g->head[v];
    if (old < 0) s->g->head[v] = first;
    else {
        int last = s->prev[first], after = s->g->rot[old];
        s->g->rot[old] = first; s->prev[first] = old;
        s->g->rot[last] = after; s->prev[after] = last;
    }
}

/* The stack interval contains only the block edges. Use local vertex
   IDs for ear decomposition and translate the resulting darts back. */
static Planarity block(Search *s, int start, int *local)
{
    int m = s->stack.n-start;
    const int *edges = s->stack.a+start;
    if (m == 1) {
        int d = edges[0];
        for (int i = 0; i < 2; ++i, d ^= 1) {
            s->g->rot[d] = s->prev[d] = d;
            merge_rotation(s, s->g->at[d], d);
        }
        return PLANAR;
    }
    List verts = {0};
    int *at = mem(2*(size_t)m, sizeof(int));
    for (int d = 0; d < 2*m; ++d) {
        int u = s->g->at[edges[d/2] ^ (d & 1)];
        if (local[u] < 0) { local[u] = verts.n; push(&verts, u); }
        at[d] = local[u];
    }
    int *hd = mem((size_t)verts.n, sizeof(int)), *nx = mem(2*(size_t)m, sizeof(int));
    for (int v = 0; v < verts.n; ++v) hd[v] = -1;
    for (int d = 0; d < 2*m; ++d) { nx[d] = hd[at[d]]; hd[at[d]] = d; }
    Plane h = {.n = verts.n, .m = m, .at = at};
    List *ears = mem((size_t)m, sizeof(List));
    int n = decompose(&h, hd, nx, ears);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < ears[i].n; ++j) {
            int d = ears[i].a[j];
            ears[i].a[j] = edges[d/2] ^ (d & 1);
        }
    for (int j = 0; j < ears[0].n; ++j) {
        int d = ears[0].a[j];
        add_at_vertex(s, d); add_at_vertex(s, d ^ 1);
    }
    Planarity r = place(s, ears, 1, n);
    for (int v = 0; v < verts.n; ++v) {
        int u = verts.a[v];
        if (r == PLANAR) merge_rotation(s, u, s->block_head[u]);
        s->block_head[u] = local[u] = -1;
    }
    for (int i = 0; i < n; ++i) free(ears[i].a);
    free(ears); free(verts.a); free(hd); free(nx); free(at);
    return r;
}

static void dfs(Search *s, int u, int incoming, int *local)
{
    s->low[u] = s->disc[u] = s->clock++;
    for (int d = s->head[u]; d >= 0 && s->result == PLANAR; d = s->next[d]) {
        if (d == (incoming ^ 1)) continue;
        int v = s->g->at[d ^ 1];
        if (s->disc[v] < 0) {
            int start = s->stack.n;
            push(&s->stack, d);
            dfs(s, v, d, local);
            if (s->result != PLANAR) return;
            if (s->low[v] < s->low[u]) s->low[u] = s->low[v];
            if (s->low[v] >= s->disc[u]) {
                s->result = block(s, start, local);
                s->stack.n = start;
            }
        } else if (s->disc[v] < s->disc[u]) {
            push(&s->stack, d);
            if (s->disc[v] < s->low[u]) s->low[u] = s->disc[v];
        }
    }
}

Planarity embed(Plane *g)
{
    size_t n = (size_t)g->n, m = 2*(size_t)g->m;
    Search s = {.g = g, .head = mem(n, sizeof(int)), .next = mem(m, sizeof(int)),
                .disc = mem(n, sizeof(int)), .low = mem(n, sizeof(int)),
                .prev = mem(m, sizeof(int)), .block_head = mem(n, sizeof(int)),
                .result = PLANAR, .budget = SEARCH_BUDGET};
    int *local = mem(n, sizeof(int));
    for (int v = 0; v < g->n; ++v)
        s.head[v] = s.disc[v] = s.block_head[v] = local[v] = g->head[v] = -1;
    for (int d = 0; d < 2*g->m; ++d) {
        int u = g->at[d];
        g->rot[d] = -1;
        if (u == g->at[d ^ 1]) continue;
        s.next[d] = s.head[u]; s.head[u] = d;
    }
    int components = 0;
    for (int v = 0; v < g->n && s.result == PLANAR; ++v)
        if (s.disc[v] < 0 && s.head[v] >= 0) {
            dfs(&s, v, -1, local);
            ++components;
        }
    if (s.result == PLANAR) {
        for (int d = 0; d < 2*g->m; d += 2) if (g->at[d] == g->at[d+1]) {
            int v = g->at[d];
            if (g->head[v] < 0) ++components;
            g->rot[d] = s.prev[d] = d+1; g->rot[d+1] = s.prev[d+1] = d;
            merge_rotation(&s, v, d);
        }
        int vertices = 0, nf;
        for (int v = 0; v < g->n; ++v) vertices += g->head[v] >= 0;
        free(faces(g, &nf));
        if (g->m && (int64_t)vertices - g->m + nf != 2*(int64_t)components)
            die("internal error: constructed embedding failed Euler validation");
    }
    free(s.head); free(s.next); free(s.disc); free(s.low); free(s.prev);
    free(s.block_head); free(s.stack.a); free(local);
    return s.result;
}
