/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#define _POSIX_C_SOURCE 200809L
#include "cut.h"
#include <ctype.h>
#include <limits.h>

typedef struct { char *key; int id; bool terminal; } Slot;
typedef struct { Slot *slots; size_t cap; int n; } Map;
typedef struct { char *id; int a, b; } Edge;

static size_t hash(const char *p)
{
    size_t h = 2166136261u;
    for (; *p; ++p) h = (h ^ (unsigned char)*p)*16777619u;
    return h;
}

static Slot *lookup(Slot *slots, size_t cap, const char *key)
{
    size_t i = hash(key) & (cap-1);
    while (slots[i].key && strcmp(slots[i].key, key)) i = (i+1) & (cap-1);
    return &slots[i];
}

static char *copy(const char *s)
{
    size_t n = strlen(s)+1;
    return memcpy(mem(n, 1), s, n);
}

static Slot *intern(Map *m, const char *key)
{
    if ((size_t)m->n >= m->cap/2) {
        size_t cap = m->cap ? 2*m->cap : 256;
        Slot *slots = mem(cap, sizeof(Slot));
        for (size_t i = 0; i < m->cap; ++i)
            if (m->slots[i].key) *lookup(slots, cap, m->slots[i].key) = m->slots[i];
        free(m->slots);
        m->slots = slots; m->cap = cap;
    }
    Slot *s = lookup(m->slots, m->cap, key);
    if (!s->key) {
        if (m->n == INT_MAX) die("too many identifiers");
        *s = (Slot){.key = copy(key), .id = m->n++};
    }
    return s;
}

static void map_free(Map *m)
{
    for (size_t i = 0; i < m->cap; ++i) free(m->slots[i].key);
    free(m->slots);
}

static FILE *open_input(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(EXIT_FAILURE); }
    return f;
}

static char *token(char **p)
{
    while (isspace((unsigned char)**p)) ++*p;
    if (!**p || **p == '#') return NULL;
    char *start = *p;
    while (**p && !isspace((unsigned char)**p)) ++*p;
    if (**p) *(*p)++ = '\0';
    return start;
}

static void *grow(void *p, size_t *cap, size_t width)
{
    size_t n = *cap ? 2*(*cap) : 256;
    if (n > SIZE_MAX/width) die("instance too large");
    p = realloc(p, n*width);
    if (!p) die("out of memory");
    *cap = n;
    return p;
}

Inst tclp_read(const char *dir, const char *edges_path, const char *od_path)
{
    char epath[4096], opath[4096];
    if (dir) {
        int e = snprintf(epath, sizeof epath, "%s/edges.txt", dir);
        int o = snprintf(opath, sizeof opath, "%s/od_nodes.txt", dir);
        if (e < 0 || (size_t)e >= sizeof epath || o < 0 || (size_t)o >= sizeof opath)
            die("instance path too long");
        edges_path = epath; od_path = opath;
    }
    Map nodes = {0}, ids = {0};
    Edge *edges = NULL;
    size_t edge_cap = 0, line_cap = 0, term_cap = 0;
    int m = 0, k = 0, *term = NULL;
    char *line = NULL;
    FILE *f = open_input(edges_path);
    while (getline(&line, &line_cap, f) >= 0) {
        char *p = line, *id = token(&p);
        if (!id) continue;
        char *a = token(&p), *b = token(&p), *lanes = token(&p);
        if (!a || !b || !lanes || token(&p))
            die("edges.txt: expected edge_id node_a node_b lanes");
        if (m == INT_MAX/2) die("too many edges");
        if (intern(&ids, id)->id != m) die("edges.txt: duplicate edge ID");
        if ((size_t)m == edge_cap) edges = grow(edges, &edge_cap, sizeof *edges);
        int u = intern(&nodes, a)->id, v = intern(&nodes, b)->id;
        edges[m++] = (Edge){copy(id), u, v};
    }
    if (!feof(f)) die("edges.txt: read failed");
    fclose(f);
    map_free(&ids);

    f = open_input(od_path);
    while (getline(&line, &line_cap, f) >= 0) {
        char *p = line, *id = token(&p);
        if (!id) continue;
        if (token(&p)) die("od_nodes.txt: expected one node ID per line");
        Slot *s = intern(&nodes, id);
        if (s->terminal) continue;
        s->terminal = true;
        if ((size_t)k == term_cap) term = grow(term, &term_cap, sizeof *term);
        term[k++] = s->id;
    }
    if (!feof(f)) die("od_nodes.txt: read failed");
    fclose(f); free(line);

    Inst ins = {.g = plane(nodes.n, m), .k = k, .term = term,
                .edge_id = mem((size_t)m, sizeof(char *))};
    for (int e = 0; e < m; ++e) {
        ins.g.at[2*e] = edges[e].a;
        ins.g.at[2*e+1] = edges[e].b;
        ins.g.w[e] = 1;
        ins.edge_id[e] = edges[e].id;
    }
    map_free(&nodes); free(edges);
    return ins;
}

void tclp_free(Inst *ins)
{
    for (int e = 0; e < ins->g.m; ++e) free(ins->edge_id[e]);
    free(ins->edge_id); free(ins->term);
    plane_free(&ins->g);
}
