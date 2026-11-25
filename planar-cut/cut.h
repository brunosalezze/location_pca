/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CUT_H
#define CUT_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint64_t Cost;
typedef unsigned Mask;
#define INF UINT64_MAX
#define MAX_TERMINALS 20

/* Edge e has outgoing darts 2e and 2e+1, pointing in opposite directions. */
typedef struct {
    int n, m, *at, *rot, *head;
    Cost *w;
} Plane;

typedef struct {
    int m;
    Cost cost;
    bool *cut;
} Answer;

/* G1: boundary copies are numbered first, in contour order from rep 0.
   Each interior star is represented by its center for each boundary copy. */
typedef struct {
    int n, b, m, k, input_m, *label, *gap, *center;
    int *u, *v, *id;
    Cost *w, *dist;
    int *first;
    /* Per-gap candidate lists for boundary positions, and the full
       vertex list used by unconstrained branch positions. */
    int *dom, *dom_n, *all;
} Open;

typedef enum { CHORD, PATH, ORDER } Kind;
typedef struct { int u, v; Kind kind; } Link;
typedef struct {
    int n, b, l;
    int *gap;
    int *eg; /* exterior group per position, set during enumeration */
    Link *link;
} Topology;

void die(const char *msg);
void *mem(size_t n, size_t s);
Cost plus(Cost a, Cost b);
Cost mincut(const Plane *g, int s, int t, bool *cut);
void mincut_max_side(const Plane *g, int s, int t, bool *side);
Plane plane(int n, int m);
void plane_free(Plane *g);
int *faces(const Plane *g, int *nf);
void open_build(const Plane *g, const int *term, int k, Open *h);
void open_free(Open *h);
void shortest(Open *h);
void enumerate(const Open *h, Answer *best, const Cost *other, Cost lb);
Cost realize(const Open *h, const Topology *t, Cost limit, int *map);
void profile_init(void);
Answer solve(const Plane *g, const int *term, int k);
bool solve_narrow(const Plane *g, const int *term, int k, Answer *out);

/* Finds and validates a planar rotation system, within a search budget. */
typedef enum { NONPLANAR, PLANAR, PLANARITY_UNKNOWN } Planarity;
Planarity embed(Plane *g);

/* Instance in the traffic counting location format of location_pca. */
typedef struct {
    Plane g;
    int k;
    int *term;
    char **edge_id;         /* the m edge identifiers, in input order */
} Inst;

Inst tclp_read(const char *dir, const char *edges_path, const char *od_path);
void tclp_free(Inst *ins);

/* Cap on interior crossings per terminal (paper: |cr| <= 24k). */
extern int g_crossings;

#endif
