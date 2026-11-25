/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
/* Development harness: builds small graphs and checks embed(). */
#include "../cut.h"

int g_crossings = 32;
static int failures;

static void check(const char *name, int n, int m, const int (*edges)[2], int expect)
{
    Plane g = plane(n, m);
    for (int e = 0; e < m; ++e) {
        g.at[2*e] = edges[e][0];
        g.at[2*e+1] = edges[e][1];
        g.w[e] = 1;
    }
    Planarity result = embed(&g);
    bool ok = result == PLANAR;
    int nf = 0;
    if (ok) {
        free(faces(&g, &nf));
        int seen = 0;
        for (int v = 0; v < n; ++v) if (g.head[v] >= 0) {
            int d = g.head[v], a = d;
            do {
                if (a < 0 || a >= 2*m || g.at[a] != v || ++seen > 2*m)
                    die("invalid dart in test embedding");
                a = g.rot[a];
            } while (a != d);
        }
        if (seen != 2*m) die("missing darts in test embedding");
    }
    bool good = expect ? ok && nf == 2 - n + m : result == NONPLANAR;
    if (!good) ++failures;
    printf("%-14s %-9s %s\n", name, ok ? "planar" : "rejected",
           good ? "ok" : "WRONG");
    plane_free(&g);
}

int main(void)
{
    const int tri[][2] = {{0,1},{1,2},{2,0}};
    check("triangle", 3, 3, tri, 1);
    const int sqd[][2] = {{0,1},{1,2},{2,3},{3,0},{0,2}};
    check("square-diag", 4, 5, sqd, 1);
    const int k4[][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    check("K4", 4, 6, k4, 1);
    const int k5[][2] = {{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}};
    check("K5", 5, 10, k5, 0);
    const int k33[][2] = {{0,3},{0,4},{0,5},{1,3},{1,4},{1,5},{2,3},{2,4},{2,5}};
    check("K33", 6, 9, k33, 0);
    const int grid[][2] = {{0,1},{1,2},{3,4},{4,5},{6,7},{7,8},{0,3},{1,4},{2,5},{3,6},{4,7},{5,8}};
    check("grid3x3", 9, 12, grid, 1);
    const int cube[][2] = {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
    check("cube", 8, 12, cube, 1);
    const int bowtie[][2] = {{0,1},{1,2},{2,0},{2,3},{3,4},{4,2}};
    check("two-triangles", 5, 6, bowtie, 1);
    const int bridge[][2] = {{0,1},{1,2},{2,0},{2,3}};
    check("bridge", 4, 4, bridge, 1);
    const int parallel[][2] = {{0,1},{0,1},{1,2}};
    check("parallel", 3, 3, parallel, 1);
    const int loop[][2] = {{0,1},{1,1},{1,2},{2,0}};
    check("self-loop", 3, 4, loop, 1);
    const int petersen[][2] = {{0,1},{1,2},{2,3},{3,4},{4,0},{5,7},{7,9},{9,6},{6,8},{8,5},
                               {0,5},{1,6},{2,7},{3,8},{4,9}};
    check("petersen", 10, 15, petersen, 0);
    /* K7: m = 21 > 3n-6 = 15, non-planar */
    int k7[21][2], nk = 0;
    for (int a = 0; a < 7; ++a) for (int b = a+1; b < 7; ++b)
        { k7[nk][0] = a; k7[nk][1] = b; ++nk; }
    check("K7", 7, nk, (const int (*)[2])k7, 0);
    /* 16-node planar grid with a diagonal in each cell */
    int big[39][2], nbig = 0;
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c) {
            if (c < 3) { big[nbig][0] = r*4+c; big[nbig][1] = r*4+c+1; ++nbig; }
            if (r < 3) { big[nbig][0] = r*4+c; big[nbig][1] = (r+1)*4+c; ++nbig; }
            if (r < 3 && c < 3) { big[nbig][0] = r*4+c; big[nbig][1] = (r+1)*4+c+1; ++nbig; }
        }
    check("grid4x4-diag", 16, nbig, (const int (*)[2])big, 1);
    /* octahedron = K_{2,2,2}: 6 vertices, 12 edges, maximal planar */
    const int octa[][2] = {{0,2},{0,4},{0,3},{0,5},{1,2},{1,4},{1,3},{1,5},
                           {2,4},{4,3},{3,5},{5,2}};
    check("octahedron", 6, 12, octa, 1);
    const int two_parallel[][2] = {{0,1},{1,0}};
    check("two-parallel", 2, 2, two_parallel, 1);
    /* Two edge-bearing components plus an isolated vertex. */
    Plane split = plane(4, 2);
    split.at[0] = 0; split.at[1] = 1;
    split.at[2] = split.at[3] = 2;
    if (embed(&split) != PLANAR || split.head[3] != -1) ++failures;
    plane_free(&split);
    printf("%s\n", failures ? "embed: FAILED" : "embed: all ok");
    return failures != 0;
}
