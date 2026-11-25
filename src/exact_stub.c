/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
/**
 * @file exact_stub.c
 * @brief Stand-in for exact_gurobi.c used when the program is built without
 *        Gurobi (make NO_GUROBI=1).
 *
 * The planar algorithm (-a planar) never calls the MIP solver, so a build
 * without Gurobi is fully useful: it computes certified minimum cuts on
 * planar networks. The methods built on the Set Covering Framework stop
 * with this explanation when selected.
 */
#include "exact.h"

#include "util.h"

#include <stdlib.h>

struct ExactSolver {
    int unused;
};

ExactSolver *exact_solver_create(int threads, bool show_solver_log)
{
    (void)threads;
    (void)show_solver_log;
    fatal("this build has no MIP solver; only '-a planar' and '-a greedy --cleaning none' work");
}

void exact_solver_free(ExactSolver *solver)
{
    free(solver);
}

const char *exact_solver_description(void)
{
    return "not available (built without Gurobi)";
}

ExactResult exact_cover(ExactSolver *solver, Cover *cover, double time_limit,
                        const uint8_t *warm_start, double known_bound)
{
    (void)solver;
    (void)cover;
    (void)time_limit;
    (void)warm_start;
    (void)known_bound;
    fatal("this build has no MIP solver");
}
