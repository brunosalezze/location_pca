/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file exact.h
 * @brief Exact covering method (Section 4.1.2): solves the unicost set covering
 *        model (1)-(3) of the paper with a MIP solver.
 *
 * Only the paths of P1 that are still uncovered become constraints, and only
 * the edges that lie on them become variables. Every edge already in the
 * cover is fixed. In the paper's notation:
 *
 *     min  Σ x_a
 *     s.t. Σ_{a ∈ p} x_a ≥ 1   for every uncovered path p ∈ P1
 *          x_a ∈ {0, 1}
 *
 * The interface is solver-independent. exact_gurobi.c implements it with
 * Gurobi. The unit tests link a brute-force implementation instead, so they
 * build without a MIP solver.
 */
#ifndef LOCATION_PCA_EXACT_H
#define LOCATION_PCA_EXACT_H

#include "cover.h"

#include <stdint.h>

/** Solver handle (license/environment). Create once, reuse for every solve. */
typedef struct ExactSolver ExactSolver;

/** Outcome of one call to exact_cover(). */
typedef struct ExactResult {
    double runtime; /**< Seconds spent by the MIP solver. */
    /**
     * Lower bound on the number of edges that must be added to the cover
     * to cover all of P1, or NAN if the solver could not provide one.
     */
    double objective_bound;
    bool has_solution;    /**< A solution was found and its edges were added to the cover. */
    int32_t solver_calls; /**< MIP models actually solved: 0 if no path was uncovered, else 1. */
} ExactResult;

/**
 * @brief Starts the MIP solver (checks the license and sets global parameters).
 *
 * @param threads         Number of solver threads (0 = solver decides).
 * @param show_solver_log Forward the solver's log to stderr.
 * @return Newly allocated handle; release it with exact_solver_free(). Stops
 *         the program with an explanatory message if the solver cannot start.
 */
[[nodiscard]] ExactSolver *exact_solver_create(int threads, bool show_solver_log);

/**
 * @brief Releases the solver handle.
 *
 * @param solver Handle to release, or nullptr.
 */
void exact_solver_free(ExactSolver *solver);

/**
 * @brief Human-readable solver name and version, e.g. "Gurobi 13.0.3".
 *
 * @return Pointer to a static string.
 */
const char *exact_solver_description(void);

/**
 * @brief Completes @p cover so that it covers every path of P1, using as few edges as possible.
 *
 * Inner workings:
 *  1. one binary variable is created per unselected edge that lies on an
 *     uncovered path (cover->uncovered_paths_through[e] > 0);
 *  2. one covering constraint is created per uncovered path;
 *  3. the model is solved within @p time_limit seconds, starting from
 *     @p warm_start if given. If @p known_bound is positive, the solver stops
 *     as soon as it finds a solution with that many edges: such a solution is
 *     optimal, and proving the bound a second time is usually the expensive
 *     part of the solve;
 *  4. if a solution is available, the edges with x_a = 1 are added to @p cover.
 *
 * If no solution is found in time, @p cover is left unchanged and
 * has_solution is false. The caller then completes it heuristically.
 * If no path is uncovered there is nothing to solve: the solver is not called
 * and the result is "0 edges added, bound 0".
 *
 * @param solver     Solver handle.
 * @param cover      Cover to complete (synchronized with its pool).
 * @param time_limit Time limit in seconds (> 0).
 * @param warm_start Optional [num_edges] flags of a known feasible
 *                   completion (MIP start), or nullptr.
 * @param known_bound Number of edges that every solution of this model is
 *                   already known to need (a valid lower bound on its
 *                   optimum, an integer), or 0 if nothing is known.
 * @return Runtime, bound, whether a solution was found and whether the solver
 *         ran. The returned bound is the solver's own and may be smaller than
 *         @p known_bound when the solve stopped early.
 */
ExactResult exact_cover(ExactSolver *solver, Cover *cover, double time_limit,
                        const uint8_t *warm_start, double known_bound);

#endif /* LOCATION_PCA_EXACT_H */
