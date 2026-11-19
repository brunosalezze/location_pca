/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file algorithms.h
 * @brief The solution methods of the paper, built on the Set Covering Framework.
 *
 *  - ALGORITHM_PCA: Progressive Covering Algorithm (Algorithm 4), the paper's main method.
 *  - ALGORITHM_EXACT: framework with the exact covering method (Algorithm 1 +
 *    Section 4.1.2); proves optimality when it finishes within the time limit.
 *    Equivalent to PCA with alpha0 = 1.
 *  - ALGORITHM_GREEDY: framework with the greedy covering heuristic
 *    (Algorithm 1 + Algorithm 2), used for the cleaning comparison of Table 2.
 */
#ifndef LOCATION_PCA_ALGORITHMS_H
#define LOCATION_PCA_ALGORITHMS_H

#include "cleaning.h"
#include "exact.h"
#include "graph.h"
#include "separation.h"

#include <stdint.h>

/** Solution method. */
typedef enum Algorithm {
    ALGORITHM_PCA,    /**< Progressive Covering Algorithm (Algorithm 4). */
    ALGORITHM_EXACT,  /**< Set Covering Framework with exact covering.   */
    ALGORITHM_GREEDY, /**< Set Covering Framework with greedy covering.  */
} Algorithm;

/** Parameters of one solution method. Fill with solver_params_defaults(). */
typedef struct SolverParams {
    Algorithm algorithm;            /**< Method to run. */
    double alpha0;                  /**< PCA: initial hybrid fraction alpha_0, in (0, 1]. */
    double time_limit_total;        /**< TL2: wall-clock limit for a whole run, seconds. */
    double time_limit_iteration;    /**< TL1: limit for one exact covering solve, seconds. */
    double time_limit_cleaning;     /**< Limit for one exact cleaning solve, seconds. */
    int32_t initial_paths_per_root; /**< Paths per O-D root in the first Cut phase (0 = all). */
    int32_t paths_per_root;         /**< Paths per O-D root in later Cut phases (0 = all). */
    SpreadCriterion spread;         /**< Spread criterion of the Cut phase. */
    CleaningMethod cleaning;        /**< Cleaning procedure applied after each covering. */
} SolverParams;

/** Outcome and statistics of one run. Release with run_result_free(). */
typedef struct RunResult {
    uint8_t *selected;       /**< [num_edges] 1 if the edge receives a counting station. */
    int32_t cost;            /**< Number of counting stations |C|. */
    int64_t connected_pairs; /**< O-D pairs C fails to separate (0 = feasible), checked afresh. */
    double lower_bound;      /**< Proven lower bound on the optimum, or NAN if none was computed. */
    int32_t num_paths;       /**< Final size of the path set P1. */
    int32_t iterations;      /**< Covering/Cut iterations performed. */
    int32_t exact_solves;    /**< Calls to the MIP solver (covering and cleaning). */
    double time_total;       /**< Wall-clock seconds of the whole run. */
    double time_solver;      /**< Seconds spent inside the MIP solver. */
    double time_separation;  /**< Seconds spent in the Cut phase. */
} RunResult;

/**
 * @brief Fills @p params with the defaults of @p algorithm.
 *
 * Defaults are the settings of the original experiments: TL2 = 7200 s,
 * TL1 = 500 s, cleaning TL = 25 s, MLC spread, exact cleaning, alpha0 = 1/16,
 * and paths per root of (all, 5) for PCA, (20, 10) for EXACT and (500, 500)
 * for GREEDY.
 *
 * @param params    Structure to fill.
 * @param algorithm Solution method.
 */
void solver_params_defaults(SolverParams *params, Algorithm algorithm);

/**
 * @brief Name of an algorithm as used on the command line ("pca", "exact", "greedy").
 *
 * @param algorithm Solution method.
 * @return Static string.
 */
const char *algorithm_name(Algorithm algorithm);

/**
 * @brief Solves the TCLP on @p graph with the method described by @p params.
 *
 * Inner workings (see algorithms.c for the step-by-step description):
 *  1. an initial Cut phase builds P1 from the empty cover (locked edges only);
 *  2. the chosen method repeatedly covers P1, cleans the cover and runs the
 *     Cut phase, until no uncovered O-D path exists, time runs out, or the
 *     lower bound reaches the cost of the best cover (which is then optimal);
 *  3. if no feasible cover was obtained within the time limit, the greedy
 *     framework completes one without time limit;
 *  4. the final cover is checked independently with
 *     cover_count_connected_od_pairs().
 *
 * @param graph  The network.
 * @param params Method and parameters.
 * @param solver MIP solver handle (needed by PCA, EXACT and exact cleaning).
 * @param seed   Seed of the run's random generator.
 * @param result Output; allocated by this function, release with run_result_free().
 */
void run_algorithm(const Graph *graph, const SolverParams *params, ExactSolver *solver,
                   uint64_t seed, RunResult *result);

/**
 * @brief Releases the memory owned by @p result.
 *
 * @param result Result filled by run_algorithm().
 */
void run_result_free(RunResult *result);

#endif /* LOCATION_PCA_ALGORITHMS_H */
