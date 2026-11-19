/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file exact_bruteforce.c
 * @brief Stand-in implementation of exact.h for the unit tests, used instead
 *        of exact_gurobi.c so that the test binary builds without a MIP solver.
 *
 * exact_cover() solves the unicost set covering model on the paths of P1 that
 * are still uncovered by enumerating the subsets of candidate edges in
 * increasing order of size, so the first feasible subset found is an optimal
 * one. Candidates are the unselected edges that lie on an uncovered path
 * (cover->uncovered_paths_through[e] > 0). The search is exponential, which is
 * fine for the small instances of the unit tests; a few bounds make it usable
 * on slightly larger ones:
 *  - a positive @p known_bound starts the enumeration at that size (a
 *    solution of exactly that size is optimal by its definition);
 *  - a feasible warm start caps the enumerated size and is applied if no
 *    smaller solution exists;
 *  - when the model is too large for this stand-in (more candidates than
 *    BF_MAX_CANDIDATES, or more uncovered paths than bits in a mask), the warm
 *    start is applied when feasible and the greedy heuristic completes the
 *    cover otherwise. The objective bound is then NAN: no bound was proved.
 *
 * The time limit is ignored (the enumeration is not interruptible); tests keep
 * the models small enough that this does not matter.
 */
#include "exact.h"

#include "cover.h"
#include "greedy.h"
#include "path_pool.h"
#include "rng.h"
#include "timer.h"
#include "util.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

/** Largest number of candidate edges the enumeration accepts. */
#define BF_MAX_CANDIDATES 22

/** Largest number of uncovered paths representable with one uint64_t mask. */
#define BF_MAX_PATHS 64

/** Solver handle of the stand-in; keeps the creation parameters for nothing. */
struct ExactSolver {
    int threads;
    bool show_solver_log;
};

ExactSolver *exact_solver_create(int threads, bool show_solver_log)
{
    ExactSolver *solver = xcalloc(1, sizeof *solver);
    solver->threads = threads;
    solver->show_solver_log = show_solver_log;
    return solver;
}

void exact_solver_free(ExactSolver *solver)
{
    free(solver);
}

const char *exact_solver_description(void)
{
    return "brute-force set covering (unit-test stand-in for a MIP solver)";
}

/* ------------------------------------------------------------------------- */
/* Brute-force model                                                         */
/* ------------------------------------------------------------------------- */

/** The covering model of one exact_cover() call, with paths as bitmask bits. */
typedef struct BruteForceModel {
    int32_t num_candidates;
    int32_t candidates[BF_MAX_CANDIDATES]; /**< Edge indices. */
    uint64_t cover_mask[BF_MAX_CANDIDATES]; /**< Bit set of uncovered paths that choosing candidate c covers. */
    uint64_t suffix_or[BF_MAX_CANDIDATES + 1]; /**< cover_mask[c] | ... | cover_mask[n-1]. */
    uint64_t required;                       /**< Bit set of all uncovered paths. */
    int32_t chosen[BF_MAX_CANDIDATES];       /**< Candidate indices of the subset being built. */
    int32_t solution[BF_MAX_CANDIDATES];     /**< Candidate indices of the best subset found. */
} BruteForceModel;

/**
 * @brief Collects the uncovered paths of @p cover as bits, and the candidate
 *        edges with, for each candidate, the set of uncovered paths it covers.
 *
 * @param cover      Cover synchronized with its pool.
 * @param model      Model to fill.
 * @param capacity   Maximum number of candidates accepted.
 * @return 0 on success; 1 if the model is too large to enumerate (the caller
 *         then falls back to the warm start or the greedy heuristic).
 */
static int build_model(const Cover *cover, BruteForceModel *model, int32_t capacity)
{
    const Graph *graph = cover->graph;
    const PathPool *pool = cover->pool;
    *model = (BruteForceModel){0};

    /* Candidate edges: unselected, on at least one uncovered path. */
    for (int32_t e = 0; e < graph->num_edges; e++) {
        if (!cover->selected[e] && cover->uncovered_paths_through[e] > 0) {
            if (model->num_candidates == capacity) {
                return 1;
            }
            model->candidates[model->num_candidates++] = e;
        }
    }

    /* One bit per uncovered path; each candidate accumulates the bits of the
     * uncovered paths it lies on. */
    int32_t bit = 0;
    for (int32_t p = 0; p < cover->num_paths; p++) {
        if (cover->path_hits[p] > 0) {
            continue;
        }
        if (bit == BF_MAX_PATHS) {
            return 1;
        }
        const int32_t *edges = path_pool_edges(pool, p);
        int32_t length = path_pool_length(pool, p);
        int32_t touched = 0;
        for (int32_t i = 0; i < length; i++) {
            for (int32_t c = 0; c < model->num_candidates; c++) {
                if (model->candidates[c] == edges[i]) {
                    model->cover_mask[c] |= (uint64_t)1 << bit;
                    touched++;
                }
            }
        }
        if (touched == 0) {
            return 1; /* an uncovered path without candidate edge: cannot happen */
        }
        bit++;
    }
    model->required = bit == 64 ? ~(uint64_t)0 : ((uint64_t)1 << bit) - 1;
    model->suffix_or[model->num_candidates] = 0;
    for (int32_t c = model->num_candidates - 1; c >= 0; c--) {
        model->suffix_or[c] = model->suffix_or[c + 1] | model->cover_mask[c];
    }
    return 0;
}

/**
 * @brief Depth-first search for a subset of @p target_size candidates covering
 *        every uncovered path.
 *
 * @param model       Model with its candidates and path masks.
 * @param first       First candidate index available to choose from.
 * @param chosen      Candidates already selected (count in @p depth).
 * @param depth       Size of @p chosen.
 * @param target_size Size the subset must reach.
 * @param covered     Union of the path masks of @p chosen.
 * @return true if a feasible subset was found (stored in model->solution).
 */
static bool search_subset(BruteForceModel *model, int32_t first, int32_t *chosen,
                          int32_t depth, int32_t target_size, uint64_t covered)
{
    if (covered == model->required) {
        for (int32_t i = 0; i < depth; i++) {
            model->solution[i] = chosen[i];
        }
        return true;
    }
    int32_t missing = target_size - depth;
    if (missing == 0 || model->num_candidates - first < missing) {
        return false;
    }
    if ((covered | model->suffix_or[first]) != model->required) {
        return false; /* the remaining candidates cannot cover all paths */
    }
    for (int32_t c = first; c <= model->num_candidates - missing; c++) {
        chosen[depth] = c;
        if (search_subset(model, c + 1, chosen, depth + 1, target_size,
                          covered | model->cover_mask[c])) {
            return true;
        }
    }
    return false;
}

/**
 * @brief Selects the edges of a feasible warm start, if it covers every
 *        uncovered path.
 *
 * @param cover      Cover to complete.
 * @param warm_start [num_edges] flags, or nullptr.
 * @return Number of edges added, or -1 if the warm start is absent or not feasible.
 */
static int32_t apply_warm_start(Cover *cover, const uint8_t *warm_start)
{
    if (warm_start == nullptr) {
        return -1;
    }
    const PathPool *pool = cover->pool;
    for (int32_t p = 0; p < cover->num_paths; p++) {
        if (cover->path_hits[p] > 0) {
            continue;
        }
        const int32_t *edges = path_pool_edges(pool, p);
        int32_t length = path_pool_length(pool, p);
        bool hit = false;
        for (int32_t i = 0; i < length && !hit; i++) {
            hit = warm_start[edges[i]] != 0 && !cover->selected[edges[i]];
        }
        if (!hit) {
            return -1;
        }
    }
    int32_t added = 0;
    for (int32_t e = 0; e < cover->graph->num_edges; e++) {
        if (warm_start[e] != 0 && !cover->selected[e]) {
            cover_add_edge(cover, e);
            added++;
        }
    }
    return added;
}

ExactResult exact_cover(ExactSolver *solver, Cover *cover, double time_limit,
                        const uint8_t *warm_start, double known_bound)
{
    (void)solver;
    (void)time_limit;
    ExactResult result = {
        .runtime = 0.0,
        .objective_bound = NAN,
        .has_solution = false,
        .solver_calls = 0,
    };
    double start = timer_now();

    if (cover->num_uncovered == 0) {
        result.objective_bound = 0.0;
        result.has_solution = true;
        return result; /* nothing to solve: the solver is not called */
    }

    BruteForceModel model;
    if (build_model(cover, &model, BF_MAX_CANDIDATES) != 0) {
        /* Too large for this stand-in: complete with the warm start, else greedily. */
        if (apply_warm_start(cover, warm_start) < 0) {
            Rng rng;
            rng_seed(&rng, 20261119);
            greedy_cover(cover, 0, &rng);
        }
        result.runtime = timer_now() - start;
        result.has_solution = true;
        result.solver_calls = 1;
        return result; /* objective_bound stays NAN: nothing was proved */
    }

    /* Enumerate sizes starting at the known lower bound, if one was given. */
    int32_t first_size = 1;
    if (known_bound > 0.0) {
        long long rounded = (long long)(known_bound + 0.5);
        if (rounded > first_size) {
            first_size = (int32_t)rounded;
        }
    }
    int32_t found_size = -1;
    for (int32_t size = first_size; size <= model.num_candidates; size++) {
        int32_t chosen[BF_MAX_CANDIDATES];
        if (search_subset(&model, 0, chosen, 0, size, 0)) {
            found_size = size;
            break;
        }
    }
    if (found_size < 0 && warm_start != nullptr) {
        found_size = apply_warm_start(cover, warm_start);
        if (found_size >= 0) {
            result.runtime = timer_now() - start;
            result.objective_bound = (double)found_size;
            result.has_solution = true;
            result.solver_calls = 1;
            return result;
        }
    }
    if (found_size >= 0) {
        for (int32_t i = 0; i < found_size; i++) {
            cover_add_edge(cover, model.candidates[model.solution[i]]);
        }
    }
    result.runtime = timer_now() - start;
    result.solver_calls = 1;
    result.objective_bound = found_size >= 0 ? (double)found_size : NAN;
    result.has_solution = found_size >= 0;
    return result;
}
