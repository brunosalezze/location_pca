/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file algorithms.c
 * @brief Set Covering Framework (Algorithm 1) and the methods built on it.
 *
 * Reading guide, from the paper's pseudo-code to this file:
 *  - Algorithm 1 (SCF) ................. covering_framework()
 *  - Covering phase (Alg. 2 / 3, exact)  covering_phase()
 *  - Cleaning phase .................... cleaning_phase()  -> cleaning.c
 *  - Cut phase (Alg. 5) ................ separation_phase() -> separation.c
 *  - Algorithm 4 (PCA) ................. run_pca()
 *
 * Bugs of the original implementation fixed here:
 *  - PCA kept the cover of an iteration interrupted by the time limit as
 *    "best" solution although it was infeasible (O-D paths still uncovered),
 *    so runs that hit TL2 could report a cost below the true optimum. Only
 *    covers for which the Cut phase finds no new path are accepted now.
 *  - PCA did not implement the final greedy completion (Algorithm 4, lines
 *    9-10) used when no feasible cover is found within the time limit.
 *  - When the solver stopped before computing a bound, the lower bound
 *    collapsed to 0 (and the gap became a division by zero).
 *  - Time spent outside the MIP solver (Cut phase, heuristics) and in the
 *    cleaning solves was not counted towards TL2. The limit is now on the
 *    total wall-clock time of the run.
 *
 * Additions to the article (see README.md, "Additions to the article"):
 *  - every exact solve gives a lower bound for the TCLP, also the ones of the
 *    hybrid covering (see ProvenBound);
 *  - a bound that is already known is handed to the solver, which then does
 *    not prove it a second time (see known_residual_bound());
 *  - a run stops as soon as its best cover is proven optimal;
 *  - the Cut phase looks only for direct paths (see separation.h).
 */
#include "algorithms.h"

#include "cover.h"
#include "greedy.h"
#include "path_pool.h"
#include "rng.h"
#include "timer.h"
#include "util.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/** Tolerance used when rounding the solver's bound up to an integer. */
#define BOUND_TOLERANCE 1e-6

/** Covering method used inside the framework. */
typedef enum CoveringMethod {
    COVERING_GREEDY, /**< Algorithm 2. */
    COVERING_HYBRID, /**< Algorithm 3 (alpha >= 1 means purely exact). */
} CoveringMethod;

/** Number of exact solves whose proven bounds are remembered (see ProvenBound). */
enum { BOUND_MEMORY_SIZE = 8 };

/**
 * A lower bound proved by one exact covering solve.
 *
 * The solve started from the cover `fixed` (the locked edges plus, in the
 * hybrid covering, the edges chosen by the greedy part) and proved that at
 * least `bound` more edges are needed to cover the paths P1 contained at that
 * time. Since P1 only grows during a run, two things follow:
 *  - Every feasible cover S of the TCLP has at least num_locked_edges + bound
 *    edges. S contains the locked edges, and the paths that `fixed` leaves
 *    uncovered contain no edge of `fixed`, so S covers them with at least
 *    `bound` edges outside `fixed`. Hence every exact solve gives a lower
 *    bound for the TCLP, not only the ones that cover all of P1.
 *  - A later solve starting from a subset of `fixed` leaves at least the same
 *    paths uncovered, so it needs at least `bound` edges as well. Near
 *    convergence the greedy part keeps choosing the same few edges, and the
 *    optimum of consecutive iterations is often the same.
 */
typedef struct ProvenBound {
    uint8_t *fixed;    /**< [num_edges] Selection flags when the solve started. */
    int32_t num_fixed; /**< Number of edges selected in `fixed`. */
    double bound;      /**< Further edges needed at least (an integer value). */
} ProvenBound;

/** Everything one run needs, bundled to keep function signatures short. */
typedef struct Run {
    const Graph *graph;
    const SolverParams *params;
    ExactSolver *solver;
    PathPool *pool;         /**< P1. */
    Cover *cover;           /**< Current cover C. */
    Separator *separator;   /**< Cut phase scratch memory. */
    Rng rng;                /**< Random generator of this run. */
    double start_time;      /**< timer_now() at the start of the run. */
    double lower_bound;     /**< Best lower bound found (NAN if none). */
    uint8_t *best_selected; /**< Best feasible cover found so far. */
    int32_t best_cost;      /**< Its cost (INT32_MAX if none yet). */
    RunResult *result;      /**< Statistics being accumulated. */

    ProvenBound proven[BOUND_MEMORY_SIZE]; /**< Most recent bounds (ring buffer). */
    int32_t num_proven;                    /**< Entries of `proven` in use. */
    int32_t next_proven;                   /**< Entry that is overwritten next. */
    uint8_t *fixed_before_solve;           /**< [num_edges] Scratch copy of the cover. */
} Run;

/* ------------------------------------------------------------------------- */
/* Parameters                                                                */
/* ------------------------------------------------------------------------- */

void solver_params_defaults(SolverParams *params, Algorithm algorithm)
{
    *params = (SolverParams){
        .algorithm = algorithm,
        .alpha0 = 1.0 / 16.0,
        .time_limit_total = 7200.0,
        .time_limit_iteration = 500.0,
        .time_limit_cleaning = 25.0,
        .spread = SPREAD_MLC,
        .cleaning = CLEANING_EXACT,
    };
    switch (algorithm) {
    case ALGORITHM_PCA:
        params->initial_paths_per_root = 0; /* all partners */
        params->paths_per_root = 5;
        break;
    case ALGORITHM_EXACT:
        params->initial_paths_per_root = 20;
        params->paths_per_root = 10;
        break;
    case ALGORITHM_GREEDY:
        params->initial_paths_per_root = 500;
        params->paths_per_root = 500;
        break;
    case ALGORITHM_PLANAR:
        params->initial_paths_per_root = 0; /* unused by the planar solver */
        params->paths_per_root = 0;
        break;
    }
}

const char *algorithm_name(Algorithm algorithm)
{
    switch (algorithm) {
    case ALGORITHM_PCA:
        return "pca";
    case ALGORITHM_EXACT:
        return "exact";
    case ALGORITHM_GREEDY:
        return "greedy";
    case ALGORITHM_PLANAR:
        return "planar";
    }
    return "?";
}

/* ------------------------------------------------------------------------- */
/* Time management                                                           */
/* ------------------------------------------------------------------------- */

/**
 * @brief Seconds left before the total time limit TL2.
 *
 * @param run Current run.
 * @return Remaining time (negative once the limit has passed).
 */
static double time_left(const Run *run)
{
    return run->params->time_limit_total - (timer_now() - run->start_time);
}

/**
 * @brief Time limit for one MIP solve: @p limit, capped by the time left if @p bounded.
 *
 * @param run     Current run.
 * @param limit   Per-solve limit (TL1 or the cleaning limit).
 * @param bounded Whether the total limit TL2 applies (false for the final fallback).
 * @return Seconds the solver may use (never negative).
 */
static double solve_time_limit(const Run *run, double limit, bool bounded)
{
    if (!bounded) {
        return limit;
    }
    double left = time_left(run);
    double capped = left < limit ? left : limit;
    return capped > 0.0 ? capped : 0.0;
}

/* ------------------------------------------------------------------------- */
/* The three phases of the framework                                         */
/* ------------------------------------------------------------------------- */

/**
 * @brief Best lower bound already known for the exact solve that is about to start.
 *
 * Inner workings: the largest of
 *  - the remembered bounds of the solves whose starting cover contains every
 *    edge selected now (second property of ProvenBound), and
 *  - lower_bound - |C|: the current cover C plus the edges the solve adds is
 *    a cover of P1 containing the locked edges, so it has at least
 *    lower_bound edges.
 * For a solve that covers all of P1 (nothing fixed but the locked edges) this
 * is simply the best bound found so far.
 *
 * @param run Current run; its cover holds the edges fixed before the solve.
 * @return Number of edges the solve has to add at least (0 if nothing is known).
 */
static double known_residual_bound(const Run *run)
{
    const Cover *cover = run->cover;
    int32_t num_edges = run->graph->num_edges;
    double known = isnan(run->lower_bound) ? 0.0 : run->lower_bound - cover->num_selected;
    for (int32_t k = 0; k < run->num_proven; k++) {
        const ProvenBound *entry = &run->proven[k];
        if (entry->bound <= known || entry->num_fixed < cover->num_selected) {
            continue; /* not better, or too few edges to contain the current cover */
        }
        bool contains_cover = true;
        for (int32_t e = 0; e < num_edges && contains_cover; e++) {
            contains_cover = cover->selected[e] <= entry->fixed[e];
        }
        if (contains_cover) {
            known = entry->bound;
        }
    }
    return known > 0.0 ? known : 0.0;
}

/**
 * @brief Records the bound proved by the exact solve that just finished.
 *
 * Inner workings: the bound is stored with the cover the solve started from
 * (run->fixed_before_solve), replacing the oldest of the BOUND_MEMORY_SIZE
 * remembered solves, or updating the entry of the same starting cover. The
 * run's lower bound rises to num_locked_edges + @p bound if that is better
 * (first property of ProvenBound).
 *
 * @param run       Current run.
 * @param num_fixed Number of edges selected when the solve started.
 * @param bound     Edges the solve had to add at least (an integer value >= 0).
 */
static void record_bound(Run *run, int32_t num_fixed, double bound)
{
    double tclp_bound = run->graph->num_locked_edges + bound;
    if (isnan(run->lower_bound) || tclp_bound > run->lower_bound) {
        run->lower_bound = tclp_bound;
    }
    if (bound <= 0.0) {
        return; /* nothing worth remembering */
    }
    size_t num_edges = (size_t)run->graph->num_edges;
    for (int32_t k = 0; k < run->num_proven; k++) {
        ProvenBound *entry = &run->proven[k];
        if (entry->num_fixed == num_fixed &&
            memcmp(entry->fixed, run->fixed_before_solve, num_edges) == 0) {
            entry->bound = fmax(entry->bound, bound);
            return;
        }
    }
    ProvenBound *entry = &run->proven[run->next_proven];
    run->next_proven = (run->next_proven + 1) % BOUND_MEMORY_SIZE;
    if (run->num_proven < BOUND_MEMORY_SIZE) {
        run->num_proven++;
    }
    memcpy(entry->fixed, run->fixed_before_solve, num_edges);
    entry->num_fixed = num_fixed;
    entry->bound = bound;
}

/**
 * @brief Covering phase: extends the (reset) cover until it covers all of P1.
 *
 * Inner workings:
 *  - COVERING_GREEDY: Algorithm 2.
 *  - COVERING_HYBRID: Algorithm 3. The greedy heuristic selects edges until at
 *    most alpha·|P1| paths are uncovered, then the MIP solver covers the rest.
 *    With alpha >= 1 the greedy part is skipped and the exact solve covers
 *    all of P1.
 *    If the solver finds no solution in time, the greedy heuristic completes
 *    the cover.
 *
 * Every exact solve yields a lower bound for the TCLP (see ProvenBound). The
 * solver is also told the best bound already known for the solve it starts
 * (see known_residual_bound()). A cover that attains that bound is optimal,
 * so the solver stops there instead of proving the same bound again. This
 * does not change the cover an iteration produces, only the time it takes.
 *
 * @param run     Current run (its cover was just reset).
 * @param method  Covering method.
 * @param alpha   Hybrid fraction (ignored by COVERING_GREEDY).
 * @param bounded Whether the total time limit applies.
 */
static void covering_phase(Run *run, CoveringMethod method, double alpha, bool bounded)
{
    Cover *cover = run->cover;
    if (method == COVERING_GREEDY) {
        greedy_cover(cover, 0, &run->rng);
        return;
    }

    if (alpha < 1.0) {
        greedy_cover(cover, (int32_t)(alpha * cover->num_paths), &run->rng);
    }
    int32_t fixed_edges = cover->num_selected;
    memcpy(run->fixed_before_solve, cover->selected, (size_t)run->graph->num_edges);
    double known_bound = known_residual_bound(run);
    double limit = solve_time_limit(run, run->params->time_limit_iteration, bounded);
    ExactResult solve = exact_cover(run->solver, cover, limit, nullptr, known_bound);
    run->result->time_solver += solve.runtime;
    run->result->exact_solves += solve.solver_calls;

    if (!isnan(solve.objective_bound)) {
        /* The objective is an integer, so the solver's bound is rounded up. It may
         * be below known_bound when the solver stopped early. */
        double solver_bound = ceil(solve.objective_bound - BOUND_TOLERANCE);
        record_bound(run, fixed_edges, fmax(known_bound, solver_bound));
    }
    if (!solve.has_solution) {
        log_message(LOG_DEBUG, "    no MIP solution within %.1fs; completing greedily", limit);
        greedy_cover(cover, 0, &run->rng);
    }
}

/**
 * @brief Cleaning phase: removes redundant counting stations (see cleaning.h).
 *
 * @param run     Current run.
 * @param bounded Whether the total time limit applies.
 */
static void cleaning_phase(Run *run, bool bounded)
{
    double limit = solve_time_limit(run, run->params->time_limit_cleaning, bounded);
    CleaningStats stats =
        clean_cover(run->cover, run->params->cleaning, run->solver, limit, &run->rng);
    run->result->time_solver += stats.solver_time;
    run->result->exact_solves += stats.solver_calls;
}

/**
 * @brief Cut phase: adds O-D paths that avoid the current cover to P1 (see separation.h).
 *
 * @param run            Current run.
 * @param paths_per_root Maximum new paths per O-D root (0 = no limit).
 * @return Number of new paths (0 = the cover is feasible).
 */
static int32_t separation_phase(Run *run, int32_t paths_per_root)
{
    double start = timer_now();
    int32_t new_paths = separator_run(run->separator, run->cover, run->pool, run->params->spread,
                                      paths_per_root, &run->rng);
    cover_sync(run->cover);
    run->result->time_separation += timer_now() - start;
    return new_paths;
}

/* ------------------------------------------------------------------------- */
/* Algorithm 1: Set Covering Framework                                       */
/* ------------------------------------------------------------------------- */

/**
 * @brief Keeps @p run's current cover as the best one if it is cheaper.
 *
 * @param run Current run; its cover must be feasible.
 */
static void offer_solution(Run *run)
{
    if (run->cover->num_selected < run->best_cost) {
        run->best_cost = run->cover->num_selected;
        memcpy(run->best_selected, run->cover->selected, (size_t)run->graph->num_edges);
    }
}

/** How a call of the Set Covering Framework ended. */
typedef enum FrameworkOutcome {
    FRAMEWORK_FEASIBLE,      /**< The Cut phase found no uncovered path: the cover is feasible. */
    FRAMEWORK_TIME_LIMIT,    /**< TL2 came first: the cover still leaves O-D paths uncovered. */
    FRAMEWORK_BOUND_REACHED, /**< The lower bound reached the best known cost; cover unfinished. */
} FrameworkOutcome;

/**
 * @brief Tells whether the best cover found so far is proven optimal.
 *
 * @param run Current run.
 * @return true if a feasible cover is known and its cost equals the lower bound.
 */
static bool is_proven_optimal(const Run *run)
{
    return run->best_cost != INT32_MAX && !isnan(run->lower_bound) &&
           (double)run->best_cost <= run->lower_bound;
}

/**
 * @brief Set Covering Framework (Algorithm 1).
 *
 * Repeats Covering -> Cleaning -> Cut on the growing path set P1 until the
 * Cut phase finds no uncovered O-D path, i.e. the cover is feasible.
 *
 * One addition to the article: the loop also stops as soon as an exact solve
 * proves a lower bound equal to the cost of the best feasible cover already
 * known. That cover is then optimal and no further iteration can improve it.
 * This typically happens in PCA's last rounds and does not change the cost a
 * run returns, only its running time.
 *
 * @param run     Current run (P1 already contains the initial paths).
 * @param method  Covering method.
 * @param alpha   Hybrid fraction for COVERING_HYBRID.
 * @param bounded Stop when the total time limit TL2 is reached.
 * @return How the loop ended. The current cover is feasible only for
 *         FRAMEWORK_FEASIBLE.
 */
static FrameworkOutcome covering_framework(Run *run, CoveringMethod method, double alpha,
                                           bool bounded)
{
    for (;;) {
        cover_reset(run->cover);
        covering_phase(run, method, alpha, bounded);
        run->result->iterations++;
        if (is_proven_optimal(run)) {
            log_message(LOG_DEBUG, "    iteration %3d: LB %.0f = best cost, proven optimal",
                        run->result->iterations, run->lower_bound);
            return FRAMEWORK_BOUND_REACHED;
        }
        cleaning_phase(run, bounded);
        int32_t cost = run->cover->num_selected;
        int32_t new_paths = separation_phase(run, run->params->paths_per_root);

        log_message(
            LOG_DEBUG, "    iteration %3d: cost %6d  new paths %7d  |P1| %8d  LB %6.0f  %8.2fs",
            run->result->iterations, cost, new_paths, run->pool->num_paths,
            isnan(run->lower_bound) ? -1.0 : run->lower_bound, timer_now() - run->start_time);

        if (new_paths == 0) {
            return FRAMEWORK_FEASIBLE;
        }
        if (bounded && time_left(run) <= 0.0) {
            return FRAMEWORK_TIME_LIMIT;
        }
    }
}

/**
 * @brief Builds the initial path set P1 with a Cut phase on the empty cover.
 *
 * @param run Current run (empty pool).
 */
static void initial_separation(Run *run)
{
    cover_reset(run->cover); /* only the locked edges */
    int32_t paths = separation_phase(run, run->params->initial_paths_per_root);
    log_message(LOG_DEBUG, "    initial Cut phase: %d paths", paths);
}

/* ------------------------------------------------------------------------- */
/* Methods                                                                   */
/* ------------------------------------------------------------------------- */

/**
 * @brief Progressive Covering Algorithm (Algorithm 4).
 *
 * The framework is solved with the hybrid covering for alpha = alpha0,
 * 2·alpha0, 4·alpha0, ... while alpha <= 1 and time remains. Each round starts
 * from the paths accumulated by the previous ones and returns a cover. The
 * cheapest feasible cover is kept. Small alphas let the greedy heuristic do
 * most of the work (fast). Larger alphas give the exact solver more freedom
 * (better covers). The rounds stop early once the best cover is proven
 * optimal (see covering_framework()).
 *
 * @param run Current run.
 */
static void run_pca(Run *run)
{
    initial_separation(run);
    for (double alpha = run->params->alpha0;
         alpha <= 1.0 && time_left(run) > 0.0 && !is_proven_optimal(run); alpha *= 2) {
        switch (covering_framework(run, COVERING_HYBRID, alpha, true)) {
        case FRAMEWORK_FEASIBLE:
            log_message(LOG_DEBUG, "  alpha %.4f: cost %d (feasible)", alpha,
                        run->cover->num_selected);
            offer_solution(run);
            break;
        case FRAMEWORK_TIME_LIMIT:
            log_message(LOG_DEBUG, "  alpha %.4f: interrupted by the time limit, cover discarded",
                        alpha);
            break;
        case FRAMEWORK_BOUND_REACHED:
            log_message(LOG_DEBUG, "  alpha %.4f: stopped, the best cover (cost %d) is optimal",
                        alpha, run->best_cost);
            break;
        }
    }
}

/**
 * @brief Set Covering Framework with the exact covering method.
 *
 * @param run Current run.
 */
static void run_exact(Run *run)
{
    initial_separation(run);
    if (covering_framework(run, COVERING_HYBRID, 1.0, true) == FRAMEWORK_FEASIBLE) {
        offer_solution(run);
    }
}

/**
 * @brief Set Covering Framework with the greedy covering method (no time limit).
 *
 * @param run Current run.
 */
static void run_greedy(Run *run)
{
    initial_separation(run);
    covering_framework(run, COVERING_GREEDY, 0.0, false);
    offer_solution(run);
}

void run_algorithm(const Graph *graph, const SolverParams *params, ExactSolver *solver,
                   uint64_t seed, RunResult *result)
{
    *result = (RunResult){.lower_bound = NAN};
    result->selected = xcalloc((size_t)graph->num_edges, sizeof *result->selected);

    PathPool *pool = path_pool_create(graph->num_edges);
    Run run = {
        .graph = graph,
        .params = params,
        .solver = solver,
        .pool = pool,
        .cover = cover_create(graph, pool),
        .separator = separator_create(graph),
        .start_time = timer_now(),
        .lower_bound = NAN,
        .best_selected = result->selected,
        .best_cost = INT32_MAX,
        .result = result,
    };
    rng_seed(&run.rng, seed);
    size_t num_edges = (size_t)graph->num_edges;
    run.fixed_before_solve = xmalloc_array(num_edges, sizeof *run.fixed_before_solve);
    for (int32_t k = 0; k < BOUND_MEMORY_SIZE; k++) {
        run.proven[k].fixed = xmalloc_array(num_edges, sizeof *run.proven[k].fixed);
    }

    switch (params->algorithm) {
    case ALGORITHM_PCA:
        run_pca(&run);
        break;
    case ALGORITHM_EXACT:
        run_exact(&run);
        break;
    case ALGORITHM_GREEDY:
        run_greedy(&run);
        break;
    case ALGORITHM_PLANAR:
        /* main.c dispatches the planar algorithm before reaching here. */
        fatal("the planar algorithm must not reach run_algorithm()");
    }

    if (run.best_cost == INT32_MAX) {
        /* Algorithm 4, lines 9-10: no feasible cover within TL2. */
        log_message(LOG_INFO, "  no feasible cover within the time limit: completing with the "
                              "greedy framework (no time limit)");
        covering_framework(&run, COVERING_GREEDY, 0.0, false);
        offer_solution(&run);
    }

    result->cost = run.best_cost;
    result->lower_bound = run.lower_bound;
    result->num_paths = pool->num_paths;
    result->connected_pairs = cover_count_connected_od_pairs(graph, result->selected);
    result->time_total = timer_now() - run.start_time;

    for (int32_t k = 0; k < BOUND_MEMORY_SIZE; k++) {
        free(run.proven[k].fixed);
    }
    free(run.fixed_before_solve);
    separator_free(run.separator);
    cover_free(run.cover);
    path_pool_free(pool);
}

void run_result_free(RunResult *result)
{
    free(result->selected);
    *result = (RunResult){0};
}
