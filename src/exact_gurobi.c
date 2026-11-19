/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file exact_gurobi.c
 * @brief Exact covering method implemented with the Gurobi C API.
 *
 * Differences from the original implementation:
 *  - the Gurobi environment (license check) is created once per process,
 *    not once per solve;
 *  - every Gurobi call is checked; errors stop the program with Gurobi's message;
 *  - variables are found through an index map instead of a name lookup
 *    (GRBgetvarbyname) for every edge of every path;
 *  - constraints are added in batches (GRBaddconstrs);
 *  - variables are binary (the paper's model) instead of unbounded integers;
 *  - the solution is read only if the solver actually found one (the original
 *    read uninitialized values when the time limit hit before any incumbent);
 *  - the relative MIP gap is 0 instead of Gurobi's default 1e-4;
 *  - a lower bound proved by an earlier solve is passed on (BestObjStop), so
 *    the solver does not prove it again when the optimum has not changed.
 */
#include "exact.h"

#include "util.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gurobi_c.h"

/** Maximum number of nonzeros sent to Gurobi in one GRBaddconstrs() call. */
enum { CONSTRAINT_BATCH_NONZEROS = 1 << 16 };

struct ExactSolver {
    GRBenv *env;          /**< Started Gurobi environment. */
    bool show_solver_log; /**< Forward Gurobi's log to stderr. */
};

/**
 * @brief Stops the program if a Gurobi call failed.
 *
 * @param env   Environment used to fetch the error message.
 * @param error Return code of the Gurobi call (0 = success).
 * @param call  Text of the call, for the message.
 */
static void check_gurobi(GRBenv *env, int error, const char *call)
{
    if (error != 0) {
        fatal("Gurobi error %d in %s: %s", error, call, GRBgeterrormsg(env));
    }
}

/** Calls a Gurobi function and stops the program on error. */
#define GRB_CHECK(env, call) check_gurobi((env), (call), #call)

ExactSolver *exact_solver_create(int threads, bool show_solver_log)
{
    GRBenv *env = nullptr;
    int error = GRBemptyenv(&env);
    if (error != 0 || env == nullptr) {
        fatal("could not create the Gurobi environment (error %d)", error);
    }
    GRB_CHECK(env, GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0));
    GRB_CHECK(env, GRBsetintparam(env, GRB_INT_PAR_THREADS, threads));
    error = GRBstartenv(env);
    if (error != 0) {
        fatal("could not start Gurobi (error %d): %s\n"
              "Check that Gurobi is installed and that a valid license is available "
              "(see README.md, section 'Requirements').",
              error, GRBgeterrormsg(env));
    }
    ExactSolver *solver = xcalloc(1, sizeof *solver);
    solver->env = env;
    solver->show_solver_log = show_solver_log;
    return solver;
}

void exact_solver_free(ExactSolver *solver)
{
    if (solver == nullptr) {
        return;
    }
    GRBfreeenv(solver->env);
    free(solver);
}

const char *exact_solver_description(void)
{
    static char description[64];
    int major = 0;
    int minor = 0;
    int technical = 0;
    GRBversion(&major, &minor, &technical);
    snprintf(description, sizeof description, "Gurobi %d.%d.%d", major, minor, technical);
    return description;
}

/**
 * @brief Gurobi callback that copies the solver's log messages to stderr.
 *
 * @param model    Model being solved (unused).
 * @param cbdata   Callback data, used to fetch the message.
 * @param where    Callback location; only GRB_CB_MESSAGE is handled.
 * @param userdata Unused.
 * @return 0 (never interrupts the solver).
 */
static int __stdcall forward_solver_log(GRBmodel *model, void *cbdata, int where, void *userdata)
{
    (void)model;
    (void)userdata;
    if (where == GRB_CB_MESSAGE) {
        char *message = nullptr;
        if (GRBcbget(cbdata, where, GRB_CB_MSG_STRING, &message) == 0 && message != nullptr) {
            fputs(message, stderr);
        }
    }
    return 0;
}

/** Buffers of one batch of covering constraints in CSR form. */
typedef struct ConstraintBatch {
    int *row_start; /**< [capacity] Index of each row's first nonzero. */
    int *columns;   /**< [capacity] Variable index of each nonzero.    */
    double *ones;   /**< [capacity] Coefficients (all 1).              */
    char *senses;   /**< [capacity] All GRB_GREATER_EQUAL.             */
    double *rhs;    /**< [capacity] All 1.                             */
    int rows;       /**< Rows currently buffered.                      */
    int nonzeros;   /**< Nonzeros currently buffered.                  */
    int capacity;   /**< Maximum rows/nonzeros per batch.              */
} ConstraintBatch;

/**
 * @brief Allocates (or enlarges) a batch so that it holds @p capacity nonzeros.
 *
 * @param batch    Batch to (re)allocate; buffered content is preserved.
 * @param capacity New capacity (rows and nonzeros).
 */
static void batch_reserve(ConstraintBatch *batch, int capacity)
{
    if (capacity <= batch->capacity) {
        return;
    }
    size_t n = (size_t)capacity;
    batch->row_start = xrealloc_array(batch->row_start, n, sizeof *batch->row_start);
    batch->columns = xrealloc_array(batch->columns, n, sizeof *batch->columns);
    batch->ones = xrealloc_array(batch->ones, n, sizeof *batch->ones);
    batch->senses = xrealloc_array(batch->senses, n, sizeof *batch->senses);
    batch->rhs = xrealloc_array(batch->rhs, n, sizeof *batch->rhs);
    for (int i = batch->capacity; i < capacity; i++) {
        batch->ones[i] = 1.0;
        batch->senses[i] = GRB_GREATER_EQUAL;
        batch->rhs[i] = 1.0;
    }
    batch->capacity = capacity;
}

/**
 * @brief Sends the buffered rows to Gurobi and empties the batch.
 *
 * @param batch Batch to flush.
 * @param model Model receiving the constraints.
 * @param env   Environment (for error messages).
 */
static void batch_flush(ConstraintBatch *batch, GRBmodel *model, GRBenv *env)
{
    if (batch->rows > 0) {
        GRB_CHECK(env,
                  GRBaddconstrs(model, batch->rows, batch->nonzeros, batch->row_start,
                                batch->columns, batch->ones, batch->senses, batch->rhs, nullptr));
    }
    batch->rows = 0;
    batch->nonzeros = 0;
}

/**
 * @brief Releases the batch buffers.
 *
 * @param batch Batch to release.
 */
static void batch_free(ConstraintBatch *batch)
{
    free(batch->row_start);
    free(batch->columns);
    free(batch->ones);
    free(batch->senses);
    free(batch->rhs);
    *batch = (ConstraintBatch){0};
}

ExactResult exact_cover(ExactSolver *solver, Cover *cover, double time_limit,
                        const uint8_t *warm_start, double known_bound)
{
    ExactResult result = {.runtime = 0.0, .objective_bound = 0.0, .has_solution = true};
    if (cover->num_uncovered == 0) {
        return result; /* nothing to cover: adding zero edges is optimal */
    }
    const Graph *graph = cover->graph;
    const PathPool *pool = cover->pool;
    GRBenv *env = solver->env;

    /* 1. Variables: unselected edges lying on at least one uncovered path. */
    int *variable_of_edge = xmalloc_array((size_t)graph->num_edges, sizeof *variable_of_edge);
    int32_t *edge_of_variable = xmalloc_array((size_t)graph->num_edges, sizeof *edge_of_variable);
    int num_variables = 0;
    for (int32_t e = 0; e < graph->num_edges; e++) {
        variable_of_edge[e] = -1;
        if (!cover->selected[e] && cover->uncovered_paths_through[e] > 0) {
            variable_of_edge[e] = num_variables;
            edge_of_variable[num_variables++] = e;
        }
    }
    double *values = xmalloc_array((size_t)num_variables, sizeof *values); /* objective, then x */
    char *types = xmalloc_array((size_t)num_variables, sizeof *types);
    for (int k = 0; k < num_variables; k++) {
        values[k] = 1.0;
        types[k] = GRB_BINARY;
    }

    GRBmodel *model = nullptr;
    GRB_CHECK(env, GRBnewmodel(env, &model, "tclp_cover", num_variables, values, nullptr, nullptr,
                               types, nullptr));
    GRBenv *model_env = GRBgetenv(model);
    GRB_CHECK(env, GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE));
    GRB_CHECK(env,
              GRBsetdblparam(model_env, GRB_DBL_PAR_TIMELIMIT, time_limit > 0 ? time_limit : 0));
    /* The objective is an integer, so "optimal" must mean a closed gap. Gurobi's default
     * relative gap of 1e-4 would accept a cover one edge above the bound once the
     * objective exceeds 10 000 (possible on country-sized networks). */
    GRB_CHECK(env, GRBsetdblparam(model_env, GRB_DBL_PAR_MIPGAP, 0.0));
    if (known_bound > 0.0) {
        /* A solution with known_bound edges is optimal: stop there. The objective is
         * an integer, so "at most known_bound + 0.5" means exactly known_bound. */
        GRB_CHECK(env, GRBsetdblparam(model_env, GRB_DBL_PAR_BESTOBJSTOP, known_bound + 0.5));
    }
    if (solver->show_solver_log) {
        GRB_CHECK(env, GRBsetintparam(model_env, GRB_INT_PAR_OUTPUTFLAG, 1));
        GRB_CHECK(env, GRBsetintparam(model_env, GRB_INT_PAR_LOGTOCONSOLE, 0));
        GRB_CHECK(env, GRBsetcallbackfunc(model, forward_solver_log, nullptr));
    }

    /* 2. One covering constraint per uncovered path, sent in batches. */
    ConstraintBatch batch = {0};
    batch_reserve(&batch, CONSTRAINT_BATCH_NONZEROS);
    for (int32_t p = 0; p < pool->num_paths; p++) {
        if (cover->path_hits[p] != 0) {
            continue;
        }
        int32_t length = path_pool_length(pool, p);
        const int32_t *edges = path_pool_edges(pool, p);
        if (batch.nonzeros + length > batch.capacity) {
            batch_flush(&batch, model, env);
            batch_reserve(&batch, length); /* only grows for absurdly long paths */
        }
        batch.row_start[batch.rows++] = batch.nonzeros;
        for (int32_t i = 0; i < length; i++) {
            batch.columns[batch.nonzeros++] = variable_of_edge[edges[i]];
        }
    }
    batch_flush(&batch, model, env);
    batch_free(&batch);

    /* 3. Optional MIP start, then solve. */
    if (warm_start != nullptr) {
        for (int k = 0; k < num_variables; k++) {
            values[k] = warm_start[edge_of_variable[k]] ? 1.0 : 0.0;
        }
        GRB_CHECK(env, GRBsetdblattrarray(model, GRB_DBL_ATTR_START, 0, num_variables, values));
    }
    GRB_CHECK(env, GRBoptimize(model));
    result.solver_calls = 1;

    /* 4. Read the outcome. */
    int solution_count = 0;
    double bound = NAN;
    GRB_CHECK(env, GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &solution_count));
    GRB_CHECK(env, GRBgetdblattr(model, GRB_DBL_ATTR_RUNTIME, &result.runtime));
    if (GRBgetdblattr(model, GRB_DBL_ATTR_OBJBOUND, &bound) != 0 || !isfinite(bound) || bound < 0) {
        bound = NAN; /* not available (e.g. stopped before the root relaxation) */
    }
    result.objective_bound = bound;
    result.has_solution = solution_count > 0;

    if (result.has_solution) {
        GRB_CHECK(env, GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, num_variables, values));
        for (int k = 0; k < num_variables; k++) {
            if (values[k] > 0.5) {
                cover_add_edge(cover, edge_of_variable[k]);
            }
        }
    }

    GRBfreemodel(model);
    free(values);
    free(types);
    free(variable_of_edge);
    free(edge_of_variable);
    return result;
}
