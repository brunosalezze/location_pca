/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file cleaning.h
 * @brief Cleaning procedures (Section 4.1.1): removing redundant counting stations.
 *
 * After a covering step some selected edges may be redundant: every path
 * they cover is also covered by another selected edge. The redundant set is
 * R_C = { a ∈ C : C \ {a} still covers P1 } (see cover_collect_redundant()).
 */
#ifndef LOCATION_PCA_CLEANING_H
#define LOCATION_PCA_CLEANING_H

#include "cover.h"
#include "exact.h"
#include "rng.h"

/** The cleaning variants compared in Table 2 of the paper. */
typedef enum CleaningMethod {
    CLEANING_NONE,   /**< Keep the cover as it is. */
    CLEANING_RANDOM, /**< RC: remove one random redundant edge at a time until none is left. */
    CLEANING_GREEDY, /**< GC: remove all of R_C, then re-cover greedily. */
    CLEANING_EXACT,  /**< EC: remove all of R_C, then re-cover exactly (paper's choice). */
} CleaningMethod;

/** MIP solver usage of one call to clean_cover(). */
typedef struct CleaningStats {
    double solver_time;   /**< Seconds spent in the MIP solver. */
    int32_t solver_calls; /**< Number of MIP solves (0 or 1). */
} CleaningStats;

/**
 * @brief Removes redundant edges from @p cover with the chosen method.
 *
 * Inner workings:
 *  - RC: repeatedly computes R_C and removes one of its edges chosen
 *    uniformly at random, until R_C is empty;
 *  - GC: removes every edge of R_C at once (which may uncover some paths)
 *    and runs greedy_cover() on what became uncovered;
 *  - EC: same, but the uncovered paths are covered with exact_cover(). R_C
 *    itself is passed as a warm start, so the result is never worse than
 *    the input. If the solver finds nothing in time, R_C is put back.
 *
 * The cover stays feasible for P1 in every case. RC and EC never make it
 * larger. GC usually shrinks it too, but has no such guarantee: the greedy
 * re-covering may need more edges than the |R_C| it replaces.
 *
 * @param cover      Cover to clean (must cover all of P1).
 * @param method     Cleaning method.
 * @param solver     MIP solver (used by CLEANING_EXACT only; may be nullptr otherwise).
 * @param time_limit Time limit in seconds for the exact re-covering.
 * @param rng        Random generator.
 * @return How much the MIP solver was used (zero for the non-exact methods).
 */
CleaningStats clean_cover(Cover *cover, CleaningMethod method, ExactSolver *solver,
                          double time_limit, Rng *rng);

#endif /* LOCATION_PCA_CLEANING_H */
