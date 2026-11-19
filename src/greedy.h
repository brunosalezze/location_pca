/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file greedy.h
 * @brief Greedy covering heuristic (Algorithm 2) and the greedy part of the
 *        hybrid covering heuristic (Algorithm 3, lines 1-5).
 */
#ifndef LOCATION_PCA_GREEDY_H
#define LOCATION_PCA_GREEDY_H

#include "cover.h"
#include "rng.h"

#include <stdint.h>

/**
 * @brief Adds edges greedily until at most @p max_uncovered paths remain uncovered.
 *
 * Inner workings: at each step the edge contained in the largest number of
 * uncovered paths, l = argmax |P*_a|, is selected (ties broken uniformly at
 * random), which covers those paths. Cost per step: O(|A|) to find the
 * maximum plus the incremental counter update of cover_add_edge().
 *
 *  - Algorithm 2 (pure greedy): call with @p max_uncovered = 0.
 *  - Algorithm 3 (hybrid): call with @p max_uncovered = alpha·|P1|, then
 *    solve the remaining paths exactly.
 *
 * @param cover         Cover to extend (synchronized with its pool).
 * @param max_uncovered Stop as soon as cover->num_uncovered <= this value.
 * @param rng           Random generator used to break ties.
 */
void greedy_cover(Cover *cover, int32_t max_uncovered, Rng *rng);

#endif /* LOCATION_PCA_GREEDY_H */
