/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file timer.h
 * @brief Wall-clock time measurement.
 */
#ifndef LOCATION_PCA_TIMER_H
#define LOCATION_PCA_TIMER_H

/**
 * @brief Returns a monotonic wall-clock time stamp in seconds.
 *
 * Only differences between two calls are meaningful. The clock is not
 * affected by changes to the system date.
 *
 * @return Seconds elapsed since an arbitrary fixed point in the past.
 */
double timer_now(void);

#endif /* LOCATION_PCA_TIMER_H */
