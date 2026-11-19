/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file rng.h
 * @brief Small, fast and reproducible pseudo-random number generator (PCG32).
 *
 * The original implementation used the C library's rand(), whose sequence
 * differs between platforms and which is shared by the whole process. Every
 * run now owns an independent generator seeded explicitly, so a run can be
 * reproduced exactly from the seed printed in the results.
 *
 * Reference: M. E. O'Neill, "PCG: A Family of Simple Fast Space-Efficient
 * Statistically Good Algorithms for Random Number Generation" (2014).
 */
#ifndef LOCATION_PCA_RNG_H
#define LOCATION_PCA_RNG_H

#include <stdint.h>

/** State of a PCG32 generator. Initialize it with rng_seed(). */
typedef struct Rng {
    uint64_t state;     /**< Internal 64-bit state. */
    uint64_t increment; /**< Stream selector (always odd). */
} Rng;

/**
 * @brief Initializes @p rng from a 64-bit seed.
 *
 * Inner workings: the seed selects both the starting state and the stream,
 * then the generator is advanced once to mix the bits.
 *
 * @param rng  Generator to initialize.
 * @param seed Any 64-bit value; equal seeds yield equal sequences.
 */
void rng_seed(Rng *rng, uint64_t seed);

/**
 * @brief Returns the next 32 uniformly distributed random bits.
 *
 * @param rng Generator to advance.
 * @return A value in [0, 2^32).
 */
uint32_t rng_next(Rng *rng);

/**
 * @brief Returns an unbiased random integer in [0, @p bound).
 *
 * Inner workings: Lemire's multiply-and-reject method. Unlike `rand() % n`
 * it has no modulo bias.
 *
 * @param rng   Generator to advance.
 * @param bound Exclusive upper bound; must be greater than zero.
 * @return A uniformly distributed value in [0, bound).
 */
uint32_t rng_below(Rng *rng, uint32_t bound);

#endif /* LOCATION_PCA_RNG_H */
