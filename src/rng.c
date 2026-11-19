/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file rng.c
 * @brief PCG32 (XSH-RR variant) pseudo-random number generator.
 */
#include "rng.h"

/** Multiplier of the underlying 64-bit linear congruential generator. */
#define PCG32_MULTIPLIER UINT64_C(6364136223846793005)

void rng_seed(Rng *rng, uint64_t seed)
{
    rng->state = 0;
    rng->increment = (seed << 1u) | 1u; /* stream must be odd */
    (void)rng_next(rng);
    rng->state += seed ^ UINT64_C(0x853c49e6748fea9b);
    (void)rng_next(rng);
}

uint32_t rng_next(Rng *rng)
{
    uint64_t old = rng->state;
    rng->state = old * PCG32_MULTIPLIER + rng->increment;
    uint32_t xorshifted = (uint32_t)(((old >> 18u) ^ old) >> 27u);
    uint32_t rotation = (uint32_t)(old >> 59u);
    return (xorshifted >> rotation) | (xorshifted << ((-rotation) & 31u));
}

uint32_t rng_below(Rng *rng, uint32_t bound)
{
    uint64_t product = (uint64_t)rng_next(rng) * bound;
    uint32_t low = (uint32_t)product;
    if (low < bound) {
        uint32_t threshold = (-bound) % bound;
        while (low < threshold) {
            product = (uint64_t)rng_next(rng) * bound;
            low = (uint32_t)product;
        }
    }
    return (uint32_t)(product >> 32u);
}
