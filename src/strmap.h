/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file strmap.h
 * @brief Hash map from strings to int32_t, used to turn node labels into indices.
 *
 * The original code searched the whole node list for every label read
 * (O(|N|) per lookup, O(|A|·|N|) to load an instance). This open-addressing
 * table makes each lookup O(1) on average.
 *
 * The map stores pointers to the keys but does not copy or free them: the
 * caller must keep every inserted key alive while the map is in use.
 */
#ifndef LOCATION_PCA_STRMAP_H
#define LOCATION_PCA_STRMAP_H

#include <stddef.h>
#include <stdint.h>

/** Open-addressing (linear probing) hash table. */
typedef struct StrMap {
    const char **keys; /**< [capacity] Key of each slot, nullptr if the slot is empty. */
    int32_t *values;   /**< [capacity] Value of each occupied slot. */
    size_t capacity;   /**< Number of slots (always a power of two). */
    size_t count;      /**< Number of occupied slots. */
} StrMap;

/**
 * @brief Initializes an empty map.
 *
 * @param map               Map to initialize.
 * @param expected_elements Hint used to size the table (it grows when needed).
 */
void strmap_init(StrMap *map, size_t expected_elements);

/**
 * @brief Releases the table (but not the keys, which belong to the caller).
 *
 * @param map Map to release.
 */
void strmap_free(StrMap *map);

/**
 * @brief Looks up @p key.
 *
 * @param map The map.
 * @param key NUL-terminated key.
 * @return The value stored for @p key, or -1 if the key is absent.
 */
int32_t strmap_get(const StrMap *map, const char *key);

/**
 * @brief Inserts @p key with @p value. The key must not already be present.
 *
 * Inner workings: doubles the table when it becomes more than half full,
 * which keeps probe sequences short.
 *
 * @param map   The map.
 * @param key   NUL-terminated key; the pointer is stored, not copied.
 * @param value Value to associate (non-negative by convention).
 */
void strmap_put(StrMap *map, const char *key, int32_t value);

#endif /* LOCATION_PCA_STRMAP_H */
