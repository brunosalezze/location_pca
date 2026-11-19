/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file strmap.c
 * @brief Open-addressing string hash map (FNV-1a hash, linear probing).
 */
#include "strmap.h"

#include "util.h"

#include <stdlib.h>
#include <string.h>

/**
 * @brief 64-bit FNV-1a hash of a NUL-terminated string.
 *
 * @param text String to hash.
 * @return Hash value.
 */
static uint64_t hash_string(const char *text)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    for (const unsigned char *c = (const unsigned char *)text; *c != '\0'; c++) {
        hash ^= *c;
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

/**
 * @brief Returns the slot holding @p key, or the empty slot where it would go.
 *
 * @param keys     Slot array.
 * @param capacity Number of slots (power of two).
 * @param key      Key to locate.
 * @return Slot index.
 */
static size_t find_slot(const char *const *keys, size_t capacity, const char *key)
{
    size_t mask = capacity - 1;
    size_t slot = (size_t)hash_string(key) & mask;
    while (keys[slot] != nullptr && strcmp(keys[slot], key) != 0) {
        slot = (slot + 1) & mask;
    }
    return slot;
}

/**
 * @brief Rehashes every entry into a table with @p new_capacity slots.
 *
 * @param map          Map to grow.
 * @param new_capacity New number of slots (power of two, larger than count).
 */
static void grow(StrMap *map, size_t new_capacity)
{
    const char **keys = xcalloc(new_capacity, sizeof *keys);
    int32_t *values = xmalloc_array(new_capacity, sizeof *values);
    for (size_t i = 0; i < map->capacity; i++) {
        if (map->keys[i] != nullptr) {
            size_t slot = find_slot(keys, new_capacity, map->keys[i]);
            keys[slot] = map->keys[i];
            values[slot] = map->values[i];
        }
    }
    free(map->keys);
    free(map->values);
    map->keys = keys;
    map->values = values;
    map->capacity = new_capacity;
}

void strmap_init(StrMap *map, size_t expected_elements)
{
    size_t capacity = 16;
    while (capacity < 2 * expected_elements) {
        capacity *= 2;
    }
    map->keys = xcalloc(capacity, sizeof *map->keys);
    map->values = xmalloc_array(capacity, sizeof *map->values);
    map->capacity = capacity;
    map->count = 0;
}

void strmap_free(StrMap *map)
{
    free(map->keys);
    free(map->values);
    *map = (StrMap){0};
}

int32_t strmap_get(const StrMap *map, const char *key)
{
    size_t slot = find_slot(map->keys, map->capacity, key);
    return map->keys[slot] != nullptr ? map->values[slot] : -1;
}

void strmap_put(StrMap *map, const char *key, int32_t value)
{
    if (2 * (map->count + 1) > map->capacity) {
        grow(map, 2 * map->capacity);
    }
    size_t slot = find_slot(map->keys, map->capacity, key);
    if (map->keys[slot] == nullptr) {
        map->count++;
    }
    map->keys[slot] = key;
    map->values[slot] = value;
}
