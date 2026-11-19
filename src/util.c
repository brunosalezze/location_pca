/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file util.c
 * @brief Implementation of logging, fatal errors, checked allocation and IntVec.
 */
#include "util.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if __has_include(<stdckdint.h>)
#include <stdckdint.h> /* C23 checked integer arithmetic */
#else
#define ckd_mul(result, a, b) __builtin_mul_overflow((a), (b), (result))
#endif

/** Current verbosity; a single process-wide setting configured once by main(). */
static LogLevel g_log_level = LOG_INFO;

void log_set_level(LogLevel level)
{
    g_log_level = level;
}

LogLevel log_get_level(void)
{
    return g_log_level;
}

void log_message(LogLevel level, const char *format, ...)
{
    if (level > g_log_level) {
        return;
    }
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fputc('\n', stderr);
}

void fatal(const char *format, ...)
{
    fflush(stdout);
    fputs("location_pca: error: ", stderr);
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fputc('\n', stderr);
    exit(EXIT_FAILURE);
}

/**
 * @brief Multiplies two sizes, stopping the program if the product overflows.
 *
 * @param count Number of elements.
 * @param size  Element size in bytes.
 * @return count * size.
 */
static size_t checked_size(size_t count, size_t size)
{
    size_t bytes;
    if (ckd_mul(&bytes, count, size)) {
        fatal("allocation size overflow (%zu x %zu bytes)", count, size);
    }
    return bytes;
}

void *xmalloc(size_t size)
{
    void *ptr = malloc(size > 0 ? size : 1);
    if (ptr == nullptr) {
        fatal("out of memory (requested %zu bytes)", size);
    }
    return ptr;
}

void *xcalloc(size_t count, size_t size)
{
    void *ptr = calloc(count > 0 ? count : 1, size > 0 ? size : 1);
    if (ptr == nullptr) {
        fatal("out of memory (requested %zu x %zu bytes)", count, size);
    }
    return ptr;
}

void *xmalloc_array(size_t count, size_t size)
{
    return xmalloc(checked_size(count, size));
}

void *xrealloc_array(void *ptr, size_t count, size_t size)
{
    size_t bytes = checked_size(count, size);
    void *resized = realloc(ptr, bytes > 0 ? bytes : 1);
    if (resized == nullptr) {
        fatal("out of memory (requested %zu bytes)", bytes);
    }
    return resized;
}

char *xstrdup(const char *text)
{
    size_t length = strlen(text) + 1;
    char *copy = xmalloc(length);
    memcpy(copy, text, length);
    return copy;
}

void intvec_reserve(IntVec *vec, size_t min_capacity)
{
    if (vec->capacity >= min_capacity) {
        return;
    }
    size_t capacity = vec->capacity > 0 ? vec->capacity : 4;
    while (capacity < min_capacity) {
        capacity *= 2;
    }
    vec->items = xrealloc_array(vec->items, capacity, sizeof *vec->items);
    vec->capacity = capacity;
}

void intvec_free(IntVec *vec)
{
    free(vec->items);
    *vec = (IntVec){0};
}
