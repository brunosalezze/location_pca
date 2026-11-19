/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file util.h
 * @brief Small utilities shared by every module: fatal errors, logging,
 *        checked memory allocation and a growable integer array.
 *
 * Every allocation in the program goes through the x* wrappers below. They
 * never return a null pointer: when memory is exhausted (or a size
 * computation would overflow) the program stops with a clear message.
 * This keeps the algorithmic code free of repetitive error handling.
 */
#ifndef LOCATION_PCA_UTIL_H
#define LOCATION_PCA_UTIL_H

#include <stddef.h>
#include <stdint.h>

/* ------------------------------------------------------------------------- */
/* Logging                                                                   */
/* ------------------------------------------------------------------------- */

/** Verbosity levels, from least to most talkative. */
typedef enum LogLevel {
    LOG_QUIET = 0, /**< Only fatal errors are printed.                    */
    LOG_INFO = 1,  /**< Instance summary and one line per run (default).  */
    LOG_DEBUG = 2, /**< One line per iteration of the covering framework. */
    LOG_TRACE = 3, /**< Also forwards the MIP solver's own log.           */
} LogLevel;

/**
 * @brief Sets the global verbosity used by log_message().
 *
 * @param level New verbosity level. Messages with a level above it are dropped.
 */
void log_set_level(LogLevel level);

/**
 * @brief Returns the current global verbosity level.
 *
 * @return The level set by the last call to log_set_level() (LOG_INFO by default).
 */
LogLevel log_get_level(void);

/**
 * @brief Prints a printf-style message to stderr if @p level is enabled.
 *
 * A newline is appended automatically.
 *
 * @param level  Verbosity level of this message.
 * @param format printf-style format string, followed by its arguments.
 */
[[gnu::format(printf, 2, 3)]]
void log_message(LogLevel level, const char *format, ...);

/**
 * @brief Prints "location_pca: error: <message>" to stderr and exits with
 *        EXIT_FAILURE.
 *
 * Used for unrecoverable conditions: malformed input files, out of memory,
 * solver/license failures.
 *
 * @param format printf-style format string, followed by its arguments.
 */
[[noreturn, gnu::format(printf, 1, 2)]]
void fatal(const char *format, ...);

/* ------------------------------------------------------------------------- */
/* Checked memory allocation                                                 */
/* ------------------------------------------------------------------------- */

/**
 * @brief malloc() that never fails.
 *
 * @param size Number of bytes (0 is allowed and yields a valid, freeable pointer).
 * @return Pointer to uninitialized memory. Exits the program on failure.
 */
[[nodiscard]] void *xmalloc(size_t size);

/**
 * @brief calloc() that never fails.
 *
 * @param count Number of elements.
 * @param size  Size of each element in bytes.
 * @return Pointer to zero-initialized memory. Exits the program on failure.
 */
[[nodiscard]] void *xcalloc(size_t count, size_t size);

/**
 * @brief Allocates an uninitialized array of @p count elements of @p size bytes.
 *
 * The multiplication is overflow-checked, unlike a plain malloc(count * size).
 *
 * @param count Number of elements.
 * @param size  Size of each element in bytes.
 * @return Pointer to uninitialized memory. Exits the program on failure.
 */
[[nodiscard]] void *xmalloc_array(size_t count, size_t size);

/**
 * @brief Resizes an array to @p count elements of @p size bytes (overflow-checked).
 *
 * @param ptr   Existing allocation, or nullptr.
 * @param count New number of elements.
 * @param size  Size of each element in bytes.
 * @return Pointer to the resized block (contents preserved up to the old size).
 *         Exits the program on failure.
 */
[[nodiscard]] void *xrealloc_array(void *ptr, size_t count, size_t size);

/**
 * @brief strdup() that never fails.
 *
 * @param text NUL-terminated string to copy.
 * @return Newly allocated copy of @p text. Release it with free().
 */
[[nodiscard]] char *xstrdup(const char *text);

/* ------------------------------------------------------------------------- */
/* Growable array of int32_t                                                 */
/* ------------------------------------------------------------------------- */

/**
 * @brief A minimal dynamic array of int32_t ("vector").
 *
 * A zero-initialized IntVec ({0}) is a valid empty vector. Memory grows
 * geometrically, so a sequence of n pushes costs O(n) amortized time.
 */
typedef struct IntVec {
    int32_t *items;  /**< Elements [0, count). */
    size_t count;    /**< Number of stored elements. */
    size_t capacity; /**< Number of elements that fit without reallocating. */
} IntVec;

/**
 * @brief Ensures @p vec can hold at least @p min_capacity elements.
 *
 * @param vec          Vector to grow (contents are preserved).
 * @param min_capacity Required capacity.
 */
void intvec_reserve(IntVec *vec, size_t min_capacity);

/**
 * @brief Appends @p value at the end of @p vec.
 *
 * Inner workings: doubles the capacity when the vector is full.
 *
 * @param vec   Vector to modify.
 * @param value Value to append.
 */
static inline void intvec_push(IntVec *vec, int32_t value)
{
    if (vec->count == vec->capacity) {
        intvec_reserve(vec, vec->count + 1);
    }
    vec->items[vec->count++] = value;
}

/**
 * @brief Releases the memory of @p vec and resets it to the empty state.
 *
 * @param vec Vector to release (may already be empty).
 */
void intvec_free(IntVec *vec);

#endif /* LOCATION_PCA_UTIL_H */
