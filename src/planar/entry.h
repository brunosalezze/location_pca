/* SPDX-FileCopyrightText: 2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later */
/**
 * @file entry.h
 * @brief Interface from main.c to the planar algorithm (planar/entry.c).
 */
#ifndef LOCATION_PCA_PLANAR_ENTRY_H
#define LOCATION_PCA_PLANAR_ENTRY_H

#include "../options.h"

/**
 * @brief Runs the planar algorithm for every configured run.
 *
 * @param options Parsed options; edges_path and od_path select the input.
 * @return EXIT_SUCCESS, or EXIT_FAILURE if the network was not proven planar.
 */
int planar_run_all(const Options *options);

#endif /* LOCATION_PCA_PLANAR_ENTRY_H */
