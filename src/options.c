/* SPDX-FileCopyrightText: 2018-2026 The location_pca authors (see AUTHORS)
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Linking with the Gurobi Optimizer is additionally permitted; see LICENSE-EXCEPTION. */
/**
 * @file options.c
 * @brief Command-line parsing and validation.
 */
#define _POSIX_C_SOURCE 200809L

#include "options.h"

#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/** Exit status for command-line errors (as used by most Unix tools). */
enum { EXIT_USAGE = 2 };

/** Identifiers of the long options that have no short form. */
enum {
    OPT_ALPHA0 = 1000,
    OPT_SPREAD,
    OPT_CLEANING,
    OPT_INITIAL_PATHS,
    OPT_PATHS,
    OPT_ITERATION_TIME_LIMIT,
    OPT_CLEANING_TIME_LIMIT,
    OPT_THREADS,
    OPT_NAME,
    OPT_NO_HEADER,
    OPT_KEEP_DEGREE_TWO,
    OPT_INFO,
};

/**
 * @brief Prints the usage text.
 *
 * @param stream   Where to print (stdout for --help, stderr after an error).
 * @param program  Program name (argv[0]).
 */
static void print_usage(FILE *stream, const char *program)
{
    fprintf(
        stream,
        "Usage: %s [OPTIONS] INSTANCE_DIR\n"
        "       %s [OPTIONS] EDGES_FILE OD_FILE\n"
        "\n"
        "Solves the Traffic Counting Location Problem: chooses the fewest road links\n"
        "that must receive a counting station so that every path between every pair of\n"
        "origin-destination (O-D) nodes passes through at least one station.\n"
        "INSTANCE_DIR must contain 'edges.txt' and 'od_nodes.txt'.\n"
        "Results are printed to stdout as tab-separated values, logs go to stderr.\n"
        "\n"
        "Method:\n"
        "  -a, --algorithm NAME        pca (default), exact or greedy\n"
        "      --alpha0 X              PCA initial hybrid fraction, 0 < X <= 1 (default 0.0625)\n"
        "      --spread NAME           Cut phase spread criterion: mlc (default) or random\n"
        "      --cleaning NAME         exact (default), greedy, random or none\n"
        "      --initial-paths N       paths per O-D node in the first Cut phase\n"
        "                              (0 = all; default: pca 0, exact 20, greedy 500)\n"
        "      --paths N               paths per O-D node in later Cut phases\n"
        "                              (0 = all; default: pca 5, exact 10, greedy 500)\n"
        "      --keep-degree2          do not contract degree-2 intersections\n"
        "\n"
        "Limits:\n"
        "  -t, --time-limit S          total wall-clock limit of a run (default 7200)\n"
        "      --iteration-time-limit S  limit of one exact covering solve (default 500)\n"
        "      --cleaning-time-limit S   limit of one exact cleaning solve (default 25)\n"
        "      --threads N             MIP solver threads, 0 = automatic (default 4)\n"
        "\n"
        "Experiment and output:\n"
        "  -r, --runs N                independent runs (default 1)\n"
        "  -s, --seed N                seed of the first run; run k uses N+k\n"
        "                              (default: derived from the clock, always printed)\n"
        "      --name NAME             instance name in the output\n"
        "  -o, --solution-dir DIR      write the stations of every run to DIR\n"
        "      --no-header             do not print the header line\n"
        "      --info                  print instance statistics and exit\n"
        "  -v, --verbose               more logging (repeat for more)\n"
        "  -q, --quiet                 no logging\n"
        "  -h, --help                  show this help\n"
        "  -V, --version               show the version\n",
        program, program);
}

/**
 * @brief Prints an error about the command line and exits with status 2.
 *
 * @param program Program name.
 * @param message Description of the problem.
 * @param value   Offending value, or nullptr.
 */
[[noreturn]] static void usage_error(const char *program, const char *message, const char *value)
{
    if (value != nullptr) {
        fprintf(stderr, "%s: %s: '%s'\n", program, message, value);
    } else {
        fprintf(stderr, "%s: %s\n", program, message);
    }
    fprintf(stderr, "Try '%s --help' for more information.\n", program);
    exit(EXIT_USAGE);
}

/**
 * @brief Parses a whole string as a base-10 integer within [@p min, @p max].
 *
 * @param program Program name (for errors).
 * @param option  Option name (for errors).
 * @param text    Text to parse.
 * @param min     Smallest accepted value.
 * @param max     Largest accepted value.
 * @return The parsed value. Exits on invalid input.
 */
static long long parse_integer(const char *program, const char *option, const char *text,
                               long long min, long long max)
{
    char *end = nullptr;
    errno = 0;
    long long value = strtoll(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value < min || value > max) {
        char message[128];
        snprintf(message, sizeof message,
                 "invalid value for %s (expected an integer in [%lld, %lld])", option, min, max);
        usage_error(program, message, text);
    }
    return value;
}

/**
 * @brief Parses a whole string as a finite number within [@p min, @p max].
 *
 * @param program Program name (for errors).
 * @param option  Option name (for errors).
 * @param text    Text to parse.
 * @param min     Smallest accepted value.
 * @param max     Largest accepted value.
 * @return The parsed value. Exits on invalid input.
 */
static double parse_number(const char *program, const char *option, const char *text, double min,
                           double max)
{
    char *end = nullptr;
    errno = 0;
    double value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !isfinite(value) || value < min ||
        value > max) {
        char message[128];
        snprintf(message, sizeof message, "invalid value for %s (expected a number in [%g, %g])",
                 option, min, max);
        usage_error(program, message, text);
    }
    return value;
}

/**
 * @brief Maps a name to its index in a list of accepted names.
 *
 * @param program Program name (for errors).
 * @param option  Option name (for errors).
 * @param text    Name given by the user.
 * @param names   Accepted names, terminated by nullptr.
 * @return Index of @p text in @p names. Exits if it is not there.
 */
static int parse_choice(const char *program, const char *option, const char *text,
                        const char *const *names)
{
    for (int i = 0; names[i] != nullptr; i++) {
        if (strcmp(text, names[i]) == 0) {
            return i;
        }
    }
    char message[128];
    snprintf(message, sizeof message, "unknown value for %s", option);
    usage_error(program, message, text);
}

/**
 * @brief Joins a directory and a file name with a '/'.
 *
 * @param directory Directory path (a trailing '/' is accepted).
 * @param file      File name.
 * @return Newly allocated path.
 */
static char *join_path(const char *directory, const char *file)
{
    size_t length = strlen(directory);
    bool has_slash = length > 0 && directory[length - 1] == '/';
    size_t size = length + strlen(file) + 2;
    char *path = xmalloc(size);
    snprintf(path, size, "%s%s%s", directory, has_slash ? "" : "/", file);
    return path;
}

/**
 * @brief Derives an instance name from a path: last component, without extension.
 *
 * Examples: "instances/AC/" -> "AC", "data/net.txt" -> "net".
 *
 * @param path      File or directory path.
 * @param strip_ext Remove everything after the last '.'.
 * @return Newly allocated name.
 */
static char *name_from_path(const char *path, bool strip_ext)
{
    size_t end = strlen(path);
    while (end > 1 && path[end - 1] == '/') {
        end--;
    }
    size_t begin = end;
    while (begin > 0 && path[begin - 1] != '/') {
        begin--;
    }
    size_t stop = end;
    if (strip_ext) {
        for (size_t i = end; i > begin; i--) {
            if (path[i - 1] == '.') {
                stop = i - 1;
                break;
            }
        }
    }
    if (stop == begin) {
        stop = end;
    }
    char *name = xmalloc(stop - begin + 1);
    memcpy(name, path + begin, stop - begin);
    name[stop - begin] = '\0';
    return name;
}

void options_parse(int argc, char **argv, Options *options)
{
    static const char *const algorithm_names[] = {"pca", "exact", "greedy", nullptr};
    static const char *const spread_names[] = {"mlc", "random", nullptr};
    static const char *const cleaning_names[] = {"none", "random", "greedy", "exact", nullptr};
    static const struct option long_options[] = {
        {"algorithm", required_argument, nullptr, 'a'},
        {"alpha0", required_argument, nullptr, OPT_ALPHA0},
        {"spread", required_argument, nullptr, OPT_SPREAD},
        {"cleaning", required_argument, nullptr, OPT_CLEANING},
        {"initial-paths", required_argument, nullptr, OPT_INITIAL_PATHS},
        {"paths", required_argument, nullptr, OPT_PATHS},
        {"keep-degree2", no_argument, nullptr, OPT_KEEP_DEGREE_TWO},
        {"time-limit", required_argument, nullptr, 't'},
        {"iteration-time-limit", required_argument, nullptr, OPT_ITERATION_TIME_LIMIT},
        {"cleaning-time-limit", required_argument, nullptr, OPT_CLEANING_TIME_LIMIT},
        {"threads", required_argument, nullptr, OPT_THREADS},
        {"runs", required_argument, nullptr, 'r'},
        {"seed", required_argument, nullptr, 's'},
        {"name", required_argument, nullptr, OPT_NAME},
        {"solution-dir", required_argument, nullptr, 'o'},
        {"no-header", no_argument, nullptr, OPT_NO_HEADER},
        {"info", no_argument, nullptr, OPT_INFO},
        {"verbose", no_argument, nullptr, 'v'},
        {"quiet", no_argument, nullptr, 'q'},
        {"help", no_argument, nullptr, 'h'},
        {"version", no_argument, nullptr, 'V'},
        {nullptr, 0, nullptr, 0},
    };
    const char *program = argc > 0 ? argv[0] : "location_pca";

    *options = (Options){
        .threads = 4, /* as in the original implementation */
        .runs = 1,
        .print_header = true,
        .contract_degree_two = true,
        .verbosity = LOG_INFO,
    };

    /* Explicit values are collected first and applied over the algorithm's defaults. */
    Algorithm algorithm = ALGORITHM_PCA;
    double alpha0 = NAN;
    double time_limit = NAN;
    double iteration_time_limit = NAN;
    double cleaning_time_limit = NAN;
    long long initial_paths = -1;
    long long paths = -1;
    int spread = -1;
    int cleaning = -1;
    bool seed_given = false;
    const char *name = nullptr;
    int verbosity = LOG_INFO;

    int option;
    while ((option = getopt_long(argc, argv, "a:t:r:s:o:vqhV", long_options, nullptr)) != -1) {
        switch (option) {
        case 'a':
            algorithm = (Algorithm)parse_choice(program, "--algorithm", optarg, algorithm_names);
            break;
        case OPT_ALPHA0:
            alpha0 = parse_number(program, "--alpha0", optarg, 1e-6, 1.0);
            break;
        case OPT_SPREAD:
            spread = parse_choice(program, "--spread", optarg, spread_names);
            break;
        case OPT_CLEANING:
            cleaning = parse_choice(program, "--cleaning", optarg, cleaning_names);
            break;
        case OPT_INITIAL_PATHS:
            initial_paths = parse_integer(program, "--initial-paths", optarg, 0, INT32_MAX);
            break;
        case OPT_PATHS:
            paths = parse_integer(program, "--paths", optarg, 0, INT32_MAX);
            break;
        case OPT_KEEP_DEGREE_TWO:
            options->contract_degree_two = false;
            break;
        case 't':
            time_limit = parse_number(program, "--time-limit", optarg, 0.0, 1e9);
            break;
        case OPT_ITERATION_TIME_LIMIT:
            iteration_time_limit =
                parse_number(program, "--iteration-time-limit", optarg, 0.0, 1e9);
            break;
        case OPT_CLEANING_TIME_LIMIT:
            cleaning_time_limit = parse_number(program, "--cleaning-time-limit", optarg, 0.0, 1e9);
            break;
        case OPT_THREADS:
            options->threads = (int)parse_integer(program, "--threads", optarg, 0, 1024);
            break;
        case 'r':
            options->runs = (int32_t)parse_integer(program, "--runs", optarg, 1, INT32_MAX);
            break;
        case 's':
            options->seed = (uint64_t)parse_integer(program, "--seed", optarg, 0, INT64_MAX);
            seed_given = true;
            break;
        case OPT_NAME:
            name = optarg;
            break;
        case 'o':
            free(options->solution_dir); /* the option may be repeated; the last one wins */
            options->solution_dir = xstrdup(optarg);
            break;
        case OPT_NO_HEADER:
            options->print_header = false;
            break;
        case OPT_INFO:
            options->info_only = true;
            break;
        case 'v':
            verbosity++;
            break;
        case 'q':
            verbosity = LOG_QUIET;
            break;
        case 'h':
            print_usage(stdout, program);
            exit(EXIT_SUCCESS);
        case 'V':
            printf("location_pca %s\n", LOCATION_PCA_VERSION);
            exit(EXIT_SUCCESS);
        default:
            fprintf(stderr, "Try '%s --help' for more information.\n", program);
            exit(EXIT_USAGE);
        }
    }

    /* Positional arguments: INSTANCE_DIR, or EDGES_FILE OD_FILE. */
    int positional = argc - optind;
    if (positional == 1) {
        options->edges_path = join_path(argv[optind], "edges.txt");
        options->od_path = join_path(argv[optind], "od_nodes.txt");
        options->instance_name = name_from_path(argv[optind], false);
    } else if (positional == 2) {
        options->edges_path = xstrdup(argv[optind]);
        options->od_path = xstrdup(argv[optind + 1]);
        options->instance_name = name_from_path(argv[optind], true);
    } else {
        usage_error(program, positional == 0 ? "missing instance" : "too many arguments", nullptr);
    }
    if (name != nullptr) {
        free(options->instance_name);
        options->instance_name = xstrdup(name);
    }

    /* Algorithm defaults, then explicit overrides. */
    SolverParams *params = &options->params;
    solver_params_defaults(params, algorithm);
    if (!isnan(alpha0)) {
        params->alpha0 = alpha0;
    }
    if (!isnan(time_limit)) {
        params->time_limit_total = time_limit;
    }
    if (!isnan(iteration_time_limit)) {
        params->time_limit_iteration = iteration_time_limit;
    }
    if (!isnan(cleaning_time_limit)) {
        params->time_limit_cleaning = cleaning_time_limit;
    }
    if (initial_paths >= 0) {
        params->initial_paths_per_root = (int32_t)initial_paths;
    }
    if (paths >= 0) {
        params->paths_per_root = (int32_t)paths;
    }
    if (spread >= 0) {
        params->spread = (SpreadCriterion)spread;
    }
    if (cleaning >= 0) {
        params->cleaning = (CleaningMethod)cleaning;
    }

    if (!seed_given) {
        options->seed = (uint64_t)time(nullptr);
    }
    options->verbosity = (LogLevel)(verbosity > LOG_TRACE ? LOG_TRACE : verbosity);
}

void options_free(Options *options)
{
    free(options->edges_path);
    free(options->od_path);
    free(options->instance_name);
    free(options->solution_dir);
    *options = (Options){0};
}
