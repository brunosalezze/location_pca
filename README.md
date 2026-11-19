# location_pca

Implementation of the **Progressive Covering Algorithm (PCA)** for the
**Traffic Counting Location Problem (TCLP)**, from the article

> B. S. Vieira, T. Ferrari, G. M. Ribeiro, L. Bahiense, R. D. Orrico Filho,
> C. A. Abramides, N. F. Rosa Campos Júnior. *A progressive hybrid set covering
> based algorithm for the traffic counting location problem.* Expert Systems
> with Applications 160 (2020) 113641.
> https://doi.org/10.1016/j.eswa.2020.113641

Given a road network and a set of origin-destination (O-D) nodes
(municipalities), the TCLP asks for the **smallest set of road links that must
receive a traffic counting station so that every path between every pair of
O-D nodes passes through at least one station**. The program solves it as a
set covering problem whose rows (O-D paths) are generated on demand. It
combines a greedy heuristic, an exact MIP model (Gurobi) and a hybrid of
both, as described in the article.

The repository contains the C code, the 26 Brazilian state instances used in the
article (plus the whole country), tests and scripts to reproduce the experiments.

---

## Contents

- [Requirements](#requirements)
- [Building](#building)
- [Quick start](#quick-start)
- [Usage](#usage)
- [Reproducing the article's experiments](#reproducing-the-articles-experiments)
- [How the code is organized](#how-the-code-is-organized)
- [Testing](#testing)
- [Changes with respect to the original code](#changes-with-respect-to-the-original-code)
- [Citing](#citing)
- [Contributing](#contributing)
- [Acknowledgements](#acknowledgements)
- [License](#license)

---

## Requirements

| Requirement | Notes |
|---|---|
| Linux or macOS | Any POSIX system should work; on Windows use WSL. |
| C compiler with **C23** support | GCC ≥ 13 or Clang ≥ 18. Tested with GCC 14.2, Clang 18.1 and Clang 23.1. |
| GNU Make | Any recent version. |
| **Gurobi** ≥ 9.0 and a license | Tested with Gurobi 13.0.3. Free licenses are available for academic use at <https://www.gurobi.com/academia/>. |
| Python 3 (optional) | Only for `scripts/summarize_results.py` (standard library only). |

Check the compiler with `gcc --version`. On Ubuntu 22.04 (GCC 11), install a
newer one with `sudo apt install gcc-13` and build with `make CC=gcc-13`.

**Gurobi setup.** Follow Gurobi's installation guide. On Linux it amounts to:

```sh
tar xvfz gurobi13.0.3_linux64.tar.gz -C /opt   # creates /opt/gurobi1303/linux64
export GUROBI_HOME=/opt/gurobi1303/linux64     # add these two lines to ~/.bashrc
export PATH="$GUROBI_HOME/bin:$PATH"
grbgetkey <your-license-key>                    # writes ~/gurobi.lic
gurobi_cl --version                             # check that it works
```

## Building

```sh
git clone <repository-url> location_pca
cd location_pca
make            # builds bin/location_pca
```

The Makefile finds Gurobi through `$GUROBI_HOME`, or else picks the newest
`/opt/gurobi*/linux64` (Linux) or `/Library/gurobi*/macos_universal2` (macOS).
It also detects the library name (`gurobi130`, `gurobi110`, ...) automatically.
Run `make help` to see what was detected.

| Command | Result |
|---|---|
| `make` | Optimized build, `bin/location_pca`. |
| `make NATIVE=1` | Same, tuned for this CPU (`-march=native`). |
| `make GUROBI_HOME=/path/to/gurobiXXXX/linux64` | Use a specific Gurobi installation. |
| `make debug` | `bin/location_pca-debug` with AddressSanitizer and UBSan. |
| `make test` | Unit tests (**no Gurobi needed**). |
| `make check` | Unit tests plus a regression run on 11 small instances (needs Gurobi). |
| `make clean` | Removes `build/` and `bin/`. |

The executable records Gurobi's library path (rpath), so it runs without
setting `LD_LIBRARY_PATH`.

### Troubleshooting

| Message | Fix |
|---|---|
| `error: Gurobi not found` | Set `GUROBI_HOME` or pass it to `make` (see above). |
| `could not start Gurobi (error 10009): No Gurobi license found` | Install a license (`grbgetkey`), or set `GRB_LICENSE_FILE` to the license file. |
| `does not support C23` | Use a newer compiler, e.g. `make CC=gcc-13` or `make CC=clang-18`. |
| `error while loading shared libraries: libgurobi130.so` | The Gurobi directory moved after building. Rebuild, or set `LD_LIBRARY_PATH=$GUROBI_HOME/lib`. |

## Quick start

```sh
bin/location_pca --seed 1 instances/AC
```

Log (stderr):

```text
instance AC: 61 nodes, 84 edges, 20 O-D nodes (190 pairs), 6 locked edges
  preprocessing: 114 edge lines; removed 0 self-loops, 0 parallel edges, 0 ring edges; contracted 30 degree-2 nodes
method: pca, alpha0 0.0625, spread mlc, cleaning exact, paths per O-D node 0 then 5 (0 = all)
limits: run 7200s, exact covering 500s, exact cleaning 25s; solver: Gurobi 13.0.3 (4 threads)
run 0 (seed 1): 30 counting stations, 0.02s
```

Results (stdout, tab-separated; aligned here with `column -t`):

```text
instance  run  seed  algorithm  cost  lower_bound  gap_pct  feasible  paths  iterations  mip_solves  time_total_s  time_mip_s  time_cut_s
AC        0    1     pca        30    30           0.00     1         153    23          24          0.022         0.017       0.000
```

AC needs 30 counting stations, and since the lower bound is also 30 the
solution is proven optimal. Logs go to stderr and the table to stdout, so
`bin/location_pca instances/AC > result.tsv` saves only the table.

To also get the selected links, add `-o solutions`. This writes
`solutions/AC_pca_seed1.txt`:

```text
# location_pca 1.0.0 solution
# instance: AC
# algorithm: pca
# seed: 1
# counting stations: 30
# feasible: yes
# lower bound: 30
# One edge identifier (first column of the edge file) per line.
13964
14005
...
```

## Usage

```text
location_pca [OPTIONS] INSTANCE_DIR
location_pca [OPTIONS] EDGES_FILE OD_FILE
```

`INSTANCE_DIR` contains `edges.txt` and `od_nodes.txt`. The format is described in
[instances/README.md](instances/README.md).

| Option | Default | Meaning |
|---|---|---|
| `-a, --algorithm NAME` | `pca` | `pca` (Algorithm 4), `exact` (framework with exact covering) or `greedy` (framework with greedy covering). |
| `--alpha0 X` | `0.0625` | PCA: initial hybrid fraction α₀ ∈ (0, 1]. `--alpha0 1` gives the pure exact setting that proves optimality. |
| `--spread NAME` | `mlc` | Cut phase spread criterion: `mlc` (Minimum Link Count) or `random` (ND). |
| `--cleaning NAME` | `exact` | Cleaning procedure: `exact` (EC), `greedy` (GC), `random` (RC) or `none`. |
| `--initial-paths N` | per method | Paths per O-D node in the first Cut phase (0 = all). PCA: 0, exact: 20, greedy: 500. |
| `--paths N` | per method | Paths per O-D node in later Cut phases. PCA: 5, exact: 10, greedy: 500. |
| `--keep-degree2` | off | Do not merge the two links of degree-2 intersections. |
| `-t, --time-limit S` | `7200` | Total wall-clock limit of one run (TL₂). |
| `--iteration-time-limit S` | `500` | Limit of each exact covering solve (TL₁). |
| `--cleaning-time-limit S` | `25` | Limit of each exact cleaning solve. |
| `--threads N` | `4` | Gurobi threads (`0` = Gurobi decides). Use `--threads 1` to measure single-thread times as in the article. |
| `-r, --runs N` | `1` | Independent runs (different random seeds). |
| `-s, --seed N` | clock | Seed of the first run. Run *k* uses *N + k*, so any run can be repeated exactly. |
| `--name NAME` | from path | Instance name used in the output. |
| `-o, --solution-dir DIR` | – | Write the selected links of every run to `DIR`. |
| `--no-header` | – | Do not print the table header (handy when appending to a file). |
| `--info` | – | Print instance statistics and exit. |
| `-v` / `-q` | – | More logging (`-v`: every iteration, `-vvv`: Gurobi's log) / no logging. |

**Result columns** (tab-separated, one line per run):

| Column | Meaning |
|---|---|
| `cost` | Number of counting stations in the best feasible solution. |
| `lower_bound` | Best proven lower bound on the optimum. Every exact covering solve of PCA and of the exact method contributes one. `NA` if there was none (the greedy method, or a time limit too short for a single solve). |
| `gap_pct` | 100 · (cost − lower_bound) / lower_bound. 0 means proven optimal. |
| `feasible` | 1 if the solution separates every O-D pair. It is checked independently with union-find at the end of each run and should always be 1. |
| `paths` | Number of O-D paths generated (rows of the set covering model). |
| `iterations`, `mip_solves` | Framework iterations and Gurobi calls (a covering or cleaning step with nothing left to cover does not call Gurobi and is not counted). |
| `time_total_s`, `time_mip_s`, `time_cut_s` | Wall-clock time of the run, inside Gurobi, and in the Cut phase. |

The exit status is 0 on success, 1 on errors (or if a run ended
infeasible, which should never happen), and 2 for invalid command-line arguments.

**Note on very large networks.** The first Cut phase of `instances/Brasil`
(5 235 O-D nodes, 13.7 million O-D pairs) finds 135 003 paths, and a short
run already returns a feasible cover. Lower `--initial-paths` if memory is
tight.

## Reproducing the article's experiments

`scripts/run_experiments.sh` runs a list of instances (all 26 by default)
and writes one TSV file. `scripts/summarize_results.py` turns it into a
per-instance summary (best, average and worst cost, gap, time). Use `-j 1`
when measuring time: parallel jobs compete for the CPU.


## How the code is organized

```text
location_pca/
├── Makefile
├── src/
│   ├── main.c            program flow: options → load → solve → report
│   ├── options.[ch]      command line
│   ├── graph.[ch]        road network (CSR adjacency), O-D pair assignment
│   ├── graph_load.c      reading the instance files and preprocessing
│   ├── strmap.[ch]       hash map from node labels to indices
│   ├── path_pool.[ch]    the set P1 of O-D paths (rows of the covering model)
│   ├── cover.[ch]        a cover C and its incremental coverage counters, feasibility check
│   ├── separation.[ch]   Cut phase: search trees with spread criteria (Algorithm 5)
│   ├── greedy.[ch]       greedy covering (Algorithm 2, first half of Algorithm 3)
│   ├── exact.h           exact covering interface (solver-independent)
│   ├── exact_gurobi.c    exact covering with Gurobi
│   ├── cleaning.[ch]     cleaning procedures RC / GC / EC
│   ├── algorithms.[ch]   set covering framework (Algorithm 1) and PCA (Algorithm 4)
│   ├── report.[ch]       result table and solution files
│   └── rng.[ch], timer.[ch], util.[ch]   random numbers, clock, allocation/logging
├── instances/            the 27 networks (see instances/README.md)
└── scripts/              batch experiments and result summaries
```

### From the article to the code

| Article | Code |
|---|---|
| Network G = (N, A), O-D nodes T | `Graph` in `graph.h` |
| Set of paths P₁ | `PathPool` in `path_pool.h` |
| Cover C, uncovered paths P\*, \|P\*ₐ\| | `Cover` in `cover.h` (`num_uncovered`, `uncovered_paths_through`) |
| Algorithm 1: Set Covering Framework | `covering_framework()` in `algorithms.c` |
| Algorithm 2: Greedy covering | `greedy_cover(cover, 0, rng)` |
| Algorithm 3: Hybrid covering | `covering_phase()` (greedy to α·\|P₁\|, then `exact_cover()`) |
| Exact covering, model (1)–(3) | `exact_cover()` in `exact_gurobi.c` |
| Cleaning RC / GC / EC, set R_C | `clean_cover()`, `cover_collect_redundant()` |
| Algorithm 5: Cut generation, spread criteria ND / MLC | `separator_run()`, `SpreadCriterion` |
| Algorithm 4: PCA | `run_pca()` in `algorithms.c` |

### Data flow of one run

```text
graph_load ──► initial Cut phase (cover = locked edges) ──► P1
                     ┌───────────────────────────────────────┐
                     ▼                                       │ new paths
   cover_reset ─► covering_phase ─► cleaning_phase ─► separation_phase
                                                            │ no new path
                                                            ▼
                                                   feasible cover C
```

PCA repeats this loop for α = α₀, 2α₀, …, 1 and keeps the cheapest feasible cover.

Every function is documented in its header (inputs, inner workings, output).
Start with `algorithms.c`, which reads like the pseudo-code of the article.

## Testing

- `make test` builds the unit tests with AddressSanitizer and UBSan and runs
  them. They replace Gurobi with a brute-force solver (`tests/exact_bruteforce.c`),
  so no license is needed (they also run in CI). Besides module tests, they solve
  150 random networks with every method and compare the result with an
  exhaustive search. They also check that preprocessing never changes the optimum.
- `make check` also runs `tests/regression.sh`, which solves 11 instances 3 times
  each and compares the result with the optimum proven in the article.

## Citing

If you use this code or the instances, please cite the article above (see also
[CITATION.cff](CITATION.cff)) and, for the instances, González et al. (2019),
listed in [instances/README.md](instances/README.md).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## Acknowledgements

The measurements reported under
[Additions to the article](#additions-to-the-article) and
[The six hardest instances](#the-six-hardest-instances) were carried out on
the SCIAMA High Performance Compute (HPC) cluster, which is supported by the
ICG and the University of Portsmouth.

## License

Copyright (C) 2018–2026 the location_pca authors, listed in [AUTHORS](AUTHORS).

- **Code** (everything outside `instances/`): GNU General Public License,
  version 3 or (at your option) any later version ([LICENSE](LICENSE)). You may
  use, study, change and redistribute it. If you distribute a modified
  version, it must also be released under the GPL, version 3 or later, with
  its source code.
- **Gurobi linking exception** ([LICENSE-EXCEPTION](LICENSE-EXCEPTION)): the
  authors additionally allow the code to be linked with the proprietary
  Gurobi Optimizer libraries and the result to be distributed. Without it,
  the GPL would not clearly allow distributing a compiled program that uses
  Gurobi.
- **Instances** (`instances/`): Creative Commons Attribution 4.0 International
  ([instances/LICENSE](instances/LICENSE)). You may share and adapt the data
  for any purpose if you credit the source as described in
  [instances/README.md](instances/README.md#license).

Each source file starts with `SPDX-FileCopyrightText` and
`SPDX-License-Identifier` lines naming its copyright holders and license.
