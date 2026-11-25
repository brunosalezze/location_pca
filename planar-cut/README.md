# Planar terminal cut

An exact C23 solver for the terminal (multiway) cut problem on **planar**
graphs, offered here as an alternative to the Progressive Covering
Algorithm. It computes a minimum-cost set of edges whose removal separates
every pair of terminals, and verifies the returned cut on the original
graph. It includes the Klein–Marx algorithm from *Solving Planar
k-Terminal Cut in O(n^(c√k)) Time* (ICALP 2012), plus a faster exact
primal dynamic program that most road instances fall into.

The solver needs no Gurobi and no graph library. When an instance is not
planar it says so and stops; unlike the covering algorithm, it never
returns an infeasible or non-optimal cut on a planar instance.

## Build and test

From this directory:

```sh
make            # builds ./planar-cut
make test       # embedding tests + seeded oracle comparisons
./planar-cut tests/triangle.txt
```

Use a C23 compiler (tested with GCC 14). The default build enables OpenMP
(`make OPENMP=` builds without it).

## Road-network instances

The repository's `instances/` directories work directly:

```sh
./planar-cut --tclp ../instances/AP --planarity-only   # planarity check only
./planar-cut --tclp ../instances/AP                    # planarity + cut
```

`--tclp DIR` reads `edges.txt` (`edge_id node_a node_b lanes`) and
`od_nodes.txt` (one terminal ID per line). Blank lines and `#` comments
at token boundaries are allowed. Roads are undirected and have unit cost
regardless of the lane field. The output uses the original edge IDs.

The program checks planarity **before** trying a cut. `planar yes` means
it found and verified an embedding; `planar no` means it ruled one out;
`planar unknown` means the search exhausted its budget. Both `no` and
`unknown` stop the program without computing a cut.

On the bundled networks this build confirms 21 planar instances and rules
out Brasil, MG, PR, RJ and SP as nonplanar.

To compare against the Gurobi-based reference solver (built with
`make` in the repository root):

```sh
python3 tests/test_location.py --instances AP AC RR AM
```

The script checks planarity first, keeps every parallel road that the
reference loader would drop (by subdividing them in a temporary reference
input), and verifies planarity, cut feasibility, and agreement with the
reference optimum. `--timeout` sets the per-command time limit.

## How the default solver works

After planarity is verified, the solver tries a budgeted **exact primal
dynamic program** on the original graph. It fixes one distinct label per
terminal, then solves each connected piece of nonterminals independently,
using only the labels of its bordering terminals. It finds the
minimum-weight set of edges between different labels by eliminating the
nonterminals. This avoids topology enumeration when the remaining
elimination bags are small, without dropping any terminals. Before
forming the pieces it also uses minimum isolating cuts to fix labels when
doing so provably preserves at least one optimum. If the budget is
exceeded, the solver falls back to Klein–Marx topology enumeration
(up to 20 terminals per connected component).

Budgets are per piece: `CUT_DP_MIB` caps the estimated table space
(default 1145 MiB, maximum 60000 MiB), `CUT_DP_WORK` caps the estimated
work (default 100 billion units), and `CUT_DP_STATS=1` prints estimates
and budget rejections. `CUT_DP_MIB` is a table budget, not a process-wide
memory limit.

Use `--klein-marx-only` (on the solver or the comparison script) to
benchmark the paper implementation without the fast path.

## Results on the bundled instances

Unit-cost cuts, all original terminals, `gcc -O2`, one thread. The
Gurobi column is `bin/location_pca -a exact --alpha0 1` on this machine
(8 threads, 90 s limit); for `MA` it runs on a graph where the duplicated
road is merged, matching this solver's model.

| Instance | Terminals | Cut | Gurobi optimum | Status |
|---|---:|---:|---:|---|
| AP | 13 | 22 | 22 | certified match |
| AC | 20 | 30 | 30 | certified match |
| RR | 13 | 19 | 19 | certified match |
| AM | 37 | 39 | 39 | certified match |
| MS | 76 | 150 | 156 in 22 s | feasible; reference bound 145 at 90 s |
| GOeDF | 241 | 464 | – | feasible; not certified |
| AL | 97 | 137 | 137 | certified match |
| SE | 74 | 112 | 112 | certified match |
| RO | 34 | 88 | 88 | certified match |
| MA | 111 | 250 | 250 | certified match |
| ES | 214 | 144 | 144 | certified match |
| RN | 93 | 233 | 233 | certified match |
| PB | 133 | 296 | 296 | certified match |
| PE | 245 | 252 | 252 | certified match |
| BA | 253 | 630 | 630 (3874 s) | certified match |
| CE | 273 | 330 | 330 (1063 s) | certified match |
| PA | 187 | 175 | 175 (153 s) | certified match |
| MT | 140 | – | – | DP budget exceeded; declined |
| PI | 267 | 316 | 316 (2028 s) | certified match |
| TO | 106 | 231 | 231 (391 s) | certified match |
| SC | 464 | 371 | 371 (1804 s) | certified match |
| RS | 433 | 544 | 544 (4718 s) | certified match |

18 of the 21 planar instances match a Gurobi-certified optimum; 20 return
a cut. `MT` exceeds the dynamic-programming budget and the solver declines
rather than guess. On the larger instances the exact primal DP is usually
one to three orders of magnitude faster than the MIP: `RS` in 24 s versus
79 minutes, `PI` in 10 s versus 34 minutes.

A verification note: a reference "optimum" of 143 for `ES` turned out to
be a Gurobi time-limit misreport (the reported objective disagreed with
the printed cover); the true certified optimum is 144, matching this
solver.

## Weighted graphs with a supplied embedding

Run `./planar-cut graph.txt` or write the following whitespace-separated
format to standard input:

1. `n m k`: vertex, edge, and terminal counts.
2. `k` distinct terminal IDs between `0` and `n-1`.
3. `m` triples `u v cost`. Costs are positive integers with sum less
   than `UINT64_MAX`; edges are numbered from `1` in input order.
4. For each vertex: its degree, then the incident darts in cyclic order.
   Use `i` at edge `i`’s first endpoint and `-i` at its second. Use both
   darts for a loop and `0` for an isolated vertex.

The solver validates the supplied rotation system and supports parallel
edges, loops, and disconnected graphs.

## Implementation

`reader.c` reads road instances; `embed.c` embeds their biconnected
blocks by ear addition and validates the result with Euler’s formula.
`narrow.c` implements the primal dynamic program and the isolating-cut
reductions. `solve.c` combines cuts by the terminal-subset recurrence.
`graph.c` and `steiner.c` build and cut open the dual along a Steiner
tree. `topology.c` enumerates candidate topologies, and `realize.c` finds
minimum-cost realizations.
