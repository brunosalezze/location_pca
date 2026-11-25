# The planar algorithm (`-a planar`)

An exact solver for the terminal (multiway) cut problem on **planar**
graphs, built into this program as `-a planar`. It computes a minimum-cost
set of roads whose removal separates every pair of O-D nodes, and verifies
the returned cut on the original graph. It includes the Klein–Marx
algorithm from *Solving Planar k-Terminal Cut in O(n^(c√k)) Time*
(ICALP 2012), plus a faster exact primal dynamic program that most road
instances fall into.

The algorithm needs no Gurobi and no graph library: it never calls a MIP
solver, and `make NO_GUROBI=1` builds a binary in which it (and
`-a greedy --cleaning none`) work fully. When an instance is not planar it
says so and stops; unlike the covering methods, it never returns an
infeasible or non-optimal cut on a planar instance.

## Usage and tests

From the repository root:

```sh
make                                            # builds bin/location_pca
bin/location_pca -a planar tests/planar/triangle.txt
make test                                       # includes these tests
```

The repository's `instances/` directories work directly:

```sh
bin/location_pca -a planar instances/AP                    # verdict on stderr, TSV row on stdout
bin/location_pca -a planar -o sol instances/AP             # cut edges written to sol/
bin/location_pca -a planar instances/AP --planarity-only   # planarity verdict only
```

The instance files are read untouched: `edges.txt` (`edge_id node_a node_b
lanes`) and `od_nodes.txt` (one terminal ID per line); blank lines and `#`
comments at token boundaries are allowed. Roads are undirected and have
unit cost regardless of the lane field, parallel roads are kept, and the
solution file written by `-o` lists the original edge identifiers. The
planarity verdict (`planar yes`/`no`/`unknown`) is logged to stderr so that
stdout keeps the TSV contract of the other algorithms; the planar row has
`lower_bound = cost` and `gap_pct = 0` because the cut is proven optimal.

The program checks planarity **before** trying a cut. `planar yes` means
it found and verified an embedding; `planar no` means it ruled one out;
`planar unknown` means the search exhausted its budget. Both `no` and
`unknown` stop the program without computing a cut. On the bundled
networks this build confirms 21 planar instances and rules out Brasil,
MG, PR, RJ and SP as nonplanar.

To compare against the Gurobi-based methods of the same program:

```sh
python3 tests/planar/test_location.py --instances AP AC RR AM
```

The script checks planarity first, keeps every parallel road that the
reference loader would drop (by subdividing them in a temporary reference
input), and verifies planarity, cut feasibility, and agreement with the
reference optimum. `--timeout` sets the per-command time limit and
`--klein-marx-only` benchmarks the paper algorithm without the fast path.

## How it works

After planarity is verified, the algorithm tries a budgeted **exact
primal dynamic program** on the original graph. It fixes one distinct
label per terminal, then solves each connected piece of nonterminals
independently, using only the labels of its bordering terminals. It finds
the minimum-weight set of edges between different labels by eliminating
the nonterminals. This avoids topology enumeration when the remaining
elimination bags are small, without dropping any terminals. Before
forming the pieces it also uses minimum isolating cuts to fix labels when
doing so provably preserves at least one optimum. If the budget is
exceeded, the solver falls back to Klein–Marx topology enumeration
(up to 20 terminals per connected component); `--klein-marx-only` forces
the fallback.

Budgets are per piece, set through the environment: `CUT_DP_MIB` caps the
estimated table space (default 1145 MiB, maximum 60000 MiB), `CUT_DP_WORK`
caps the estimated work (default 100 billion units), and `CUT_DP_STATS=1`
prints estimates and budget rejections. `CUT_DP_MIB` is a table budget,
not a process-wide memory limit.

## Results on the bundled instances

Unit-cost cuts, all original terminals, `gcc -O2`, one thread. The Gurobi
column is `bin/location_pca -a exact --alpha0 1` (8 threads, 90 s limit)
on a reference input that preserves parallel roads, so both columns solve
the same model.

`Terminals` is the count of distinct O-D nodes in `od_nodes.txt`.

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
| RO | 50 | 88 | 88 | certified match |
| MA | 163 | 251 | 251 | certified match |
| ES | 75 | 144 | 144 | certified match |
| RN | 160 | 233 | 233 | certified match |
| PB | 213 | 296 | 296 | certified match |
| PE | 172 | 252 | 252 | certified match |
| BA | 395 | 630 | 630 (3874 s) | certified match |
| CE | 177 | 330 | 330 (1063 s) | certified match |
| PA | 122 | 175 | 175 (153 s) | certified match |
| MT | 140 | – | – | DP budget exceeded; declined |
| PI | 212 | 316 | 316 (2028 s) | certified match |
| TO | 134 | 231 | 231 (391 s) | certified match |
| SC | 266 | 371 | 371 (1804 s) | certified match |
| RS | 391 | 544 | 544 (4718 s) | certified match |

18 of the 21 planar instances match a Gurobi-certified optimum; 20 return
a cut. `MT` exceeds the dynamic-programming budget and the solver declines
rather than guess. On the larger instances the exact primal DP is usually
one to three orders of magnitude faster than the MIP: `RS` in 24 s versus
79 minutes, `PI` in 10 s versus 34 minutes. The default `location_pca`
loader removes the duplicated `MA` road, which lowers the optimum by one
station; both columns above preserve it.

A verification note: a reference "optimum" of 143 for `ES` turned out to
be a Gurobi time-limit misreport (the reported objective disagreed with
the printed cover); the true certified optimum is 144, matching this
algorithm.

## Weighted graphs with a supplied embedding

For tests and experiments, `-a planar` also reads a general weighted
planar graph: `bin/location_pca -a planar graph.txt`, or the same format
on standard input with `-`. The whitespace-separated format is:

1. `n m k`: vertex, edge, and terminal counts.
2. `k` distinct terminal IDs between `0` and `n-1`.
3. `m` triples `u v cost`. Costs are positive integers with sum less
   than `UINT64_MAX`; edges are numbered from `1` in input order.
4. For each vertex: its degree, then the incident darts in cyclic order.
   Use `i` at edge `i`’s first endpoint and `-i` at its second. Use both
   darts for a loop and `0` for an isolated vertex.

The algorithm validates the supplied rotation system and supports parallel
edges, loops, and disconnected graphs.

## Implementation

`src/planar/entry.c` drives the algorithm and the input formats;
`reader.c` reads road instances; `embed.c` embeds the biconnected blocks
by ear addition and validates the result with Euler’s formula.
`narrow.c` implements the primal dynamic program and the isolating-cut
reductions. `solve.c` combines cuts by the terminal-subset recurrence.
`graph.c` and `steiner.c` build and cut open the dual along a Steiner
tree. `topology.c` enumerates candidate topologies, and `realize.c` finds
minimum-cost realizations.
