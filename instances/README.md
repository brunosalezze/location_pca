# Instances

Real road networks of the Brazilian states, built from the georeferenced
database of the National Department of Transport Infrastructure (DNIT) and
introduced by González et al. (2019) (see [Source](#source)). Each O-D node is a municipality and each
edge a segment of state or federal highway. `Brasil` (the whole country) is
included for completeness; it was not used in the paper.

Each instance is a directory with two plain-text files:

| File           | One line per | Columns                                     |
|----------------|--------------|---------------------------------------------|
| `edges.txt`    | road link    | `edge_id  node_a  node_b  lanes`            |
| `od_nodes.txt` | O-D node     | `node_id`                                   |

Columns are separated by whitespace. Links are undirected. The `lanes` column
is not used by the algorithms. Lines starting with `#` are comments. Identifiers
are arbitrary strings, and solution files report the `edge_id` of each selected
link.

To use your own network, create a directory with these two files and run
`bin/location_pca path/to/directory`, or pass the two files explicitly:
`bin/location_pca my_edges.txt my_od_nodes.txt`.

## Size after preprocessing

`raw edges` counts the lines of `edges.txt`. |N|, |A|, |T| and |W| are the
nodes, edges, O-D nodes and O-D pairs after preprocessing (see `graph_load()`
in `src/graph.h`). `locked` counts the edges joining two O-D nodes, which
belong to every solution. Regenerate with `bin/location_pca --info instances/<NAME>`.

| Instance | raw edges | \|N\| | \|A\| | \|T\| | \|W\| | locked |
|----------|----------:|------:|------:|------:|---------:|-------:|
| AC       |       114 |    61 |    84 |    20 |      190 |      6 |
| AL       |       242 |   169 |   219 |    97 |     4656 |     52 |
| AM       |       210 |    74 |    77 |    37 |      666 |      9 |
| AP       |       100 |    51 |    74 |    13 |       78 |      5 |
| BA       |      1271 |   810 |  1108 |   395 |    77815 |    242 |
| CE       |       725 |   399 |   606 |   177 |    15576 |    108 |
| ES       |       463 |   283 |   393 |    75 |     2775 |     19 |
| GOeDF    |      1295 |   797 |  1162 |   241 |    28920 |     65 |
| MA       |       406 |   257 |   354 |   163 |    13203 |    137 |
| MG       |      2548 |  1467 |  1905 |   803 |   322003 |    414 |
| MS       |       686 |   343 |   497 |    76 |     2850 |     16 |
| MT       |      1697 |   709 |  1066 |   140 |     9730 |     23 |
| PA       |       542 |   286 |   363 |   122 |     7381 |     49 |
| PB       |       527 |   384 |   479 |   213 |    22578 |    137 |
| PE       |       536 |   362 |   472 |   172 |    14706 |     73 |
| PI       |       631 |   405 |   550 |   212 |    22366 |    132 |
| PR       |      1321 |   779 |  1077 |   381 |    72390 |    209 |
| RJ       |      1069 |   502 |   719 |    86 |     3655 |     11 |
| RN       |       478 |   333 |   433 |   160 |    12720 |     83 |
| RO       |       324 |   183 |   254 |    50 |     1225 |     12 |
| RR       |       175 |    75 |    97 |    13 |       78 |      3 |
| RS       |      1019 |   675 |   859 |   391 |    76245 |    244 |
| SC       |       654 |   480 |   602 |   266 |    35245 |    175 |
| SE       |       284 |   183 |   245 |    74 |     2701 |     35 |
| SP       |      1861 |  1280 |  1683 |   606 |   183315 |    240 |
| TO       |       584 |   372 |   524 |   134 |     8911 |     56 |
| Brasil   |     20853 | 11927 | 16719 |  5235 | 13699995 |   2706 |

The node and edge counts can differ slightly from Table 1 of the paper, which
lists the sizes reported for the original database.

## Source

P. H. González, G. Clímaco, G. R. Mauri, B. S. Vieira, G. M. Ribeiro,
R. D. Orrico Filho, L. Simonetti, L. R. Perim, I. C. S. Hoffmann,
*New approaches for the traffic counting location problem*, Expert Systems
with Applications 132 (2019) 189–198. https://doi.org/10.1016/j.eswa.2019.04.068

The files are the originals with Windows line endings converted to Unix line
endings and a final newline added. The data are unchanged.

## License

The instances are licensed under the
[Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/)
(CC BY 4.0); the full text is in [LICENSE](LICENSE). You may copy,
redistribute and adapt them for any purpose, including commercially, provided
that you give appropriate credit, link to the license and indicate whether you
made changes. To give credit, cite González et al. (2019) above and the
article this repository implements (see [CITATION.cff](../CITATION.cff)).

The program code in this repository is licensed separately, under
GPL-3.0-or-later with a Gurobi linking exception (see [../LICENSE](../LICENSE)
and [../LICENSE-EXCEPTION](../LICENSE-EXCEPTION)).
