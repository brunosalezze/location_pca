import math
import random
import subprocess
import sys
from pathlib import Path

BIN = sys.argv[1]
FLAGS = sys.argv[2:]


def instance(n, edges, terminals, rotations):
    return "\n".join([
        f"{n} {len(edges)} {len(terminals)}",
        " ".join(map(str, terminals)),
        *(f"{u} {v} {w}" for u, v, w in edges),
        *(" ".join(map(str, [len(row), *row])) for row in rotations),
    ]) + "\n"


def embedded(points, edges):
    rows = [[] for _ in points]
    for i, (u, v, _) in enumerate(edges, 1):
        x, y = points[u]
        a, b = points[v]
        rows[u].append((math.atan2(b-y, a-x), i))
        rows[v].append((math.atan2(y-b, x-a), -i))
    return [[i for _, i in sorted(row)] for row in rows]


def separates(n, edges, terminals, cut):
    parent = list(range(n))

    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v

    for i, (u, v, _) in enumerate(edges, 1):
        if i not in cut:
            parent[root(u)] = root(v)
    return len({root(t) for t in terminals}) == len(terminals)


def oracle(n, edges, terminals):
    best = sum(w for _, _, w in edges)
    for mask in range(1 << len(edges)):
        cut = {i+1 for i in range(len(edges)) if mask >> i & 1}
        cost = sum(edges[i-1][2] for i in cut)
        if cost < best and separates(n, edges, terminals, cut):
            best = cost
    return best


def run(text):
    p = subprocess.run([BIN, *FLAGS], input=text, capture_output=True, text=True, timeout=30)
    assert p.returncode == 0, (p.stderr, text)
    header, size, *ids = p.stdout.splitlines()
    cut = list(map(int, ids))
    assert size == f"cut {len(ids)}" and cut == sorted(set(cut)), p.stdout
    assert header.startswith("cost "), p.stdout
    return int(header.split()[1]), cut


def check(n, edges, terminals, rows, expected=None):
    text = instance(n, edges, terminals, rows)
    cost, cut = run(text)
    if expected is None:
        expected = oracle(n, edges, terminals)
    assert all(1 <= i <= len(edges) for i in cut), (text, cut)
    assert cost == sum(edges[i-1][2] for i in cut), (text, cut, cost)
    assert separates(n, edges, terminals, cut), (text, cut)
    assert cost == expected, (text, cost, expected, cut)


def examples():
    fixture = Path(__file__).with_name("triangle.txt")
    p = subprocess.run([BIN, *FLAGS, str(fixture)], capture_output=True, text=True, check=True)
    assert p.stdout == "cost 9\ncut 3\n1\n2\n3\n"
    check(1, [], [0], [[]], 0)
    check(2, [], [1, 0], [[], []], 0)
    check(2, [(0, 1, 5)], [1, 0], [[1], [-1]], 5)
    check(3, [(0, 1, 2), (1, 2, 9)], [0, 2], [[1], [-1, 2], [-2]], 2)
    check(4, [(0, 1, 8), (0, 2, 3), (0, 3, 5)],
          [1, 2, 3], [[1, 2, 3], [-1], [-2], [-3]], 8)
    check(5, [(0, 1, 2), (2, 3, 7)], [3, 4, 0, 2, 1],
          [[1], [-1], [2], [-2], []], 9)
    check(2, [(0, 1, 3), (0, 1, 5), (0, 0, 19)], [0, 1],
          [[1, 3, -3, 2], [-1, -2]], 8)
    check(3, [(0, 1, 2**63), (1, 2, 2**63-2)], [0, 2],
          [[1], [-1, 2], [-2]], 2**63-2)
    square = [(0, 0), (1, 0), (1, 1), (0, 1)]
    cycle = [(0, 1, 7), (1, 2, 2), (2, 3, 5), (3, 0, 3)]
    check(4, cycle, [0, 2], embedded(square, cycle), 5)
    check(4, cycle, [0, 1, 2], embedded(square, cycle), 12)
    check(4, cycle, [], embedded(square, cycle), 0)
    wheel = cycle + [(4, v, 1+v) for v in range(4)]
    points = square + [(0.5, 0.5)]
    check(5, wheel, [0, 2, 4], embedded(points, wheel))
    check(5, wheel, [0, 1, 2, 3], embedded(points, wheel))
    # Each terminal keeps its paired nonterminal. The unique optimum
    # cuts both triangle rings, including edges incident to no terminal.
    prism_points = [(-3, -2), (3, -2), (0, 4), (-1, -1), (1, -1), (0, 1)]
    prism = [(i, (i+1) % 3, 1+i) for i in range(3)]
    prism += [(i+3, (i+1) % 3+3, 2+i) for i in range(3)]
    prism += [(i, i+3, 100) for i in range(3)]
    check(6, prism, [0, 1, 2], embedded(prism_points, prism), 15)
    # The maximal isolating sides overlap at the center of this unit star.
    # Pinning them sequentially must not merge any two terminal labels.
    check(4, [(0, v, 1) for v in (1, 2, 3)], [1, 2, 3],
          [[1, 2, 3], [-1], [-2], [-3]], 2)
    if not FLAGS:
        # The primal fast path labels a free hub with one of 21 fixed
        # terminals, beyond the Klein–Marx subset-mask limit.
        star = [(0, v, v) for v in range(1, 22)]
        check(22, star, list(range(1, 22)), [
            list(range(1, 22)), *([-v] for v in range(1, 22)),
        ], sum(range(1, 21)))


def random_grids():
    rng = random.Random(1047)
    points = [(i % 3, i // 3) for i in range(9)]
    possible = [(u, v) for u, (x, y) in enumerate(points)
                for v, (a, b) in enumerate(points)
                if u < v and abs(x-a)+abs(y-b) == 1]
    for i in range(100):
        edges = [(u, v, rng.randint(1, 12)) for u, v in possible
                 if i % 5 == 0 or rng.random() < 0.8]
        rng.shuffle(edges)
        terms = rng.sample(range(9), 2+i % 3)
        rotations = embedded(points, edges)
        if i % 2:
            rotations = [list(reversed(row)) for row in rotations]
        check(9, edges, terms, rotations)


def rejects():
    bad = [
        ("2 1 2\n0 0\n0 1 1\n1 1\n1 -1\n", "duplicate terminal"),
        ("2 1 2\n0 1\n0 1 0\n1 1\n1 -1\n", "positive"),
        ("2 1 2\n0 1\n0 1 1\n1 1\n0\n", "missing"),
        ("2 1 2\n0 1\n0 1 1\n2 1 1\n1 -1\n", "repeated"),
        ("2 1 2\n0 1\n0 1 1\n1 -1\n1 1\n", "wrong vertex"),
        ("1 0 0\n0\nx\n", "trailing"),
        (instance(2, [(0, 1, 1)]*3, [0, 1], [[1, 2, 3], [-1, -2, -3]]),
         "not planar"),
    ]
    for text, msg in bad:
        p = subprocess.run([BIN, *FLAGS], input=text, capture_output=True, text=True, timeout=5)
        assert p.returncode != 0 and msg in p.stderr, (p.stderr, text)


examples()
random_grids()
rejects()
print("passed fixed graphs, 100 seeded planar oracle cases, and invalid-input cases")
