"""Compare unmodified location_pca instances with its exact Gurobi solver."""

import argparse
import os
from pathlib import Path
import subprocess
import tempfile


def run(args, timeout, **kwargs):
    return subprocess.run(args, text=True, capture_output=True,
                          timeout=timeout, **kwargs)


def verify(edges, terminals, cut):
    parent = {}

    def root(v):
        parent.setdefault(v, v)
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    ids = set()
    for edge, a, b, _ in records(edges):
        assert edge not in ids, edge
        ids.add(edge)
        if edge not in cut:
            parent[root(a)] = root(b)
    assert cut <= ids
    assert len({root(v) for v in terminals}) == len(terminals)


def records(path):
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields and not fields[0].startswith("#"):
            yield fields[:next((i for i, f in enumerate(fields) if f.startswith("#")), len(fields))]


def reference_run(repo, source, timeout):
    """Keep the multigraph objective despite the reference loader's deduplication.

    Subdivide each extra parallel edge with its own fresh nonterminal.
    The unit-cost path can be cut for one station, just like the original edge.
    """
    rows = list(records(source / "edges.txt"))
    keys = set()
    expanded = []
    used = {key for edge, a, b, _ in rows for key in (edge, a, b)}
    used.update(t for (t,) in records(source / "od_nodes.txt"))
    serial = 0

    def fresh():
        nonlocal serial
        while True:
            key = f"__multicut_parallel_{serial}"
            serial += 1
            if key not in used:
                used.add(key)
                return key

    for edge, a, b, lanes in rows:
        key = tuple(sorted((a, b)))
        if a != b and key in keys:
            mid = fresh()
            expanded.extend([(edge, a, mid, lanes), (fresh(), mid, b, lanes)])
        else:
            expanded.append((edge, a, b, lanes))
        keys.add(key)
    cmd = [repo / "bin/location_pca", "-a", "exact", "--alpha0", "1"]
    extra = len(expanded) - len(rows)
    if not extra:
        return run([*cmd, source], timeout, cwd=repo)
    print(f"  reference: preserving {extra} parallel edge(s) by subdivision", flush=True)
    with tempfile.TemporaryDirectory(prefix=f"{source.name}-multigraph-") as tmp:
        case = Path(tmp)
        (case / "od_nodes.txt").symlink_to(source / "od_nodes.txt")
        (case / "edges.txt").write_text("".join(" ".join(row) + "\n" for row in expanded))
        return run([*cmd, case], timeout, cwd=repo)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, default=Path("planar-cut"))
    parser.add_argument("--repository", type=Path, default=Path(".."))
    parser.add_argument("--instances", nargs="+", default=["AP", "RR", "AC", "AM"])
    parser.add_argument("--timeout", type=float, default=60, help="seconds per invocation")
    parser.add_argument("--klein-marx-only", action="store_true", help="disable the primal DP fast path")
    args = parser.parse_args()
    solver = args.solver.resolve()
    repo = args.repository.resolve()
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "1"))
    incomplete = False
    for name in args.instances:
        source = repo / "instances" / name
        terminals = {t for (t,) in records(source / "od_nodes.txt")}
        print(f"{name}: {len(terminals)} original terminals", flush=True)
        stage = "planarity check"
        try:
            plan = run([solver, "--tclp", source, "--planarity-only"], args.timeout)
            verdict = plan.stdout.strip()
            if verdict in {"planar no", "planar unknown"}:
                print(f"  {verdict}; no cut attempted", flush=True)
                incomplete |= verdict == "planar unknown"
                continue
            assert plan.returncode == 0 and verdict == "planar yes", plan.stderr
            print("  planar yes", flush=True)

            stage = "exact reference solver"
            reference = reference_run(repo, source, args.timeout)
            assert reference.returncode == 0, reference.stderr
            lines = reference.stdout.splitlines()
            header = next(i for i, row in enumerate(lines) if row.startswith("instance\t"))
            ref = dict(zip(lines[header].split("\t"), lines[header + 1].split("\t")))
            assert ref["feasible"] == "1" and int(ref["cost"]) == int(ref["lower_bound"]), ref
            print(f"  reference optimum: {ref['cost']}", flush=True)

            stage = "planar cut solver"
            flags = ["--klein-marx-only"] if args.klein_marx_only else []
            result = run([solver, *flags, "--tclp", source], args.timeout, env=env)
            if result.returncode and "terminals exceed the limit" in result.stderr:
                print(f"  {result.stderr.splitlines()[-1]}; no cut computed", flush=True)
                incomplete = True
                continue
            assert result.returncode == 0, result.stderr
            verdict, cost, count, *cut = result.stdout.splitlines()
            assert verdict == "planar yes"
            assert cost.startswith("cost ") and count.startswith("cut ")
            value = int(cost.split()[1])
            assert value == int(count.split()[1]) == len(cut) == len(set(cut))
            verify(source / "edges.txt", terminals, set(cut))
            assert value == int(ref["cost"]), (name, value, ref)
            print(f"  feasible cut; exact cost {value}", flush=True)
        except subprocess.TimeoutExpired:
            print(f"  {stage} timed out after {args.timeout:g}s", flush=True)
            incomplete = True
    return int(incomplete)


if __name__ == "__main__":
    raise SystemExit(main())
