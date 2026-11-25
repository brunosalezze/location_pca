"""Compare unmodified location_pca instances with -a planar and the exact solver.

Every instance keeps its complete, fixed O-D terminal set. The planar
algorithm must reproduce the reference optimum on planar networks and
refuse the nonplanar ones; cuts are re-verified for feasibility with the
original edge identifiers, parallel roads included.
"""

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


def reference_run(binary, source, timeout):
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
    cmd = [binary, "-a", "exact", "--alpha0", "1"]
    extra = len(expanded) - len(rows)
    if not extra:
        return run([*cmd, source], timeout)
    print(f"  reference: preserving {extra} parallel edge(s) by subdivision", flush=True)
    with tempfile.TemporaryDirectory(prefix=f"{source.name}-multigraph-") as tmp:
        case = Path(tmp)
        (case / "od_nodes.txt").symlink_to(source / "od_nodes.txt")
        (case / "edges.txt").write_text("".join(" ".join(row) + "\n" for row in expanded))
        return run([*cmd, case], timeout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path,
                        default=Path(__file__).resolve().parents[2] / "bin" / "location_pca")
    parser.add_argument("--instances", nargs="+", default=["AP", "RR", "AC", "AM"])
    parser.add_argument("--timeout", type=float, default=60, help="seconds per invocation")
    parser.add_argument("--klein-marx-only", action="store_true", help="disable the primal DP fast path")
    args = parser.parse_args()
    binary = args.binary.resolve()
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "1"))
    incomplete = False
    for name in args.instances:
        source = binary.parents[1] / "instances" / name
        terminals = {t for (t,) in records(source / "od_nodes.txt")}
        print(f"{name}: {len(terminals)} original terminals", flush=True)
        stage = "planarity check"
        try:
            plan = run([binary, "-a", "planar", "--planarity-only", source], args.timeout)
            verdict = next((line for line in plan.stderr.splitlines()
                            if line in {"planar yes", "planar no", "planar unknown"}), "")
            if verdict in {"planar no", "planar unknown"}:
                print(f"  {verdict}; no cut attempted", flush=True)
                incomplete |= verdict == "planar unknown"
                continue
            assert plan.returncode == 0 and verdict == "planar yes", plan.stderr
            print("  planar yes", flush=True)

            stage = "exact reference solver"
            reference = reference_run(binary, source, args.timeout)
            assert reference.returncode == 0, reference.stderr
            lines = reference.stdout.splitlines()
            header = next(i for i, row in enumerate(lines) if row.startswith("instance\t"))
            ref = dict(zip(lines[header].split("\t"), lines[header + 1].split("\t")))
            assert ref["feasible"] == "1" and int(ref["cost"]) == int(ref["lower_bound"]), ref
            print(f"  reference optimum: {ref['cost']}", flush=True)

            stage = "planar algorithm"
            flags = ["--klein-marx-only"] if args.klein_marx_only else []
            with tempfile.TemporaryDirectory(prefix=f"{name}-planar-") as tmp:
                result = run([binary, "-a", "planar", *flags, "--solution-dir", tmp, source],
                             args.timeout, env=env)
                if result.returncode and "terminals exceed the limit" in result.stderr:
                    print(f"  {result.stderr.splitlines()[-1]}; no cut computed", flush=True)
                    incomplete = True
                    continue
                assert result.returncode == 0, result.stderr
                assert "planar yes" in result.stderr.splitlines(), result.stderr
                lines = result.stdout.splitlines()
                header = next(i for i, row in enumerate(lines) if row.startswith("instance\t"))
                row = dict(zip(lines[header].split("\t"), lines[header + 1].split("\t")))
                assert row["algorithm"] == "planar" and row["feasible"] == "1", row
                value = int(row["cost"])
                assert int(row["lower_bound"]) == value, row
                solution = next(Path(tmp).glob("*.txt"))
                cut = {line for line in solution.read_text().splitlines() if not line.startswith("#")}
            assert len(cut) == value, (name, value, len(cut))
            verify(source / "edges.txt", terminals, cut)
            assert value == int(ref["cost"]), (name, value, ref)
            print(f"  feasible cut; exact cost {value}", flush=True)
        except subprocess.TimeoutExpired:
            print(f"  {stage} timed out after {args.timeout:g}s", flush=True)
            incomplete = True
    return int(incomplete)


if __name__ == "__main__":
    raise SystemExit(main())
