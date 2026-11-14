"""
census.py -- analysis of the GAP census logs (manul.log, ali.log, ...).

Each line "U|n|i|e|K|M|sym|lam|dims|eig|kers|sizes" written by the GAP
scripts describes a uniform structure on G = SmallGroup(n, i): the e-th
epimorphism Z_K^{*M} -> G up to automorphisms.  For each structure we
compute, with flatlin.FlatProblem:

  * whether all irreducibles are admissible (condition B),
  * whether the diagonal kernel is abelian (condition A, from GAP),
  * the number of flat covers,
  * flat linearity (kernel criterion),
  * stationarity (exact LP, or by the necessary condition B),

and we check the theorems and test the questions:

  Theorem:  A => stationary,  stationary => B.
  Question A: stationary <=> A ?
  Question B: stationary <=> B ?
  Question F: is every finite uniform group flat linear ?

Usage:  python3 census.py manul.log ali.log.gz ...   -> census.txt
(logs may be gzipped, as census/ali.log.gz, which is read as is)
"""

import ast
import gzip
import sys
import time

from flatlin import FlatProblem


def parse(path):
    recs, skips = [], []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith("SKIP|"):
                skips.append(line.split("|")[1:])
                continue
            if not line.startswith("U|"):
                continue
            p = line.split("|")
            n, i, e, K, M, sym, lam = (int(x) for x in p[1:8])
            dims = ast.literal_eval(p[8])
            eig = ast.literal_eval(p[9])
            kers = [[c - 1 for c in k] for k in ast.literal_eval(p[10])]
            sizes = ast.literal_eval(p[11])
            recs.append(dict(n=n, i=i, e=e, K=K, M=M, sym=sym, lam=lam,
                             dims=dims, eig=eig, kers=kers, sizes=sizes,
                             source=path))
    return recs, skips


def aggregated_lp(P):
    """Conjecture C test: the LP that only remembers dimensions.  Variables
    are the dimension types of the flat covers, and the equations say that
    the dimension-d covers carry the mass K d n_d / |Gamma|."""
    import numpy as np
    from scipy.optimize import linprog
    covers = P.flat_covers()
    if not P.covers_complete or not covers:
        return None
    types = sorted(set(tuple(sorted(int(P.dims[t]) for t in S))
                       for S in covers))
    ds = sorted(set(int(d) for d in P.dims))
    A = np.array([[T.count(d) for T in types] for d in ds], float)
    b = np.array([sum(P.K * d / P.order for x in P.dims if x == d)
                  for d in ds])
    res = linprog(np.zeros(len(types)), A_eq=A, b_eq=b, bounds=(0, None),
                  method="highs")
    return int(res.status == 0)


def analyze(rec):
    P = FlatProblem(rec["dims"], rec["eig"], rec["kers"], rec["sizes"],
                    rec["K"], rec["M"])
    out = dict(rec)
    out["B"] = int(P.all_admissible())
    out["ncovers"] = P.ncovers()
    out["flatlin"] = int(P.is_flat_linear())
    out["kerorder"] = P.flat_kernel_order()
    if not out["B"]:
        out["stat"] = 0
        out["how"] = "not B"
    else:
        verdict, _ = P.stationarity()
        out["stat"] = {"stationary": 1, "not stationary": 0}.get(verdict, -1)
        out["how"] = "LP" if verdict != "undecided" else "undecided"
    out["C"] = aggregated_lp(P) if out["B"] else 0
    if out["stat"] == -1 and out["lam"]:
        out["stat"] = 1
        out["how"] = "theorem (A)"
    degs = {}
    for d in rec["dims"]:
        degs[d] = degs.get(d, 0) + 1
    out["degs"] = " ".join("%d^%d" % (d, m) for d, m in sorted(degs.items()))
    return out


def main(paths):
    t0 = time.time()
    recs, skips = [], []
    for p in paths:
        r, s = parse(p)
        recs += r
        skips += s
    results = []
    for r in recs:
        results.append(analyze(r))
    results.sort(key=lambda r: (r["K"], r["M"], r["n"], r["i"], r["e"]))
    lines = []
    lines.append("census of finite uniform groups, %d structures from %s"
                 % (len(results), ", ".join(paths)))
    lines.append("columns: K M | SmallGroup(n,i) epi e | sym | A=diag.kernel "
                 "abelian | B=all admissible | covers | flat linear | "
                 "stationary (how) | degrees")
    lines.append("")
    for r in results:
        lines.append("K=%d M=%d | (%d,%d) e=%d | sym=%d | A=%d | B=%d | "
                     "covers=%s | flatlin=%d | stat=%d (%s) | C=%s | %s"
                     % (r["K"], r["M"], r["n"], r["i"], r["e"], r["sym"],
                        r["lam"], r["B"], r["ncovers"], r["flatlin"],
                        r["stat"], r["how"], r["C"], r["degs"]))
    lines.append("")
    # checks of the theorems
    bad = [r for r in results if r["lam"] and r["stat"] == 0]
    lines.append("THEOREM CHECK  A => stationary: %d violations" % len(bad))
    bad2 = [r for r in results if r["stat"] == 1 and not r["B"]]
    lines.append("THEOREM CHECK  stationary => B: %d violations" % len(bad2))
    # questions, restricted to symmetric structures
    sym = [r for r in results if r["sym"]]
    lines.append("symmetric (uniform) structures: %d of %d"
                 % (len(sym), len(results)))
    qa = [r for r in sym if r["stat"] == 1 and not r["lam"]]
    lines.append("QUESTION A  (stationary but kernel non-abelian): %d" % len(qa))
    for r in qa:
        lines.append("   (%d,%d) e=%d K=%d M=%d  %s"
                     % (r["n"], r["i"], r["e"], r["K"], r["M"], r["degs"]))
    qb = [r for r in sym if r["B"] and r["stat"] == 0]
    lines.append("QUESTION B  (all admissible but not stationary): %d" % len(qb))
    for r in qb:
        lines.append("   (%d,%d) e=%d K=%d M=%d  %s"
                     % (r["n"], r["i"], r["e"], r["K"], r["M"], r["degs"]))
    qc = [r for r in sym if r["B"] and r["C"] is not None
          and r["C"] != (r["stat"] == 1)]
    nc = len([r for r in sym if r["B"] and r["C"] is not None])
    lines.append("CONJECTURE C (stationary <=> admissible + dimension LP): "
                 "%d disagreements out of %d admissible structures"
                 % (len(qc), nc))
    for r in qc:
        lines.append("   (%d,%d) e=%d K=%d M=%d  stat=%d C=%d  %s"
                     % (r["n"], r["i"], r["e"], r["K"], r["M"], r["stat"],
                        r["C"], r["degs"]))
    qf = [r for r in sym if not r["flatlin"]]
    lines.append("QUESTION F  (not flat linear): %d" % len(qf))
    for r in qf:
        lines.append("   (%d,%d) e=%d K=%d M=%d  kernel order %d  %s"
                     % (r["n"], r["i"], r["e"], r["K"], r["M"],
                        r["kerorder"], r["degs"]))
    und = [r for r in results if r["stat"] == -1]
    lines.append("undecided stationarity: %d" % len(und))
    for r in und:
        lines.append("   (%d,%d) e=%d K=%d M=%d  %s"
                     % (r["n"], r["i"], r["e"], r["K"], r["M"], r["degs"]))
    if skips:
        lines.append("skipped orders (too many groups): %s"
                     % sorted(set((int(s[0]), int(s[1])) for s in skips)))
    lines.append("time: %.1f s" % (time.time() - t0))
    with open("census.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines[-12:]))


if __name__ == "__main__":
    main(sys.argv[1:])
