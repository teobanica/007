"""
heights.py -- flatness at height m.

A representation rho : Gamma -> U_{mK} is flat at height m when each
rho(g_i) has all the K-th roots of unity as eigenvalues, each with
multiplicity m (the rank landscape is flat, at height m).  For a finite
group these are, up to conjugation, the sums sum_tau x_tau tau with
x in N^{Irr} and sum_tau x_tau e_tau = m (all-m matrix): the m-covers.

  * flat linear at height m: the irreducibles occurring in m-covers have
    trivial common kernel,
  * stationary at height m: there are weights c_x >= 0 with
        sum_x c_x x_tau = m K dim(tau) / |Gamma|   for all tau.

The regular representation is flat at height |Gamma|/K, faithful and
stationary, so the minimal heights h_lin <= h_stat <= |Gamma|/K exist.

Usage:  python3 heights.py census.log [maxheight]   -> heights.log
(census.log = a GAP census log; only the structures that are not flat
linear or not stationary at height 1 are examined)
"""

import sys
import time
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog

from flatlin import FlatProblem, exact_solve


def m_covers(P, m, cap=200000):
    """All x in N^{Irr} with sum_tau x_tau e_tau = m * ones."""
    width = P.M * P.K
    allowed = [r for r in range(P.nirr)
               if P.dims[r] <= m * P.K and P.eig[r].max() <= m]
    allowed.sort(key=lambda r: -P.dims[r])
    E = {r: P.eig[r] for r in allowed}
    out = []
    x = {}

    def rec(k, rem):
        if len(out) > cap:
            return
        if not rem.any():
            out.append(dict(x))
            return
        if k == len(allowed):
            return
        r = allowed[k]
        e = E[r]
        # maximal multiplicity of r
        nz = e > 0
        tmax = int(min(rem[nz] // e[nz])) if nz.any() else 0
        for t in range(tmax, -1, -1):
            if t:
                x[r] = t
            elif r in x:
                del x[r]
            rec(k + 1, rem - t * e)
        x.pop(r, None)

    rec(0, np.full(width, m, dtype=int))
    complete = len(out) <= cap
    return out[:cap], complete


def height_data(P, m):
    covers, complete = m_covers(P, m)
    used = set()
    for x in covers:
        used.update(x)
    common = set(range(len(P.sizes)))
    for r in used:
        common &= P.kers[r]
    kerorder = int(sum(P.sizes[c] for c in common))
    stat = None
    if complete and covers:
        if len(used) < P.nirr:
            stat = False
        else:
            n = P.nirr
            A = np.zeros((n, len(covers)))
            for s, x in enumerate(covers):
                for r, t in x.items():
                    A[r, s] = t
            lam = [Fraction(m * P.K * int(d), P.order) for d in P.dims]
            b = np.array([float(v) for v in lam])
            res = linprog(np.zeros(len(covers)), A_eq=A, b_eq=b,
                          bounds=(0, None), method="highs-ds")
            if res.status == 0:
                sup = [s for s in range(len(covers)) if res.x[s] > 1e-9]
                sol = exact_solve([[Fraction(int(A[r, s])) for s in sup]
                                   for r in range(n)], lam)
                stat = sol is not None and all(v >= 0 for v in sol)
                if not stat:
                    stat = None          # could not certify
            else:
                res2 = linprog(-b, A_ub=A.T, b_ub=np.zeros(len(covers)),
                               bounds=(-1, 1), method="highs-ds")
                y = [Fraction(v).limit_denominator(10 ** 6) for v in res2.x]
                ok = all(sum(y[r] * t for r, t in x.items()) <= 0
                         for x in covers)
                val = sum(yi * li for yi, li in zip(y, lam))
                stat = False if (ok and val > 0) else None
    elif complete and not covers:
        stat = False
    return dict(m=m, ncovers=len(covers) if complete else None,
                kerorder=kerorder, flatlin=(kerorder == 1), stat=stat)


def heights(P, maxm=4):
    """Minimal heights for flat linearity and stationarity, up to maxm."""
    hlin = hstat = None
    rows = []
    for m in range(1, maxm + 1):
        d = height_data(P, m)
        rows.append(d)
        if hlin is None and d["flatlin"]:
            hlin = m
        if hstat is None and d["stat"]:
            hstat = m
        if hlin is not None and hstat is not None:
            break
    return hlin, hstat, rows


def main(path, maxm=4):
    from census import parse
    recs, _ = parse(path)
    out = open("heights.log", "w")
    t0 = time.time()
    for rec in recs:
        if not rec["sym"]:
            continue
        P = FlatProblem(rec["dims"], rec["eig"], rec["kers"], rec["sizes"],
                        rec["K"], rec["M"])
        d1 = height_data(P, 1)
        if d1["flatlin"] and d1["stat"]:
            continue
        hlin, hstat, rows = heights(P, maxm)
        degs = {}
        for dd in rec["dims"]:
            degs[dd] = degs.get(dd, 0) + 1
        line = ("(%d,%d) e=%d K=%d M=%d | %s | h_lin=%s h_stat=%s | %s"
                % (rec["n"], rec["i"], rec["e"], rec["K"], rec["M"],
                   " ".join("%d^%d" % t for t in sorted(degs.items())),
                   hlin if hlin else ">%d" % maxm,
                   hstat if hstat else ">%d" % maxm,
                   "; ".join("m=%d: covers=%s ker=%d stat=%s"
                             % (r["m"], r["ncovers"], r["kerorder"], r["stat"])
                             for r in rows)))
        print(line, flush=True)
        out.write(line + "\n")
        out.flush()
    out.write("time %.0fs\n" % (time.time() - t0))
    out.close()


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 4)
