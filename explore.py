"""
explore.py -- random search for small uniform groups among permutation groups.

Picks random M-tuples of elements of order K in S_n (as products of disjoint
K-cycles), keeps the tuples generating a group Gamma with Gamma/Gamma' of
order K^M (so Gamma^ab = Z_K^M) and with the S_M symmetry, and runs the
flat linearity / stationarity analysis.  Groups are identified up to the
invariant (order, degrees, number of flat covers, diagonal kernel abelian).
Results go to explore.log.  Usage: python3 explore.py K M n tries seed
"""

import random
import sys
import time

from flatlin import UniformGroup, FiniteGroup


def random_order_K(n, K, rng):
    pts = list(range(n))
    rng.shuffle(pts)
    ncyc = rng.randint(1, n // K)
    p = list(range(n))
    for c in range(ncyc):
        cyc = pts[c * K:(c + 1) * K]
        for t in range(K):
            p[cyc[t]] = cyc[(t + 1) % K]
    return tuple(p)


def derived_order(G):
    gens = G.gens
    seeds = set()
    for a in gens:
        for b in gens:
            seeds.add(G.mul(G.mul(a, b), G.mul(int(G.inv[a]), int(G.inv[b]))))
    conj = set()
    for s in seeds:
        for x in range(G.order):
            conj.add(G.mul(G.mul(x, s), int(G.inv[x])))
    conj = list(conj)
    elems, frontier = {0}, [0]
    while frontier:
        new = []
        for x in frontier:
            for s in conj:
                y = G.mul(x, s)
                if y not in elems:
                    elems.add(y)
                    new.append(y)
        frontier = new
    return len(elems)


def main(K, M, n, tries, seed, maxorder=3000):
    rng = random.Random(seed)
    seen = {}
    out = open("explore-%d-%d-%d.log" % (K, M, n), "w")
    t0 = time.time()
    for t in range(tries):
        gens = [random_order_K(n, K, rng) for _ in range(M)]
        try:
            G = FiniteGroup(gens, max_order=maxorder)
        except ValueError:
            continue
        if G.order % (K ** M):
            continue
        if derived_order(G) * K ** M != G.order:
            continue
        U = UniformGroup(gens, K, "random")
        if not U.has_symmetric_action():
            continue
        P = U.problem
        lam = U.diagonal_kernel_is_abelian()
        B = P.all_admissible()
        key = (P.order, tuple(sorted(P.dims)), lam, B)
        if key in seen:
            continue
        verdict, _ = P.stationarity()
        rec = (P.order, " ".join("%d^%d" % (d, list(P.dims).count(d))
                                 for d in sorted(set(P.dims))),
               "A=%d" % lam, "B=%d" % B, "covers=%s" % P.ncovers(),
               "flatlin=%d" % P.is_flat_linear(), verdict)
        seen[key] = (rec, gens)
        line = " | ".join(str(x) for x in rec)
        flag = ""
        if verdict == "stationary" and not lam:
            flag = "   <-- QUESTION A counterexample"
        if B and verdict == "not stationary":
            flag = "   <-- QUESTION B counterexample"
        if not P.is_flat_linear():
            flag += "   <-- NOT FLAT LINEAR"
        out.write(line + flag + "\n")
        out.write("     gens = %s\n" % [list(g) for g in gens])
        out.flush()
        print(line + flag, flush=True)
    out.write("tries=%d, distinct=%d, time=%.0fs\n"
              % (tries, len(seen), time.time() - t0))
    out.close()


if __name__ == "__main__":
    K, M, n, tries, seed = (int(x) for x in sys.argv[1:6])
    main(K, M, n, tries, seed)
