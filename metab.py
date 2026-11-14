"""metab.py -- the groups Gamma_{K,q} = (Z_K*Z_K)_met / qA, for any q (including q | K),
built from the normal form c^v g1^x g2^y, where c_{ab} = [g1^a, g2^b], 1<=a,b<=K-1,
is the free basis of the derived subgroup of Z_K*Z_K. Writes metab.log."""
import itertools, sys
import numpy as np
from flatlin import UniformGroup, perm_from_function

def build(K, q):
    d = (K - 1) ** 2
    idx = {(a, b): (a - 1) * (K - 1) + (b - 1) for a in range(1, K) for b in range(1, K)}
    def e(a, b):
        v = np.zeros(d, dtype=np.int64)
        a %= K; b %= K
        if a and b:
            v[idx[(a, b)]] = 1
        return v
    T1 = np.zeros((d, d), dtype=np.int64); T2 = np.zeros((d, d), dtype=np.int64)
    for (a, b), j in idx.items():
        T1[:, j] = e(a + 1, b) - e(1, b)     # g1 c_ab g1^-1 = c_{a+1,b} c_{1,b}^-1
        T2[:, j] = e(a, b + 1) - e(a, 1)     # g2 c_ab g2^-1 = c_{a,1}^-1 c_{a,b+1}
    P1 = [np.eye(d, dtype=np.int64)]; P2 = [np.eye(d, dtype=np.int64)]
    for _ in range(K - 1):
        P1.append(T1.dot(P1[-1]) % q); P2.append(T2.dot(P2[-1]) % q)
    T = lambda x, y: P1[x % K].dot(P2[y % K]) % q
    def mul(s, t):
        v, x, y = s; w_, x2, y2 = t
        w = -e(-x2, -y)                      # g2^y g1^x2 = g1^x2 g2^y [g2^-y, g1^-x2]
        vec = (np.array(v) + T(x, y).dot(w_) + T(x, 0).dot(T(x2, y)).dot(w)) % q
        return (tuple(int(z) for z in vec), (x + x2) % K, (y + y2) % K)
    pts = [(v, x, y) for v in itertools.product(range(q), repeat=d) for x in range(K) for y in range(K)]
    g1 = ((0,) * d, 1, 0); g2 = ((0,) * d, 0, 1)
    # sanity: associativity on random triples
    rng = np.random.default_rng(1)
    for _ in range(200):
        a, b, c = (pts[i] for i in rng.integers(len(pts), size=3))
        assert mul(mul(a, b), c) == mul(a, mul(b, c))
    return UniformGroup([perm_from_function(lambda p: mul(g1, p), pts),
                         perm_from_function(lambda p: mul(g2, p), pts)],
                        K, "Gamma_{%d,%d}" % (K, q))

if __name__ == "__main__":
    out = open("metab.log", "a")
    for arg in sys.argv[1:]:
        K, q = map(int, arg.split(","))
        r = build(K, q).report()
        print(r); out.write(r + "\n"); out.flush()
