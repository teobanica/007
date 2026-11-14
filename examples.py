"""
examples.py -- the basic families of finite uniform groups.

Runs the flat linearity and stationarity tests of flatlin.py on:
  * the abelian groups Z_K^M,
  * the dihedral groups D_n, n even (K = 2, M = 2),
  * the Heisenberg groups H_K over Z_K, of order K^3 (M = 2),
  * the metabelian groups (Z_3*Z_3)_met / qA, q = 2, 5, of orders 144, 5625.
Writes everything to examples.log.
"""

import itertools
import sys

import numpy as np

from flatlin import UniformGroup, perm_from_function


def abelian(K, M):
    pts = list(itertools.product(range(K), repeat=M))
    gens = []
    for i in range(M):
        def f(p, i=i):
            q = list(p)
            q[i] = (q[i] + 1) % K
            return tuple(q)
        gens.append(perm_from_function(f, pts))
    return UniformGroup(gens, K, "Z_%d^%d" % (K, M))


def dihedral(n):
    pts = list(range(n))
    s1 = perm_from_function(lambda i: (-i) % n, pts)
    s2 = perm_from_function(lambda i: (1 - i) % n, pts)
    return UniformGroup([s1, s2], 2, "D_%d (order %d)" % (n, 2 * n))


def heisenberg(p):
    pts = list(itertools.product(range(p), repeat=3))
    x = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    y = np.array([[1, 0, 0], [0, 1, 1], [0, 0, 1]])
    def act(A):
        return lambda v: tuple(int(t) % p for t in A.dot(v))
    return UniformGroup([perm_from_function(act(x), pts),
                         perm_from_function(act(y), pts)],
                        p, "Heisenberg H_%d (order %d)" % (p, p ** 3))


def companion(K, q):
    """Matrix of t on Z_q[t]/(1+t+...+t^{K-1}), size K-1."""
    m = np.zeros((K - 1, K - 1), dtype=int)
    for i in range(1, K - 1):
        m[i, i - 1] = 1
    m[:, K - 2] = -1
    return m % q


def metabelian_quotient(K, q, a0=None):
    """(Z_K * Z_K)_met / q A, realized as an affine group on A = Z_q^{(K-1)^2}
    with x : v -> (m (x) 1) v and y : v -> (1 (x) m) v + a0."""
    m = companion(K, q)
    I = np.eye(K - 1, dtype=int)
    X = np.kron(m, I) % q
    Y = np.kron(I, m) % q
    d = (K - 1) ** 2
    pts = list(itertools.product(range(q), repeat=d))
    if a0 is None:
        a0 = np.zeros(d, dtype=int)
        a0[0] = 1
    def fx(v):
        return tuple(int(t) % q for t in X.dot(v))
    def fy(v):
        return tuple(int(t) % q for t in (Y.dot(v) + a0))
    return UniformGroup([perm_from_function(fx, pts),
                         perm_from_function(fy, pts)],
                        K, "(Z_%d*Z_%d)_met / %dA" % (K, K, q))


def main():
    out = open("examples.log", "w")
    def emit(s):
        print(s)
        out.write(s + "\n")
        out.flush()
    tests = []
    for K, M in [(2, 2), (2, 3), (3, 2), (2, 4), (3, 3), (4, 2), (5, 2)]:
        tests.append(lambda K=K, M=M: abelian(K, M))
    for n in [4, 6, 8, 10, 12, 16, 20, 24]:
        tests.append(lambda n=n: dihedral(n))
    for p in [2, 3, 4, 5, 6, 7, 8, 9]:
        tests.append(lambda p=p: heisenberg(p))
    tests.append(lambda: metabelian_quotient(3, 2))
    tests.append(lambda: metabelian_quotient(3, 5))
    for t in tests:
        U = t()
        emit(U.report())
    out.close()


if __name__ == "__main__":
    main()
