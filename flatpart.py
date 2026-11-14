"""
flatpart.py -- the flat part of the reflection group of H_N^{[4]}.

All flat representations of Gamma_N^{[4]} send the generators to sigma_p with
p among the 6 vertices +-e_1, +-e_2, +-e_3 of an orthonormal frame (up to a
global rotation).  The flat part F_N = Gamma_N^{[4]} / N_flat is the image of
the direct sum of all these representations.  We compute |F_N|, and test
which relations hold in F_N.  Writes flatpart.log.
"""
import itertools
import numpy as np

S = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex),
     np.array([[1, 0], [0, -1]], complex)]
VERT = [s * S[k] for k in range(3) for s in (1, -1)]


def image_group(N, cap=10 ** 6):
    assigns = list(itertools.product(range(6), repeat=N))
    gens = [np.array([VERT[a[i]] for a in assigns]) for i in range(N)]  # (A,2,2)
    I = np.array([np.eye(2, dtype=complex)] * len(assigns))
    def key(X):
        return np.round(X, 6).tobytes()
    seen = {key(I): I}
    frontier = [I]
    while frontier:
        new = []
        for X in frontier:
            for g in gens:
                Y = np.einsum('aij,ajk->aik', g, X)
                k = key(Y)
                if k not in seen:
                    seen[k] = Y
                    new.append(Y)
                    if len(seen) > cap:
                        return None, gens
        frontier = new
    return len(seen), gens


def holds(word_fn, gens, N):
    """word_fn(list of generator arrays) -> array; check equals identity."""
    W = word_fn(gens)
    I = np.eye(2)
    return np.allclose(W, I[None, :, :])


def mul(*Xs):
    R = Xs[0]
    for X in Xs[1:]:
        R = np.einsum('aij,ajk->aik', R, X)
    return R


def main():
    out = open("flatpart.log", "w")
    def emit(s):
        print(s); out.write(s + "\n"); out.flush()
    for N in [2, 3, 4]:
        order, g = image_group(N)
        emit("N=%d: |F_N| = %s" % (N, order))
        if N >= 3:
            a, b, c = g[0], g[1], g[2]
            tests = {
                "(ab)^4 = 1": mul(a, b, a, b, a, b, a, b),
                "(ab)^2 central: [(ab)^2, c] = 1":
                    mul(a, b, a, b, c, b, a, b, a, c),
                "abc = cba  (half-liberation)": mul(a, b, c, a, b, c),
                "(abc)^2 = 1": mul(a, b, c, a, b, c),
                "(abc)^4 = 1": mul(*([a, b, c] * 4)),
                "(ab)^2 = (ba)^2": mul(a, b, a, b, a, b, a, b),
                "abab = baba": mul(a, b, a, b, a, b, a, b),
                "[ab, ac]^2 = 1": None,
            }
            for name, W in tests.items():
                if W is None:
                    ab = mul(a, b); ac = mul(a, c)
                    abi = mul(b, a); aci = mul(c, a)
                    comm = mul(ab, ac, abi, aci)
                    W = mul(comm, comm)
                emit("   %-40s %s" % (name, np.allclose(W, np.eye(2)[None])))
    out.close()


if __name__ == "__main__":
    main()
