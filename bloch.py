"""
bloch.py -- flat representations of the reflection groups
Gamma_N^{[s]} = Z_2^{*N} / <(g_i g_j)^s>  (Raum-Weber's H_N^{[s]} series).

For K = 2, a flat representation sends g_i to p_i . sigma (Pauli matrices),
with p_i on the 2-sphere, and (g_i g_j)^s = 1 iff the angle between p_i, p_j
is a multiple of 2 pi / s.  For N = 3 we list the realizable angle patterns,
build the 3 matrices, and compute the order of the generated subgroup of U_2
(capped).  If all images are finite, the flat representations factor through
finite groups, while Gamma_3^{[s]} (the (s,s,s) Coxeter group) is infinite
for s >= 4: then Gamma is not flat linear.
Writes bloch.log.
"""
import itertools
import numpy as np

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def refl(p):
    return p[0] * SX + p[1] * SY + p[2] * SZ


def key(A):
    return tuple(np.round(A.ravel(), 6))


def group_order(gens, cap=20000):
    I = np.eye(2, dtype=complex)
    seen = {key(I): I}
    frontier = [I]
    while frontier:
        new = []
        for A in frontier:
            for g in gens:
                B = g @ A
                k = key(B)
                if k not in seen:
                    seen[k] = B
                    new.append(B)
                    if len(seen) > cap:
                        return None
        frontier = new
    return len(seen)


def configs(s):
    angles = [2 * np.pi * k / s for k in range(s // 2 + 1)]
    out = []
    for a12, a13, a23 in itertools.product(angles, repeat=3):
        p1 = np.array([0, 0, 1.0])
        p2 = np.array([np.sin(a12), 0, np.cos(a12)])
        # p3 = (x, y, z) with p3.p1 = cos a13, p3.p2 = cos a23
        z = np.cos(a13)
        if abs(np.sin(a12)) < 1e-12:
            if abs(np.cos(a23) - np.cos(a12) * z) > 1e-9:
                continue
            x = np.sqrt(max(0, 1 - z * z)); y = 0.0
        else:
            x = (np.cos(a23) - np.cos(a12) * z) / np.sin(a12)
            r = 1 - z * z - x * x
            if r < -1e-9:
                continue
            y = np.sqrt(max(r, 0))
        p3 = np.array([x, y, z])
        out.append(((a12, a13, a23), [p1, p2, p3]))
    return out


def main():
    f = open("bloch.log", "w")
    for s in [4, 6, 8, 10, 12]:
        res = []
        for ang, ps in configs(s):
            o = group_order([refl(p) for p in ps])
            res.append((tuple(int(round(a * s / (2 * np.pi))) for a in ang), o))
        orders = sorted(set(o for _, o in res if o is not None))
        inf = [a for a, o in res if o is None]
        line = ("s=%d: %d angle patterns, finite image orders %s, "
                "infinite images for patterns (in units 2pi/s): %s"
                % (s, len(res), orders, inf))
        print(line)
        f.write(line + "\n")
    f.close()


if __name__ == "__main__":
    main()
