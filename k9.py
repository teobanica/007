"""
k9.py -- the 4 uniform structures at K=9, M=2 with more than 400000 flat covers
(Proposition 8.6).  Reads liston.log and checks the spectral facts used in the
proof:

  * all irreducibles are admissible;
  * on each generator, the spectrum of each irreducible of dimension 3 is a
    residue class modulo 3 inside {1,...,9};
  * on the first generator, the spectrum of each irreducible of dimension 2
    consists of 2 elements of a residue class modulo 3;
  * the 9 pairs of residue classes (generator 1, generator 2) occur equally
    often among the irreducibles of dimension 3.

Usage:  python3 k9.py census/liston.log
"""
import sys
from collections import Counter
from census import parse

CASES = {(648, 90, 1), (648, 90, 2), (729, 26, 4), (972, 173, 2)}

recs, _ = parse(sys.argv[1] if len(sys.argv) > 1 else "census/liston.log")
for r in recs:
    key = (r["n"], r["i"], r["e"])
    if key not in CASES:
        continue
    K = r["K"]
    admissible = True
    structure = True
    cells = Counter()
    for d, row in zip(r["dims"], r["eig"]):
        rows = [row[i * K:(i + 1) * K] for i in range(2)]
        if any(max(rw) > 1 for rw in rows):
            admissible = False
        slots = [[j + 1 for j in range(K) if rw[j]] for rw in rows]
        classes = [sorted(set(s % 3 for s in sl)) for sl in slots]
        if d == 3:
            if all(len(sl) == 3 and len(c) == 1 for sl, c in zip(slots, classes)):
                cells[(classes[0][0], classes[1][0])] += 1
            else:
                structure = False
        if d == 2 and not (len(slots[0]) == 2 and len(classes[0]) == 1):
            structure = False
    dims = Counter(r["dims"])
    print("SmallGroup%s e=%d  degrees %s  admissible=%s  spectral structure=%s"
          % ((r["n"], r["i"]), r["e"], dict(sorted(dims.items())), admissible, structure))
    print("   pairs of residue classes for dimension 3:", dict(sorted(cells.items())),
          " equal:", len(set(cells.values())) == 1 and len(cells) == 9)
