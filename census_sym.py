"""
census_sym.py -- faster variant of census.py, used for liston.log and canelo.log.

Same analysis as census.py (flatlin.FlatProblem: admissibility, flat covers,
flat linearity, exact stationarity LP, dimension-type LP), but restricted to
the symmetric epimorphisms, i.e. to the uniform structures.  Writes one JSON
line per uniform structure.

Usage:  python3 census_sym.py liston.log liston.json
"""
import sys, json
from census import parse, analyze

path, outpath = sys.argv[1], sys.argv[2]
out = open(outpath, "w")
recs, skips = parse(path)
for r in recs:
    if not r["sym"]:
        continue
    o = analyze(r)
    keys = ["n", "i", "e", "K", "M", "sym", "lam", "B", "ncovers", "flatlin",
            "kerorder", "stat", "how", "C", "degs"]
    out.write(json.dumps({k: o.get(k) for k in keys}) + "\n")
    out.flush()
