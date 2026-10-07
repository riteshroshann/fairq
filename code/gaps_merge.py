"""Merge chunks from results/parts/ into results/gaps_<family>.json and fit exponents."""
import glob, json, sys
from gaps import fit, SIZES

kind = sys.argv[1]
rows = []
for f in sorted(glob.glob(f"results/parts/gaps_{kind}_N*_s*.json")):
    rows += json.load(open(f))
rows.sort(key=lambda r: (r["N"], r["seed"], r["beta"]))
expected = len(SIZES[kind]) * 5 * 2
assert len(rows) == expected, (len(rows), expected)
json.dump(rows, open(f"results/gaps_{kind}.json", "w"), indent=1)
json.dump(fit(rows), open(f"results/exponents_{kind}.json", "w"), indent=1)
print(json.dumps(fit(rows), indent=1))
