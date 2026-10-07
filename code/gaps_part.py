"""Run one chunk of the E2 sweep and store it in results/parts/ (used to stay within time limits).
Usage: python gaps_part.py <family> <N> <seed> [<seed> ...];  then python gaps_merge.py <family>"""
import json, sys
from gaps import sweep

kind, N, seeds = sys.argv[1], int(sys.argv[2]), [int(s) for s in sys.argv[3:]]
rows = sweep(kind, (N,), seeds=tuple(seeds))
for sd in seeds:
    json.dump([r for r in rows if r["seed"] == sd], open(f"results/parts/gaps_{kind}_N{N}_s{sd}.json", "w"), indent=1)
