"""Summarize sc_sim result logs (results/*.jsonl) into per-configuration threshold tables.

For each (scale, k, matrix, window, offsets, decoder) group we report every measured l/d with its
success rate, and the smallest l/d at which all trials succeeded ("threshold"), together with the
communication cost there relative to the SetX lower bound log2 C(|B|, |B \\ A|):
  H/LB    ideal per-row Slepian-Wolf cost  sum_j H(Poisson(lambda_j))
  CSg/LB  CommonSense's mod-q + BCH scheme with one parameter set (tuned for the bulk rows)
  CSa/LB  same scheme with parameters per row-load class (cheaper coupling boundary rows)
Usage: python3 analyze.py [results/*.jsonl]
"""
import glob
import json
import sys
from collections import defaultdict


def load(paths):
    rows = []
    for p in paths:
        with open(p) as f:
            for line in f:
                line = line.strip()
                if line.startswith("{"):
                    rows.append(json.loads(line))
    return rows


def main():
    paths = sys.argv[1:] or sorted(glob.glob("results/*.jsonl"))
    rows = load(paths)
    groups = defaultdict(list)
    for r in rows:
        key = (r["n"], r["d"], r["k"], r["matrix"], round(r["wfrac"], 3), r["offsets"], r["decoder"])
        groups[key].append(r)
    total_trials = sum(r["trials"] for r in rows)
    total_silent = sum(r["silent_errors"] for r in rows)
    total_ok = sum(r["success"] for r in rows)
    print(f"{len(rows)} configurations, {total_trials} trials, {total_ok} successes, "
          f"{total_silent} silent errors (zero residual but wrong set)\n")
    hdr = f"{'|B|':>8} {'d':>6} {'k':>2} {'matrix':7} {'W/l':>5} {'offsets':7} {'decoder':5} | " \
          f"{'measured l/d: successes/trials':45} | {'thresh':>6} {'H/LB':>5} {'CSg/LB':>6} {'CSa/LB':>6}"
    print(hdr)
    print("-" * len(hdr))
    for key in sorted(groups, key=lambda k: (k[0], k[1], k[2], k[3] != "uniform", k[6], k[4], k[5])):
        g = sorted(groups[key], key=lambda r: r["ratio"])
        cells = " ".join(f"{r['ratio']:.2f}:{r['success']}/{r['trials']}" for r in g)
        full = [r for r in g if r["success"] == r["trials"]]
        if full:
            t = full[0]
            tail = f"{t['ratio']:6.2f} {t['entropy_over_lb']:5.2f} {t['cs_global_over_lb']:6.2f} {t['cs_adaptive_over_lb']:6.2f}"
        else:
            tail = f"{'>max':>6}"
        n, d, k, matrix, wf, off, dec = key
        print(f"{n:>8} {d:>6} {k:>2} {matrix:7} {wf:5.2f} {off:7} {dec:5} | {cells:45} | {tail}")


if __name__ == "__main__":
    main()
