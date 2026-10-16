"""Peeling DE thresholds vs k, converted to communication cost / lower bound (n/d = 101).
Ideal cost uses per-row entropy H(Poisson(k d / l)); counting bound is the l/d where it equals LB."""
import math
from peeling_de import uniform_converges, coupled_converges, bisect

def Hpois(lam):
    H, lp = 0.0, -lam
    for x in range(0, 500):
        if x > 0: lp += math.log(lam) - math.log(x)
        p = math.exp(lp)
        if p > 0: H -= p * lp / math.log(2)
        if x > lam and p < 1e-18: break
    return H

nd, n, d = 101, 1010000, 10000
LB = (math.lgamma(n + 1) - math.lgamma(d + 1) - math.lgamma(n - d + 1)) / math.log(2) / d
print(f"LB/d = {LB:.3f} bits per element of B\\A")
print(f"{'k':>3} {'count.bnd':>9} {'unif.peel':>9} {'H/LB':>6} {'SC.peel':>8} {'H/LB':>6}")
for k in (4, 6, 8, 10, 12, 16, 20, 24):
    cb = bisect(lambda r: r * Hpois(k / r) >= LB, 0.5, 30)
    ru = bisect(lambda r: uniform_converges(r, k, nd), 1.0, 30.0)
    rc = bisect(lambda r: coupled_converges(r, k, nd, 300, 15)[0], 1.0, ru)
    print(f"{k:>3} {cb:9.2f} {ru:9.2f} {ru*Hpois(k/ru)/LB:6.2f} {rc:8.2f} {rc*Hpois(k/rc)/LB:6.2f}")
