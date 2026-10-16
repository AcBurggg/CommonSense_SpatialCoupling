"""Density evolution (DE) for exact peeling on CommonSense's counting sketch (A subset of B).

Peeling rules on a row with residual r and u undecided candidates:
  r == 0  -> all undecided candidates are 0   (resolves false candidates)
  r == u  -> all undecided candidates are 1   (resolves true candidates)
Let p_t (q_t) be the probability that an edge message from a true (false) candidate is still
undecided. With Poisson row degrees, lam_T = k d / l true and lam_F = k (n - d) / l false edges
per row, a false candidate is resolved through a row iff the row has no other undecided true
element, and a true one iff the row has no other undecided false element:
  q_{t+1} = (1 - exp(-lam_T p_t))^(k-1),    p_{t+1} = (1 - exp(-lam_F q_t))^(k-1).
Spatially coupled version (Walzer-style, positions 0..L-1, window w positions): a row at position j
sees the window average of p, q over variable positions (j-w, j], with zero padding outside, and a
variable at position i averages the row-resolution probability over rows in [i, i+w).

Usage: python3 peeling_de.py
"""
import numpy as np


def uniform_converges(r, k, nd, iters=20000):
    lam_t, lam_f = k / r, k * (nd - 1) / r
    p = q = 1.0
    for _ in range(iters):
        p_new = (1 - np.exp(-lam_f * q)) ** (k - 1)
        q_new = (1 - np.exp(-lam_t * p)) ** (k - 1)
        if p_new < 1e-12 and q_new < 1e-12:
            return True
        if abs(p_new - p) < 1e-14 and abs(q_new - q) < 1e-14:
            return False
        p, q = p_new, q_new
    return False


def coupled_converges(r_bulk, k, nd, L, w, iters=200000):
    """r_bulk = bulk rows per true element (l_bulk / d): rows per position are uniform."""
    lam_t, lam_f = k / r_bulk, k * (nd - 1) / r_bulk
    p = np.ones(L)
    q = np.ones(L)
    kernel = np.ones(w) / w
    for it in range(iters):
        # row positions 0..L+w-2; window average of variables in (j-w, j]
        pbar = np.convolve(p, kernel)  # length L+w-1
        qbar = np.convolve(q, kernel)
        res_f = np.exp(-lam_t * pbar)  # row resolves a false candidate
        res_t = np.exp(-lam_f * qbar)  # row resolves a true candidate
        # variable at i averages rows in [i, i+w)
        avg_res_f = np.convolve(res_f, kernel, mode="valid")  # length L
        avg_res_t = np.convolve(res_t, kernel, mode="valid")
        q_new = (1 - avg_res_f) ** (k - 1)
        p_new = (1 - avg_res_t) ** (k - 1)
        if p_new.max() < 1e-10 and q_new.max() < 1e-10:
            return True, it
        if np.abs(p_new - p).max() < 1e-13 and np.abs(q_new - q).max() < 1e-13:
            return False, it
        p, q = p_new, q_new
    return False, iters


def bisect(pred, lo, hi, tol=0.005):
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return hi


if __name__ == "__main__":
    nd = 101  # |B| / |B \ A| in the paper's first unidirectional group
    for k in (6, 8, 10):
        ru = bisect(lambda r: uniform_converges(r, k, nd), 1.0, 20.0)
        print(f"k={k}: uniform peeling threshold l/d = {ru:.3f}")
        for L, w in ((200, 10), (400, 20), (800, 40)):
            rb = bisect(lambda r: coupled_converges(r, k, nd, L, w)[0], 1.0, ru)
            # total rows include the w-1 extra positions of the coupled chain
            total = rb * (L + w - 1) / L
            print(f"   coupled L={L:4d} w={w:3d}: bulk l/d = {rb:.3f}   total l/d incl. boundary = {total:.3f}")
