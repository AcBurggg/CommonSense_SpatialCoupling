# Spatial coupling for CommonSense: applicability, correctness, implementation

Scope: §1–5 cover unidirectional SetX (A ⊆ B), which is where the paper's information-theoretic
comparisons are cleanest. §6 simulates bidirectional SetX (CommonSense's ping-pong protocol).

Notation follows the CommonSense paper: `n = |B|`, `d = |B \ A|`, `M` is `l × n` with `k` ones per
column, Bob decodes the `d`-sparse binary `x = 1_{B\A}` from `y = M1_B − M1_A = Mx`. The SetX lower
bound is `LB = log2 C(n, d)` bits (≈ 8.09 bits per element of `B \ A` when `n/d ≈ 101`).

## Summary

**(a) Applicable: yes, but only together with a message-passing decoder.** CommonSense's sketch
`M1_A` is a counting Bloom filter, and Walzer's spatial coupling transfers to it directly: hash
each element to a position and confine its `k` counters to a window. METTLE's *touch-less*
leading edge does not transfer, because SetX encoding is blind. The paper's MP decoder gets
*worse* on a coupled matrix (its threshold rises from 7.8d to 8.4d at paper scale, ≈9d at small
scale), because localized columns are much more coherent. With exact peeling or BP, coupling produces threshold saturation, as density
evolution predicts.

**Correctness is preserved and strengthened.** The paper's RIP-1 certificate does not actually
apply at its own operating point (it needs ≥ 13.3d rows; even the weakest expansion argument
needs > 8d), so CommonSense already relies on average-case uniqueness. With coupling:
* exact peeling, when it completes, is a *deterministic* proof that the decoded set is the unique
  solution, whatever the matrix;
* BP's zero-residue output is correct w.h.p. above the MAP threshold;
* a 64-bit set hash (8 bytes) makes either verifiable.

There were 0 silent errors (zero residue with a wrong set) in 1,816 logged trials.

**(b) Implemented** (opt-in, defaults unchanged): `DoroCode::enable_spatial_coupling`
(`lib/doro.hpp`), a peeling + BP decoder (`lib/sc_decoder.hpp`), and config plumbing in
`bin/two_party.cpp`.

**Results, paper scale** (`|B| = 1.01M`, `|B \ A| = 10K`, the paper's first unidirectional group):

| Scheme | rows `l/d` at 100% success | 1st-message cost / LB (CommonSense coder) | ideal coder / LB |
|---|---|---|---|
| CommonSense as published (uniform, MP, k=8) | 7.8 | **2.31** | 1.84 |
| uniform + exact peeling, k=8 | 7.6 | 2.26 | 1.81 |
| uniform + BP, k=8 (decoder upgrade only) | 5.8 | 1.86 | 1.54 |
| **coupled (5% window, TLE) + exact peeling, k=8** | **5.2** | **1.74** | 1.42 |
| coupled + exact peeling, k=6 | 5.4 | 1.64 | 1.30 |
| **coupled (5% window, TLE) + BP, k=8** | **4.4** | **1.56** | 1.28 |
| density-evolution limit, coupled peeling, k=8 (5% window) | 4.46* | 1.57* | 1.29* |
| counting bound, k=8 | 3.0 | — | 1.0 |

\* asymptotic (`n → ∞`). Exact peeling at `d = 10K` sits ~0.7d above it; that finite-size gap
shrinks as `d` grows. BP at finite size already matches the asymptotic coupled-peeling limit.

So at the paper's own scale, coupling + BP needs **44% fewer counters** and **~33% less
communication** than CommonSense (2.31× → 1.56× of the SetX lower bound). The gain splits about
evenly. The decoder upgrade alone (MP → BP on the same uniform matrix) gives 7.8d → 5.8d, and
coupling gives 5.8d → 4.4d. Coupling needs the new decoder: it makes MP worse (§2.3). Coupling + exact
peeling saves 25–29% while also giving a deterministic correctness certificate. The rest of the
gap to the lower bound splits into the finite-size gap, the MAP-vs-counting-bound gap, and the
CommonSense coder's 17–32% overhead over entropy.

At a second operating point (|B \ A| / |B| ≈ 1/21, run at 1/10 scale) the same picture holds:
2.48× (paper's MP) → 2.07× (BP) → **1.91×** (coupled + BP) of the lower bound (§5.1).

**Bidirectional SetX is different (§6).** In CommonSense's ping-pong protocol, coupling did not
reduce communication. Symmetric differences: best coupled 11.71 KB vs 11.19–11.35 KB uniform, at
1/10 scale; asymmetric 10:1 differences: 8.28 KB vs 7.71 KB. Coupling also roughly doubles the
round count wherever it lowers the counter count.
The other party's unknown elements remove the exact "zero counter ⇒ absent" test that makes
coupling's boundary easy, and the decoding wave only advances with alternations.

---

## 1. CommonSense through Walzer's lens

* `M1_A` is a counting Bloom filter (paper §3.3), i.e. exactly the hashing-based data structure
  Walzer's scheme improves: elements are balls, counters are bins, each ball is thrown into `k`
  bins.
* Bob's decoding problem is *quantitative group testing* over the candidate set `B`: each counter
  reports how many elements of `B \ A` hit it.
* **The paper's MP decoder operates at the threshold of a plain peeling decoder.** Exact peeling
  (row residual 0 ⇒ every undecided candidate in it is 0; residual = #undecided ⇒ all are 1) has
  density-evolution threshold `l/d = 7.38` for `k = 8`, `n/d = 101`, and simulates at ≈ 7.5; the
  paper runs MP at `l/d = 7.8` (`config/threshold.json`), and our faithful port of the MP decoder
  has its threshold at 7.8–8.0. So at its operating point CommonSense behaves like an
  IBLT-style peeling structure, and peeling is precisely the decoder that spatial coupling is
  proven to help (Walzer; Kudekar–Richardson–Urbanke threshold saturation).
* **There is a lot of room.** A counting bound — the `l` at which `l` counters can even carry
  `log2 C(n,d)` bits, `l · H(Poisson(kd/l)) ≥ LB` — is `l ≈ 3.0d` for `k = 8`. CommonSense uses
  `7.8d`, and its communication cost is 2.3× LB (22.9 KB vs 10.1 KB in the paper's Example 3; our
  cost model reproduces 2.31×).

## 2. Is spatial coupling applicable? (task a)

### 2.1 What carries over from Walzer / METTLE

| Ingredient | Walzer (IBLT) | METTLE (stream code) | CommonSense |
|---|---|---|---|
| Ball position | hash of element | sequence number | **hash of element** (encoding is blind: Alice and Bob must agree on an element's column without knowing each other's sets) |
| Window coupling | k bins uniform in `[pos, pos+w)` | same | **applies verbatim** — change `hash_all` only |
| Two-sided decoding ("both coasts") | yes | no (left-to-right) | **yes**: the decoder sees all counters at once |
| TLE (touch-less leading edge) | n/a | needs distinct integer positions | the *touch-less* property is impossible (hashed positions collide), but a deterministic leading edge at `pos` is still definable; slightly helpful for peeling at small scale, no effect on BP at paper scale |
| MET (non-identical edge offset laws) | n/a | exponentially decaying offsets | applicable; METTLE's law (tuned for one-sided decoding) was not better than TLE alone here |
| Tail compression | n/a | compress east coast (cannot decode from the right) | not useful (we decode from both ends); boundary rows are instead *cheap* to transmit because they are lightly loaded |

So the right model is **Walzer's scheme on the CBF sketch**, optionally with a deterministic
leading edge. METTLE's touch-less property depends on non-blind encoding, which SetX lacks.

### 2.2 Why it helps: threshold saturation (density evolution)

For exact peeling the DE is a two-dimensional recursion. Let `p_t` (`q_t`) be the probability
that a true (false) candidate's edge message is still undecided, `λ_T = kd/l`,
`λ_F = k(n−d)/l`:

```
q_{t+1} = (1 − exp(−λ_T p_t))^{k−1}      (a row resolves a false candidate if it has no other undecided true one)
p_{t+1} = (1 − exp(−λ_F q_t))^{k−1}      (a row resolves a true candidate if it has no other undecided false one)
```

The coupled version averages over the window (`peeling_de.py`). Thresholds at `n/d = 101`:

| k | counting bound l/d | uniform peeling l/d | coupled peeling l/d (bulk) | ideal cost / LB, uniform → coupled | CommonSense-coder cost / LB, uniform → coupled |
|---|---|---|---|---|---|
| 4 | 4.56 | 7.23 | 5.64 | 1.26 → **1.12** | 1.77 → 1.48 |
| 6 | 3.46 | 7.02 | 4.58 | 1.52 → **1.19** | 1.93 → **1.44** |
| 8 | 3.01 | 7.38 | 4.25 | 1.78 → 1.27 | 2.21 → 1.49 |
| 10 | 2.76 | 7.89 | 4.12 | 2.03 → 1.33 | — |
| 12 | 2.60 | 8.43 | 4.07 | 2.27 → 1.39 | — |

(Costs here are for a vanishing window fraction `w/L → 0`; with a 5% window the boundary adds
~5% rows, e.g. k = 8 becomes 4.46d and 1.57× with CommonSense's coder. Baseline for comparison:
CommonSense today, `k = 8`, `l/d = 7.8`: 1.84× LB ideal, **2.31×** with its coder.)

Two classic spatial-coupling signatures appear: the coupled threshold is far below the uncoupled one, and raising `k` *hurts* the uncoupled decoder but *helps* the coupled one, so
`k` must be re-optimised once coupling is in. BP (sum–product) is strictly stronger than peeling,
and coupling lets it approach the true MAP threshold, which lies between the counting bound and the
coupled-peeling threshold.

### 2.3 The catch: the existing MP decoder gets *worse* with coupling

Swapping in a coupled `M` but keeping the paper's MP decoder raised the threshold from ≈7.8d to
≈9d (|B| = 101K, d = 1K, k = 8; 0/24 successes at 8.0d with a 10% window, 0/24 at 8.5d with a 3%
window). The reason is coherence. Greedy pursuit needs columns that rarely overlap, and
localizing columns makes nearby columns overlap often. The expected number of false candidates
sharing ≥ 2 rows with a given column is

```
uniform:  n k^4 / (2 l^2)          coupled:  n k^4 / (3 W (l − W))      ratio ≈ (2/3)(l / W)
```

That is 13× more at `W = 5% l` and 22× at 3% (checked by Monte Carlo). Peeling and BP reason per
row and do not suffer from this. **So "spatial coupling for CommonSense" means coupling + a
message-passing decoder; the matrix change alone backfires.**

### 2.4 Correctness

**(i) What the paper's guarantee covers.** The paper argues "zero residue + RIP-1 of `M_B` ⇒
correct" (Theorems 7–8). But the RIP-1 condition of Theorem 7 needs every set `S` of ≤ 2d columns
to have ≥ (5/6)·k·|S| distinct neighbours, i.e. ≥ 13.3d rows at `k = 8`, while the experiments use
7.8d. Even the weaker unique-neighbour expansion (the minimum that rules out a binary/ternary
kernel vector) needs `l > k·d = 8d` for 2d-sized sets. Neither holds at the operating point, for
the uniform matrix or any other with `k` ones per column (pure counting: `2d·k/2 > l`). The
paper's correctness at practical sizes is therefore already **empirical / average-case**:
zero residue plus uniqueness of the sparse solution with high probability. That is the bar to
preserve, and the coupled scheme meets it and exceeds it as follows.

**(ii) Peeling gives a deterministic certificate, independent of the matrix.** Each peeling
decision is a logical consequence of `y`, `x ∈ {0,1}^B` and `x ≥ 0`. If peeling decides every
candidate, *every* binary solution of `Mx' = y` equals the peeled one, so it is unique and equal
to the truth. No RIP, expansion or probability is involved. Coupling is what makes peeling
complete at small `l`.

**(iii) BP output: zero residue + uniqueness.** When BP finishes the job, its output satisfies
`Mx' = y`; it is correct iff the solution is unique. For coupled matrices two facts keep this
local:
* *Local counts are pinned by `y`.* For any position interval `I = [a, b)` with `b − a > W`, every
  zero-residue solution has between `Σ_{j∈[a+W,b)} y_j / k` and `Σ_{j∈[a,b+W)} y_j / k` ones
  among columns positioned in `I`. Columns in `I` own all their rows inside `[a, b+W)`, and only
  they reach rows `[a+W, b)`. So a wrong solution must be locally balanced: it cannot pile up
  extra ones anywhere.
* *Small confusable structures are negligible.* The smallest kernel structure is a true and a false
  column with identical row sets. The expected number is ≈ `2·d·n·k!/((k+1)·l·W^{k−1})`, about
  6·10^-15 at paper scale with `W = 5% l`. Larger structures are rarer still. Above the MAP threshold,
  uniqueness holds w.h.p., as in channel coding.

**(iv) Make it checkable: a 64-bit set hash.** Alice appends `h(A) = Σ_{a∈A} g(a) mod 2^64` (`g`
a keyed hash independent of `M`). Bob accepts only if `h(B) − Σ_{x∈X̂} g(x) = h(A)`; a wrong
`X̂` passes with probability ≤ 2^-64. This costs 8 bytes and makes the protocol verifiable no
matter which matrix or decoder is used. We recommend it for CommonSense as is, too.

**(v) Empirically**, across every simulated configuration (uniform and coupled; MP, peeling, BP)
there were **0 silent errors**: no run ever reached a zero residue with a wrong set.
Counted over the logged runs: 114 configurations + 8 in `paper_sc_w05_tle.txt`, **1,816 trials, 0 silent errors** (`python3 analyze.py` prints the tally).

## 3. How to implement it (task b)

Done in this repository, opt-in, defaults unchanged:

1. **Matrix**: `DoroCode::enable_spatial_coupling(window, tle, seed)` in `lib/doro.hpp`. Every
   element is hashed to `s(x) ∈ [0, l − W]` and its `k` distinct counters are drawn from
   `[s(x), s(x) + W)` (optionally the first one fixed at `s(x)`). Only `hash_all` changes, so
   encoding, peeling, `copy()` and the streaming use are untouched.
2. **Decoder**: `lib/sc_decoder.hpp` (`ScBpDecoder`): exact peeling, then damped sum–product BP
   with exact Poisson-binomial check messages, then BP-guided decimation if needed. `|B \ A|` is
   read off the counters (`Σ y = k·|B \ A|` exactly). It only commits when the residue is zero; on
   failure the code is left untouched so the existing MP / IBLT fallback still runs.
   Damping matters: undamped flooding BP oscillates on this graph (confidence swinging between 4%
   and 99% of variables); with damping 0.5 it converges in tens of iterations.
3. **Compression parameters**: with coupling the bulk counters are more loaded
   (`λ_bulk = kd/(l − W + 1)`), so `two_party.cpp` now tunes `lb/ub/BCH` for the bulk. Boundary
   rows are lightly loaded; tuning per load class (`CSa` below) saves a few % more and is not yet
   implemented in `two_party.cpp`.
4. **Config keys** (`bin/two_party.cpp`): `"sc window"` (rows) or `"sc window fraction"`,
   `"sc tle"`, `"decoder": "bp"` (unidirectional). See `config/sc_unidirectional.json`.

## 4. Results in detail (paper scale: |B| = 1.01M, |B \ A| = 10K)

Success counts per measured `l/d` (16 trials for MP/peeling, 8 for coupled BP), from
`results/paper_uniform.jsonl`, `results/paper_sc.jsonl`, `results/paper_sc_w05_tle.txt` and
`results/peel_scan.jsonl`:

| k | matrix | decoder | successes by l/d | threshold |
|---|---|---|---|---|
| 8 | uniform | MP (paper) | 7.6: 4/16 · 7.8: 16/16 · 8.0: 16/16 | 7.8 |
| 8 | uniform | peeling | 7.2: 0 · 7.4: 11 · 7.6: 16/16 | 7.6 |
| 8 | uniform | BP | 5.6: 5 · 5.8: 16/16 · 6.0: 16/16 | 5.8 |
| 8 | coupled 5%, TLE | **MP (paper)** | 7.8: **0/8** · 8.4: 8/8 | **8.4 (worse)** |
| 8 | coupled 5%, TLE | peeling | 4.8: 1 · 5.0: 13 · 5.2: 16/16 | 5.2 |
| 8 | coupled 10%, TLE | peeling | 5.0: 7 · 5.2: 16/16 | 5.2 |
| 8 | coupled 2%, TLE | peeling | 5.0: 0 · 5.2: 11/16 | > 5.2 |
| 8 | coupled 5%, TLE | BP | 4.2: 0/8 · 4.4: 8/8 | 4.4 |
| 8 | coupled 5%, uniform offsets | BP | 4.2: 0/8 · 4.4: 8/8 | 4.4 |
| 8 | coupled 10%, TLE | BP | 4.2: 0/8 · 4.4: 8/8 · 4.6: 8/8 | 4.4 |
| 8 | coupled 2%, TLE | BP (1000 iters) | 4.4: 1/8 · 4.6: 8/8 | ≤ 4.6 (iteration-limited) |
| 6 | uniform | peeling | 7.0: 3 · 7.2: 16/16 | 7.2 |
| 6 | coupled 5%, TLE | peeling | 5.2: 2 · 5.4: 16/16 | 5.4 |
| 4 | uniform | peeling | 7.2: 2 · 7.4: 16/16 | 7.4 |
| 4 | coupled 5%, TLE | peeling | 6.4: 13/16 | > 6.4 |

## 5. Parameter findings

* **Window.** For BP, 5% and 10% of `l` give the same threshold (4.4d). A 2% window makes the
  decoding wave slower in iterations (≈700 iterations at 4.6d vs ≈200–330 for 5–10%) and hit the
  1000-iteration cap at 4.4d. For peeling, 2% is clearly worse at this size: the chain is about 50
  windows long, and the wave stalls on local fluctuations. Use 5%.
* **Leading edge (TLE).** No measurable effect on BP at paper scale; slightly helpful for peeling
  at small scale. Plain Walzer coupling is enough.
* **METTLE's MET offset law.** Not better than TLE alone at small scale. Its design goal
  (one-sided, left-to-right decoding) is not ours.
* **k.** In rows, larger `k` helps the coupled decoder. In bits, smaller `k` (sparser counters)
  is cheaper, so cost is minimised around `k = 6` for coupled peeling (1.64× vs 1.74× for k = 8).
  The paper's `k = 8` was tuned for the uncoupled MP decoder.
* **Runtime.** Peeling is ~10× faster than the paper's MP (0.1–0.8 s vs ~5 s per instance at
  paper scale, 8 threads). Coupled BP near threshold takes ~1–3.5 min single-threaded per instance in
  this unoptimised implementation, mostly iterations spent in already-decoded regions. A
  sliding-window schedule (standard for SC-LDPC) would cut that substantially.

### 5.1 Second operating point: |B \ A| / |B| ≈ 1/21

At 1/10 of the paper's third group (|B| = 105K, |B \ A| = 5K, k = 8, `results/nd21.jsonl`;
LB = 5.8 bits per element):

| Scheme | successes by l/d | threshold | cost / LB (CommonSense coder) |
|---|---|---|---|
| uniform + MP (paper) | 4.8: 0/16 · 5.1: 11/16 · 5.4: 16/16 | 5.4 | 2.48 |
| uniform + BP | 3.6: 0/8 · 3.9: 6/8 · 4.2: 16/16 | 4.2 | 2.07 |
| coupled 10% + BP | 3.0: 0/16 · 3.3: 14/16 · 3.6: 16/16 | 3.6 | **1.91** |

Here the decoder upgrade does more of the work (5.4d → 4.2d) and coupling adds 4.2d → 3.6d.
This point was run at a tenth of the paper's size, where coupling's finite-size penalty is larger.

### 5.2 k sweep (|B| = 101K, |B \ A| = 1K; small scale, so finite-size effects are large)

| k | uniform + MP threshold | coupled 10% + BP threshold | cost / LB at coupled threshold |
|---|---|---|---|
| 6 | > 8.5 | 5.5 | **1.69** |
| 8 | ≈ 7.8–8.2 | 5.0 | 1.72 |
| 10 | 8.0 | > 5.0 (12/16 at 5.0) | — |
| 12 | ≤ 8.0 (16/16 at the lowest tested, 8.0) | > 4.5 (0/16 at 3.5–4.5; higher not tested) | — |

The paper's MP decoder falls apart at small `k`, which is why it uses `k = 7–8`. With coupling +
BP, `k = 6` is at least as cheap as `k = 8` in bits, as density evolution predicts.

## 6. Bidirectional SetX (simulated): coupling does not pay in the ping-pong protocol

### 6.1 Setup

`sc_study/bidir_sim.cpp` ports the bidirectional loop of `bin/two_party.cpp`:
* Alice's first message (mod-q + BCH tuning, generalised to the Skellam difference);
* the shared residual ping-pong, Bob decoding first;
* the SMF fingerprints against common hallucinations, plus the collision-resolving stage;
* `two_party.cpp`'s status machine and signature lengths;
* the paper's bidirectional `k = 5`.

Two per-party decoders are compared:
* `mp`, a port of `DoroDecoder` including the k = 5 "5-in-3" repair;
* `bp`, sum–product BP over the party's own candidates. It treats the other party's undecoded
  elements as Poisson noise, with the rate moment-fitted per block of rows from the residual.
  It re-decides every round, with hysteresis: a new commitment needs LLR > 3, an existing one
  only LLR > 0.

Costs follow `two_party.cpp`: residuals are rANS-coded under a fitted Skellam, and fingerprints
are rANS-coded Bernoulli. Each scheme is also offered better residual coders and gets credit
for the cheapest per message: block-adaptive coding, and coding the residual's *change* since
the previous message.

**Validation.** At the paper's Example 11 setup (|A ∩ B| = 1M, |A \ B| = |B \ A| = 10K,
`k = 5`, `l = 8.5 |A △ B|`), the port decodes 8/8 in **7.0 rounds** (paper: 7.0–8.6) at
**111.6 KB** (paper: 129.2 KB). The ~14% gap is plausibly the paper's larger safety margin
(10,000 trials) and real encoder overheads; all comparisons below use the same cost model.

### 6.2 Results (|A ∩ B| = 100K, |A \ B| = |B \ A| = 1K, k = 5; `results/bidir_small_v3.jsonl`)

Cost = total of all messages with the cheapest residual coder (lower bound 2.02 KB).
Ratio = `l / |A △ B|`.

| Scheme | 5.5 | 6.0 | 7.0 | 8.5–9.0 |
|---|---|---|---|---|
| CommonSense MP, uniform | fail | fail | fail | **11.32–11.35 KB**, 7.0 rounds |
| BP, uniform | fail | 13.11 KB, 8.1 rounds | **11.19 KB**, 7.0 rounds | — |
| coupled 5% + BP | 2/8; 17.99 KB, 36 rounds | 14.99 KB, 13.9 rounds | 11.71 KB, 7.4 rounds | — |
| coupled 10% + BP | 3/8; 18.16 KB, 36 rounds | 15.02 KB, 14.5 rounds | 12.49 KB, 7.2 rounds | — |
| coupled 10% + MP | — | — | — | 8.5: 5/8 · 10.0: 11.95 KB |

With CommonSense's own residual coder the coupled runs cost more again (19.7 KB at 6.0, 48.7 KB
at 5.5). A 25% window was worse still (36.4 KB, 22.8 rounds at 6.0). There were 0 silent errors
in every configuration.

### 6.3 Why coupling does not help here

* **The other party's elements destroy the exact zero test.** One-directionally, a counter of 0
  proves every candidate in it is absent; that is what makes the lightly loaded boundary easy to
  peel. Bidirectionally a 0 can be `+1 − 1`. The other side's elements act as noise with the
  *same* spatial profile as your own signal, so the boundary is only mildly easier. A zero row
  carries an LLR of about `log λ_noise ≈ −1` instead of `−∞`.
* **The decoding wave needs alternations.** A party can only decode where the other party has
  already cleared its noise. Traces (`bidir_sim verbose=3`) show the front advancing about one
  coupling window every 2–4 rounds. At 5.0× the edge deciles go from 30% to 100% decoded over 13
  rounds while the bulk stays stuck. With a 5% window the chain is 20 windows long.
* **Rounds are expensive.** Every round is a residual message and a network round trip. Coding
  only the change in the residual cuts the coupled runs' cost by 25–60%, but not below the
  uniform matrix, and it does nothing for latency.
* **The decoder upgrade alone is a wash here.** The sketch is only ~45% of the bidirectional
  cost. BP lowers the threshold from 8.0× to 6.0×, but near that threshold it needs more rounds,
  so its best cost (11.19 KB) is within 1.5% of MP's (11.35 KB).

**Conclusion.** Spatial coupling substantially improves one-directional SetX (§4: 2.31× → 1.56×
of the lower bound). For symmetric bidirectional SetX with CommonSense's ping-pong protocol it
does not reduce communication and roughly doubles the round count. Making coupling work in two
directions needs a different protocol, one that does not rely on strict alternation to clear
noise. We did not run the bidirectional BP variants at full paper scale: they take 10–60
minutes per trial in this unoptimised simulator. The mechanism above is structural (the number
of windows the wave must cross does not shrink with scale), so we expect the same conclusion.

### 6.4 Asymmetric differences, and a control

Asymmetric differences, `|A \ B| = 200`, `|B \ A| = 2,000` (`results/bidir_asym*.jsonl`), are
the regime closest to one-directional: Bob decodes first, against light noise.

| Scheme | threshold (`l / |A △ B|`) | cost, best coder | rounds |
|---|---|---|---|
| CommonSense MP, uniform | 8.0 (fails 4–7) | 7.74 KB | 7.0 |
| BP, uniform | 5.0 (fails 3–4.5) | **7.71 KB** | 7.0 |
| coupled 5% + BP | 5.0 (fails 2.5–4.5) | 8.28 KB | 8.2 |

As a control, the same decoder with `|A \ B| = 0` (pure one-directional, no noise) does gain
from coupling: threshold 5.0 → 4.5 and cost 2.99 → 2.76 KB. So it is the other party's elements
that cancel coupling's advantage. Even 10% noise is enough at this scale. (This decoder has no
exact-peeling stage, so its one-directional gain is smaller than `sc_sim`'s in §4.)

## 7. Next steps

1. Re-optimise `k` and `W` at paper scale and above. DE says smaller `k` minimises cost, and at
   |B \ A| = 1K `k = 6` was cheapest (§5.2); finite-size gaps shrink with |B \ A|, so the
   Ethereum-sized instances (|A △ B| ≈ 0.9M) should land closer to the DE limits.
2. A better Slepian–Wolf coder for the first message. CommonSense's mod-q + BCH coder costs
   17–32% above the per-counter entropy, worst at the low loads coupling favours.
3. A bidirectional protocol that does not rely on strict alternation. In the ping-pong protocol
   coupling did not help (§6); asymmetric differences are the most promising regime to start from.
4. DE for BP (Krzakala et al. and Donoho–Javanmard–Montanari give the machinery for coupled
   CS matrices + message passing).
5. Streaming use (paper §4): the in-memory sketch shrinks by the same factor as `l`.

## 8. Related work worth reading (verified)

* Walzer, *Peeling close to the orientability threshold — spatial coupling in hashing-based data
  structures*, ACM TALG 21(3), 2025 (the professor's reading).
* Krzakala, Mézard, Sausset, Sun, Zdeborová, *Statistical-physics-based reconstruction in
  compressed sensing*, Phys. Rev. X 2, 021005 (2012) — "seeded" (spatially coupled) CS matrices +
  BP reach the information-theoretic limit; cited by METTLE as [26].
* Donoho, Javanmard, Montanari, *Information-theoretically optimal compressed sensing via spatial
  coupling and approximate message passing*, IEEE T-IT 2013 (arXiv 1112.0708) — the proof.
* Mashauri, Graell i Amat, Lentmaier, *Low-density parity-check codes and spatial coupling for
  quantitative group testing*, ISIT 2023 (arXiv 2404.11463) — spatial coupling + peeling for
  quantitative group testing with DE analysis; the closest prior work to what we do here.
  Follow-up: *Threshold saturation for quantitative group testing with low-density parity-check
  codes*, ISIT 2024. It proves threshold saturation (vector admissible systems / potential
  functions) and is the template for a rigorous version of §2.2.
* Karimi, Kazemi, Heidarzadeh, Narayanan, Sprintson, *Sparse graph codes for non-adaptive
  quantitative group testing*, arXiv 1901.07635.
* Coja-Oghlan, Gebhard, Hahn-Klimroth, Loick, *Optimal group testing*, COLT 2020 / CPC 2021 —
  spatially coupled test design (SPIV) reaching the information-theoretic bound.
* Hahn-Klimroth, Müller, *Near optimal efficient decoding from pooled data*, COLT 2022.
* Kudekar, Richardson, Urbanke, *Threshold saturation via spatial coupling*, IEEE T-IT 2011.

## Appendix: reproducing

```bash
cd sc_study
clang++ -O3 -std=c++17 -o sc_sim sc_sim.cpp -lpthread        # or g++ -O3 -std=c++17 -pthread
./sc_sim n=1000000 d=10000 k=8 ratio=7.8 matrix=uniform decoder=mp trials=16          # paper baseline
./sc_sim n=1000000 d=10000 k=8 ratio=4.8 matrix=sc wfrac=0.05 offsets=tle decoder=peel trials=16
./sc_sim n=1000000 d=10000 k=8 ratio=4.4 matrix=sc wfrac=0.05 offsets=tle decoder=bp trials=8
bash run_batch.sh paper_uniform   # etc.; see run_batch.sh for the grid
python3 analyze.py                # tables from results/*.jsonl
python3 peeling_de.py; python3 de_cost_sweep.py
```

Each `sc_sim` run prints one JSON line: successes, silent errors (zero residue but wrong set),
and the cost of the first message as `entropy_bits` (Σ_j H(Poisson(λ_j))), `cs_global_bits`
(CommonSense's mod-q + BCH coder ported from `doro_parameter_tuning`, one parameter set),
`cs_adaptive_bits` (same coder, parameters per row-load class), each also divided by
`lower_bound_bits = log2 C(|B|, |B \ A|)`.

The real-code path (`two_party` with `config/sc_unidirectional.json`) needs the original build
(g++-14, `cmake`). `test_sc_decoder.cpp` exercises the patched `DoroCode` and `ScBpDecoder`
directly and builds with Apple clang (`-std=c++2b`) given the wyhash headers and a stub for the
Eigen include in `libONIAK/global.h`.

## Appendix: files

| File | What it is |
|---|---|
| `sc_study/sc_sim.cpp` | standalone simulator: uniform / coupled matrices; MP (faithful port of `DoroDecoder`), peeling, BP; cost model (entropy, CommonSense coder) |
| `sc_study/peeling_de.py`, `de_cost_sweep.py` | density evolution for peeling, uncoupled and coupled |
| `sc_study/run_batch.sh`, `analyze.py`, `results/` | experiment grid and summaries |
| `sc_study/test_sc_decoder.cpp` | end-to-end test of the real `DoroCode` + `ScBpDecoder` |
| `sc_study/bidir_sim.cpp`, `run_bidir.sh` | bidirectional simulator: port of `two_party.cpp`'s ping-pong (MP, SMF, resolving, costs) plus a BP party decoder |
| `lib/doro.hpp` (modified) | `enable_spatial_coupling`, `hash_all_coupled` |
| `lib/sc_decoder.hpp` (new) | peeling + BP decoder |
| `bin/two_party.cpp` (modified) | config plumbing, bulk-load parameter tuning, `"decoder": "bp"` |
| `config/sc_unidirectional.json` | example config |
