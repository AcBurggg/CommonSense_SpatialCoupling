// sc_sim.cpp -- standalone simulator for unidirectional CommonSense (A subset of B)
// with uniform vs. spatially-coupled (SC) sparse binary CS matrices.
//
// Build:  clang++ -O3 -std=c++17 -o sc_sim sc_sim.cpp -lpthread
// Usage:  ./sc_sim key=value ...   (see parse_args / main for keys)
//
// What is simulated
//   * Bob holds B (n = |B| elements), Alice holds A subset of B, d = |B \ A|.
//   * Every element x gets k distinct rows of an l-row counter array (column m_x of M).
//   * Bob obtains y = M 1_{B\A} (= M1_B - M1_A) exactly (we assume the first-round
//     compression/BCH layer succeeded, as CommonSense does when it reports success).
//   * Bob decodes the binary d-sparse signal supported on B.
//
// Matrices
//   uniform : k distinct rows uniform in [0, l)                          (CommonSense, Def. 6)
//   sc      : position s uniform in {0..l-W}; k distinct rows s + o_t,  (Walzer-style coupling)
//             o_t drawn from an offset law on {0..W-1}:
//               offsets=uni    : all o_t ~ U{0..W-1}
//               offsets=tle    : o_0 = 0, others uniform                  (METTLE "TLE")
//               offsets=mettle : o_0 = 0, o_t = W-1 - Bin(W-1, 2^-t)     (METTLE "TLE+MET")
//
// Decoders
//   mp   : faithful re-implementation of CommonSense's DoroDecoder for the unidirectional
//          one-round case used in the paper's threshold experiments (L2 MP pass, then L1 pass,
//          max_num_peels = 2.5 d, same delta/clamp rules, same lazy priority queue semantics).
//   peel : exact deduction peeling (row residual 0 -> all unknowns are 0;
//          residual == #unknowns -> all unknowns are 1). Never makes a wrong decision.
//   bp   : peel, then sum-product BP on the residual graph (exact Poisson-binomial check
//          updates), BP-guided decimation if needed.
//   bpmp : bp, and if it does not reach zero residual, fall back to mp on the original problem.
//
// Success = decoded set equals B \ A. We also record "silent" errors: zero residual but wrong.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

// ----------------------------------------------------------------------------- RNG
static inline uint64_t splitmix64(uint64_t &x) {
  uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}
struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed) {}
  uint64_t next() { return splitmix64(s); }
  double unif() { return (next() >> 11) * 0x1.0p-53; }
  uint64_t below(uint64_t n) { return (uint64_t)(((__uint128_t)next() * n) >> 64); }
};

// ----------------------------------------------------------------------------- params
struct Params {
  long n_common = 1000000;  // |A cap B| = |A|
  long d = 10000;           // |B \ A|
  int k = 8;
  long l = 78000;           // rows
  std::string matrix = "uniform";
  long W = 0;               // window (rows) for sc
  double wfrac = 0.0;       // alternatively W = wfrac * l
  std::string offsets = "uni";
  std::string decoder = "mp";
  int trials = 10;
  int threads = 8;
  uint64_t seed = 1;
  int bp_iters = 400;
  int bp_rounds = 30;     // max BP-guided decimation rounds
  double bp_damp = 0.5;
  int verbose = 0;
};

// Binomial sample (small helper; W up to ~1e5 so direct sum of bits is too slow; use std).
static long binom_sample(Rng &rng, long n, double p) {
  std::mt19937_64 g(rng.next());
  std::binomial_distribution<long> bd(n, p);
  return bd(g);
}

// ----------------------------------------------------------------------------- matrix
struct Matrix {
  long n, l;
  int k;
  std::vector<int32_t> cols;      // n*k rows of each element
  std::vector<int64_t> row_start; // CSR row -> elements
  std::vector<int32_t> row_elems;
};

// Generate the k rows of one element. Deterministic in (key, matrix seed) like a hash.
static void element_rows(const Params &P, long W, uint64_t key, int32_t *out) {
  Rng rng(key * 0x9E3779B97F4A7C15ULL ^ 0xD1B54A32D192ED03ULL);
  const int k = P.k;
  if (P.matrix == "uniform") {
    for (int t = 0; t < k; ++t) {
      while (true) {
        int32_t r = (int32_t)rng.below(P.l);
        bool dup = false;
        for (int u = 0; u < t; ++u) dup |= (out[u] == r);
        if (!dup) { out[t] = r; break; }
      }
    }
    return;
  }
  // spatially coupled
  long npos = P.l - W + 1;
  long s = (long)rng.below(npos);
  for (int t = 0; t < k; ++t) {
    for (int attempt = 0;; ++attempt) {
      long o;
      if (P.offsets == "uni" || attempt > 50) {
        o = (long)rng.below(W);
      } else if (P.offsets == "tle") {
        o = (t == 0) ? 0 : (long)rng.below(W);
      } else if (P.offsets == "mettle") {
        if (t == 0) o = 0;
        else o = (W - 1) - binom_sample(rng, W - 1, std::ldexp(1.0, -t));
      } else {
        fprintf(stderr, "unknown offsets %s\n", P.offsets.c_str());
        exit(1);
      }
      int32_t r = (int32_t)(s + o);
      bool dup = false;
      for (int u = 0; u < t; ++u) dup |= (out[u] == r);
      if (!dup) { out[t] = r; break; }
    }
  }
}

static long window_of(const Params &P) {
  if (P.matrix == "uniform") return P.l;
  long W = P.W > 0 ? P.W : (long)std::llround(P.wfrac * P.l);
  if (W < P.k) W = P.k;
  if (W > P.l) W = P.l;
  return W;
}

// Exact per-row probability that edge t of a random element lands in row j (ignores the
// tiny distortion from duplicate rejection). Returns sum over t (i.e. expected #edges/elem).
static std::vector<double> row_edge_mass(const Params &P) {
  std::vector<double> mass(P.l, 0.0);
  if (P.matrix == "uniform") {
    for (auto &m : mass) m = (double)P.k / P.l;
    return mass;
  }
  long W = window_of(P);
  long npos = P.l - W + 1;
  for (int t = 0; t < P.k; ++t) {
    std::vector<double> pmf(W, 0.0);
    if (P.offsets == "uni" || (P.offsets != "uni" && t > 0 && P.offsets == "tle")) {
      for (auto &p : pmf) p = 1.0 / W;
    } else if (t == 0) {
      pmf[0] = 1.0;
    } else { // mettle: o = W-1 - Bin(W-1, 2^-t)
      double q = std::ldexp(1.0, -t);
      long nn = W - 1;
      for (long e = 0; e <= nn; ++e) {
        double lp = std::lgamma(nn + 1.0) - std::lgamma(e + 1.0) - std::lgamma(nn - e + 1.0) +
                    e * std::log(q) + (nn - e) * std::log1p(-q);
        pmf[nn - e] = std::exp(lp);
      }
    }
    std::vector<double> pre(W + 1, 0.0);
    for (long o = 0; o < W; ++o) pre[o + 1] = pre[o] + pmf[o];
    for (long j = 0; j < P.l; ++j) {
      long lo = std::max(0L, j - (npos - 1)), hi = std::min(W - 1, j);
      if (lo <= hi) mass[j] += (pre[hi + 1] - pre[lo]) / npos;
    }
  }
  return mass;
}

// ----------------------------------------------------------------------------- instance
struct Instance {
  long n, d, l;
  int k;
  Matrix M;
  std::vector<uint64_t> key;     // random element keys (priority-queue tie breaks)
  std::vector<uint8_t> truth;    // 1 if element in B \ A
  std::vector<int32_t> y;        // measurement M 1_{B\A}
};

static void build_instance(const Params &P, uint64_t trial_seed, Instance &I) {
  I.n = P.n_common + P.d;
  I.d = P.d;
  I.l = P.l;
  I.k = P.k;
  long W = window_of(P);
  Rng rng(trial_seed);
  I.key.resize(I.n);
  for (auto &k : I.key) k = rng.next();
  I.truth.assign(I.n, 0);
  // choose d true elements uniformly (keys are random anyway, so take a random subset)
  {
    std::vector<int32_t> idx(I.n);
    std::iota(idx.begin(), idx.end(), 0);
    for (long i = 0; i < I.d; ++i) {
      long j = i + (long)rng.below(I.n - i);
      std::swap(idx[i], idx[j]);
      I.truth[idx[i]] = 1;
    }
  }
  Matrix &M = I.M;
  M.n = I.n; M.l = I.l; M.k = I.k;
  M.cols.resize((size_t)I.n * I.k);
  for (long e = 0; e < I.n; ++e) element_rows(P, W, I.key[e], &M.cols[(size_t)e * I.k]);
  M.row_start.assign(I.l + 1, 0);
  for (auto r : M.cols) M.row_start[r + 1]++;
  for (long j = 0; j < I.l; ++j) M.row_start[j + 1] += M.row_start[j];
  M.row_elems.resize(M.cols.size());
  std::vector<int64_t> fill(M.row_start.begin(), M.row_start.end() - 1);
  for (long e = 0; e < I.n; ++e)
    for (int t = 0; t < I.k; ++t) M.row_elems[fill[M.cols[(size_t)e * I.k + t]]++] = (int32_t)e;
  I.y.assign(I.l, 0);
  for (long e = 0; e < I.n; ++e)
    if (I.truth[e])
      for (int t = 0; t < I.k; ++t) I.y[M.cols[(size_t)e * I.k + t]]++;
}

struct Result {
  bool zero_residual = false;
  bool correct = false;
  long peels = 0;
  long iters = 0;
  long undecided = 0;  // for peel-type decoders
  double secs = 0;
  std::string stage;   // which decoder stage finished
};

static bool check_correct(const Instance &I, const std::vector<int8_t> &x) {
  for (long e = 0; e < I.n; ++e)
    if ((x[e] != 0) != (I.truth[e] != 0)) return false;
  return true;
}

static bool check_zero_residual(const Instance &I, const std::vector<int8_t> &x) {
  std::vector<int32_t> r = I.y;
  for (long e = 0; e < I.n; ++e)
    if (x[e])
      for (int t = 0; t < I.k; ++t) r[I.M.cols[(size_t)e * I.k + t]] -= x[e];
  for (auto v : r)
    if (v) return false;
  return true;
}

// ----------------------------------------------------------------------------- MP (CommonSense)
// Mirrors Doro::DoroDecoder (lib/doro_decoder.hpp) for: CBF matrix (all signs +1), lb=0, ub=1,
// collision_resolving_ from the start (config "resolving round": -1), no quantization recentering
// (exact first-round transmission), max_num_peels = 2.5 d.
struct MPDecoder {
  const Instance &I;
  int k;
  std::vector<int32_t> r;
  std::vector<int8_t> val;
  std::vector<int32_t> str;
  struct Entry {
    int32_t prio;
    uint64_t key;
    int32_t idx;
    bool operator<(const Entry &o) const {
      return prio != o.prio ? prio < o.prio : key < o.key;
    }
  };
  std::priority_queue<Entry> pq;
  bool l2 = true;
  long peels = 0;
  long max_peels;
  long hard_cap;
  std::vector<int32_t> acc;
  std::vector<int32_t> touched;
  std::vector<uint8_t> is_touched;

  MPDecoder(const Instance &inst, std::vector<int32_t> residual, std::vector<int8_t> start)
      : I(inst), k(inst.k), r(std::move(residual)), val(std::move(start)), str(inst.n, 0),
        acc(inst.n, 0), is_touched(inst.n, 0) {
    max_peels = (long)(2.5 * I.d);
    hard_cap = 60L * std::max(1L, I.d) + 1000;
  }

  const int32_t *rows(long e) const { return &I.M.cols[(size_t)e * k]; }

  int32_t sense(long e) const {
    int32_t s = 0;
    const int32_t *rw = rows(e);
    for (int t = 0; t < k; ++t) s += r[rw[t]];
    return s;
  }
  // CommonSense's get_median, including its even-k convention.
  int32_t get_median(long e) const {
    int32_t sig[64];
    const int32_t *rw = rows(e);
    for (int t = 0; t < k; ++t) sig[t] = r[rw[t]];
    int32_t median = 0;
    int32_t *mid = sig + k / 2;
    std::nth_element(sig, mid, sig + k);
    if (k % 2 == 0) {
      int32_t *mid2 = mid - 1;
      std::nth_element(sig, mid2, mid);
      if ((long)(*mid) * (long)(*mid2) < 0) median = 0;
      else if (std::abs(*mid) > std::abs(*mid2)) median = *mid2;
    } else median = *mid;
    return median;
  }
  int32_t sense_l1(long e) const {
    const int32_t *rw = rows(e);
    int32_t med = get_median(e);
    double old_l1 = 0, new_l1 = 0;
    for (int t = 0; t < k; ++t) {
      old_l1 += std::abs(r[rw[t]]);
      new_l1 += std::abs(r[rw[t]] - med);
    }
    return (int32_t)(old_l1 - new_l1);
  }
  int32_t strength_to_delta(long e, int32_t strength) const {
    double normalized = l2 ? (double)strength / (double)k : (double)get_median(e);
    double cur = val[e];
    double nv = cur + normalized;
    nv = std::max(0.0, std::min(1.0, nv));
    double delta = nv - cur;
    if (std::fabs(delta) < 0.50001) return 0;
    return (int32_t)std::lround(delta);
  }
  void add_pq(long e, int32_t strength, int32_t delta) {
    str[e] = strength;
    if (delta != 0) pq.push({std::abs(strength), I.key[e], (int32_t)e});
  }
  bool pq_empty() {
    while (!pq.empty()) {
      const Entry &t = pq.top();
      if (t.prio == std::abs(str[t.idx])) return false;
      pq.pop();
    }
    return true;
  }
  void scan() {
    for (long e = 0; e < I.n; ++e) {
      int32_t s = l2 ? sense(e) : sense_l1(e);
      add_pq(e, s, strength_to_delta(e, s));
    }
  }
  void peel(long e, int32_t delta) {
    ++peels;
    const int32_t *rw = rows(e);
    for (int t = 0; t < k; ++t) r[rw[t]] -= delta;
    val[e] += (int8_t)delta;
    for (int t = 0; t < k; ++t) {
      int32_t j = rw[t];
      for (int64_t p = I.M.row_start[j]; p < I.M.row_start[j + 1]; ++p) {
        int32_t nb = I.M.row_elems[p];
        if (!is_touched[nb]) { is_touched[nb] = 1; touched.push_back(nb); acc[nb] = 0; }
        acc[nb] -= delta;
      }
    }
    for (int32_t nb : touched) {
      int32_t ns = l2 ? str[nb] + acc[nb] : sense_l1(nb);
      add_pq(nb, ns, strength_to_delta(nb, ns));
      is_touched[nb] = 0;
    }
    touched.clear();
  }
  // returns false if the safety cap was hit
  bool peel_until_empty() {
    while (!pq_empty()) {
      Entry t = pq.top();
      int32_t strength = str[t.idx];
      int32_t delta = strength_to_delta(t.idx, strength);
      if (delta == 0) { pq.pop(); continue; }
      peel(t.idx, delta);
      if (peels > hard_cap) return false;
    }
    return true;
  }
  void decode(int max_rounds = 100) {
    l2 = true;
    scan();
    for (int rnd = 0; rnd < max_rounds; ++rnd) {
      if (!peel_until_empty()) return;
      if (l2) { l2 = false; scan(); }
      if (peels >= max_peels) break;
      if (pq_empty()) break;
    }
  }
  bool zero_residual() const {
    for (auto v : r)
      if (v) return false;
    return true;
  }
};

static Result run_mp(const Instance &I) {
  Result R;
  MPDecoder D(I, I.y, std::vector<int8_t>(I.n, 0));
  D.decode();
  R.peels = D.peels;
  R.zero_residual = D.zero_residual();
  R.correct = check_correct(I, D.val);
  R.stage = "mp";
  return R;
}

// ----------------------------------------------------------------------------- exact peeling
// state: -1 unknown, 0 known zero, 1 known one
struct Peeler {
  const Instance &I;
  int k;
  std::vector<int8_t> state;
  std::vector<int32_t> r;  // residual: y - #known ones
  std::vector<int32_t> u;  // #unknown in row
  std::vector<int32_t> q;  // row queue
  std::vector<uint8_t> inq;
  long n_unknown;
  long n_ones = 0;

  explicit Peeler(const Instance &inst)
      : I(inst), k(inst.k), state(inst.n, -1), r(inst.y), u(inst.l, 0), inq(inst.l, 0),
        n_unknown(inst.n) {
    for (long j = 0; j < I.l; ++j) u[j] = (int32_t)(I.M.row_start[j + 1] - I.M.row_start[j]);
  }
  void maybe_push(int32_t j) {
    if (!inq[j] && u[j] > 0 && (r[j] == 0 || r[j] == u[j])) { inq[j] = 1; q.push_back(j); }
  }
  void set(long e, int8_t v) {
    state[e] = v;
    --n_unknown;
    if (v) ++n_ones;
    const int32_t *rw = &I.M.cols[(size_t)e * k];
    for (int t = 0; t < k; ++t) {
      int32_t j = rw[t];
      u[j]--;
      r[j] -= v;
      maybe_push(j);
    }
  }
  void run() {
    for (long j = 0; j < I.l; ++j) maybe_push((int32_t)j);
    while (!q.empty()) {
      int32_t j = q.back();
      q.pop_back();
      inq[j] = 0;
      if (u[j] == 0) continue;
      int8_t v;
      if (r[j] == 0) v = 0;
      else if (r[j] == u[j]) v = 1;
      else continue;
      for (int64_t p = I.M.row_start[j]; p < I.M.row_start[j + 1]; ++p) {
        int32_t e = I.M.row_elems[p];
        if (state[e] < 0) set(e, v);
      }
    }
  }
};

static Result run_peel(const Instance &I, const Params &P) {
  Result R;
  Peeler Pe(I);
  Pe.run();
  R.undecided = Pe.n_unknown;
  if (P.verbose > 1) {
    // where along the row axis (by the element's first row) do undecided candidates sit?
    std::vector<long> und(20, 0), tot(20, 0);
    for (long e = 0; e < I.n; ++e) {
      int32_t lo = I.l;
      for (int t = 0; t < I.k; ++t) lo = std::min(lo, I.M.cols[(size_t)e * I.k + t]);
      int b = (int)((long)lo * 20 / I.l);
      tot[b]++;
      und[b] += Pe.state[e] < 0;
    }
    fprintf(stderr, "  undecided by 5%% slice:");
    for (int b = 0; b < 20; ++b) fprintf(stderr, " %.2f", tot[b] ? (double)und[b] / tot[b] : 0.0);
    fprintf(stderr, "\n");
  }
  std::vector<int8_t> x(I.n, 0);
  for (long e = 0; e < I.n; ++e) x[e] = Pe.state[e] == 1;
  R.zero_residual = Pe.n_unknown == 0 && check_zero_residual(I, x);
  R.correct = R.zero_residual && check_correct(I, x);
  R.stage = "peel";
  return R;
}

// ----------------------------------------------------------------------------- BP
// Sum-product on the residual graph left by exact peeling. Checks: sum of unknown x's in
// row j equals residual r_j. Exact leave-one-out Poisson-binomial via prefix/suffix DP
// truncated to [0, r_j]. Hard decisions are verified against all residuals each iteration.
// If BP stalls, BP-guided decimation fixes the most confident variables and re-peels.
static Result run_bp(const Instance &I, const Params &P) {
  Result R;
  Peeler Pe(I);
  Pe.run();
  const int k = I.k;
  long total_iters = 0;
  for (int round = 0; round < P.bp_rounds && Pe.n_unknown > 0; ++round) {
    // collect unknown vars & active rows
    std::vector<int32_t> var_id(I.n, -1);
    std::vector<int32_t> vars;
    for (long e = 0; e < I.n; ++e)
      if (Pe.state[e] < 0) { var_id[e] = (int32_t)vars.size(); vars.push_back((int32_t)e); }
    const long nv = (long)vars.size();
    std::vector<int32_t> rows_act;
    std::vector<int32_t> row_id(I.l, -1);
    for (long j = 0; j < I.l; ++j)
      if (Pe.u[j] > 0) { row_id[j] = (int32_t)rows_act.size(); rows_act.push_back((int32_t)j); }
    // edges: var-major, slot t of var v is edge v*k+t (or -1 if row inactive; cannot be since
    // an unknown var's rows all have u>0)
    std::vector<int64_t> rstart(rows_act.size() + 1, 0);
    for (long v = 0; v < nv; ++v)
      for (int t = 0; t < k; ++t) rstart[row_id[I.M.cols[(size_t)vars[v] * k + t]] + 1]++;
    for (size_t a = 0; a < rows_act.size(); ++a) rstart[a + 1] += rstart[a];
    std::vector<int64_t> redge(rstart.back());
    {
      std::vector<int64_t> fill(rstart.begin(), rstart.end() - 1);
      for (long v = 0; v < nv; ++v)
        for (int t = 0; t < k; ++t) {
          int32_t a = row_id[I.M.cols[(size_t)vars[v] * k + t]];
          redge[fill[a]++] = (int64_t)v * k + t;
        }
    }
    const double pi = std::min(0.999, std::max(1e-9, (double)(I.d - Pe.n_ones) / (double)nv));
    const double prior = std::log(pi / (1 - pi));
    std::vector<float> q2c((size_t)nv * k, (float)pi);  // P(x=1) var->check
    std::vector<float> c2v((size_t)nv * k, 0.f);        // LLR check->var
    std::vector<double> tot(nv, prior);
    std::vector<double> pre, suf;
    bool solved = false;
    double LCLAMP = 30.0;
    int it = 0;
    for (; it < P.bp_iters; ++it) {
      // check update
      for (size_t a = 0; a < rows_act.size(); ++a) {
        int32_t rj = Pe.r[rows_act[a]];
        int64_t s0 = rstart[a], s1 = rstart[a + 1];
        long deg = s1 - s0;
        int S = rj + 1;
        pre.assign((size_t)(deg + 1) * S, 0.0);
        suf.assign((size_t)(deg + 1) * S, 0.0);
        pre[0] = 1.0;
        for (long t = 0; t < deg; ++t) {
          double qq = q2c[redge[s0 + t]];
          double *cur = &pre[t * S], *nxt = &pre[(t + 1) * S];
          double mx = 0;
          for (int s = 0; s < S; ++s) {
            nxt[s] = cur[s] * (1 - qq) + (s ? cur[s - 1] * qq : 0.0);
            mx = std::max(mx, nxt[s]);
          }
          if (mx > 0) for (int s = 0; s < S; ++s) nxt[s] /= mx;
        }
        suf[deg * S] = 1.0;
        for (long t = deg - 1; t >= 0; --t) {
          double qq = q2c[redge[s0 + t]];
          double *cur = &suf[(t + 1) * S], *nxt = &suf[t * S];
          double mx = 0;
          for (int s = 0; s < S; ++s) {
            nxt[s] = cur[s] * (1 - qq) + (s ? cur[s - 1] * qq : 0.0);
            mx = std::max(mx, nxt[s]);
          }
          if (mx > 0) for (int s = 0; s < S; ++s) nxt[s] /= mx;
        }
        for (long t = 0; t < deg; ++t) {
          const double *A = &pre[t * S], *B = &suf[(t + 1) * S];
          double m0 = 0, m1 = 0;  // P(S_-t = rj), P(S_-t = rj-1)
          for (int s = 0; s <= rj; ++s) m0 += A[s] * B[rj - s];
          for (int s = 0; s <= rj - 1; ++s) m1 += A[s] * B[rj - 1 - s];
          double L;
          if (m1 <= 0 && m0 <= 0) L = 0;
          else if (m1 <= 0) L = -LCLAMP;
          else if (m0 <= 0) L = LCLAMP;
          else L = std::max(-LCLAMP, std::min(LCLAMP, std::log(m1) - std::log(m0)));
          c2v[redge[s0 + t]] = (float)L;
        }
      }
      // variable update
      for (long v = 0; v < nv; ++v) {
        double s = prior;
        for (int t = 0; t < k; ++t) s += c2v[(size_t)v * k + t];
        tot[v] = s;
        for (int t = 0; t < k; ++t) {
          double ext = s - c2v[(size_t)v * k + t];
          double p = 1.0 / (1.0 + std::exp(-ext));
          float &qq = q2c[(size_t)v * k + t];
          qq = (float)(P.bp_damp * qq + (1 - P.bp_damp) * p);
        }
      }
      // hard decision check
      bool ok = true;
      for (size_t a = 0; a < rows_act.size() && ok; ++a) {
        int32_t cnt = 0;
        for (int64_t p = rstart[a]; p < rstart[a + 1]; ++p) cnt += tot[redge[p] / k] > 0;
        ok = cnt == Pe.r[rows_act[a]];
      }
      if (ok) { solved = true; ++it; break; }
      if (P.verbose > 2 && it % 25 == 0) {
        // progress profile: fraction of confident variables per tenth of the row range
        long conf = 0, unsat = 0;
        for (long v = 0; v < nv; ++v) conf += std::fabs(tot[v]) > 8;
        std::vector<long> bad(10, 0);
        for (size_t a = 0; a < rows_act.size(); ++a) {
          int32_t cnt = 0;
          for (int64_t p = rstart[a]; p < rstart[a + 1]; ++p) cnt += tot[redge[p] / k] > 0;
          if (cnt != Pe.r[rows_act[a]]) { ++unsat; bad[(long)rows_act[a] * 10 / I.l]++; }
        }
        fprintf(stderr, "    it %4d confident=%ld/%ld unsat=%ld  by-decile:", it, conf, nv, unsat);
        for (auto b : bad) fprintf(stderr, " %ld", b);
        fprintf(stderr, "\n");
      }
    }
    total_iters += it;
    if (P.verbose > 1) {
      long unsat = 0, ones = 0;
      for (size_t a = 0; a < rows_act.size(); ++a) {
        int32_t cnt = 0;
        for (int64_t p = rstart[a]; p < rstart[a + 1]; ++p) cnt += tot[redge[p] / k] > 0;
        unsat += cnt != Pe.r[rows_act[a]];
      }
      for (long v = 0; v < nv; ++v) ones += tot[v] > 0;
      fprintf(stderr, "  bp round %d: unknown=%ld rows=%zu iters=%d unsat_rows=%ld hard_ones=%ld "
              "need_ones=%ld solved=%d\n", round, nv, rows_act.size(), it, unsat, ones,
              I.d - Pe.n_ones, solved);
    }
    if (solved) {
      for (long v = 0; v < nv; ++v) Pe.set(vars[v], tot[v] > 0 ? 1 : 0);
      break;
    }
    // decimation: fix the most confident ~2% of unknowns (at least 1) to their BP decision,
    // preferring confident zeros (cheap, rarely wrong), then re-peel.
    std::vector<std::pair<double, int32_t>> conf;
    conf.reserve(nv);
    for (long v = 0; v < nv; ++v) conf.push_back({-std::fabs(tot[v]), (int32_t)v});
    long nfix = std::max(1L, nv / 50);
    std::partial_sort(conf.begin(), conf.begin() + nfix, conf.end());
    for (long i = 0; i < nfix; ++i) {
      int32_t e = vars[conf[i].second];
      if (Pe.state[e] < 0) Pe.set(e, tot[conf[i].second] > 0 ? 1 : 0);
    }
    Pe.run();
    // inconsistent rows -> give up
    bool bad = false;
    for (long j = 0; j < I.l && !bad; ++j) bad = Pe.r[j] < 0 || Pe.r[j] > Pe.u[j];
    if (bad) break;
  }
  R.iters = total_iters;
  R.undecided = Pe.n_unknown;
  std::vector<int8_t> x(I.n, 0);
  for (long e = 0; e < I.n; ++e) x[e] = Pe.state[e] == 1;
  R.zero_residual = Pe.n_unknown == 0 && check_zero_residual(I, x);
  R.correct = R.zero_residual && check_correct(I, x);
  R.stage = "bp";
  return R;
}

static Result run_decoder(const Instance &I, const Params &P) {
  if (P.decoder == "mp") return run_mp(I);
  if (P.decoder == "peel") return run_peel(I, P);
  if (P.decoder == "bp") return run_bp(I, P);
  if (P.decoder == "bpmp") {
    Result R = run_bp(I, P);
    if (R.zero_residual) return R;
    Result M = run_mp(I);
    M.stage = "bp->mp";
    M.iters = R.iters;
    return M;
  }
  fprintf(stderr, "unknown decoder %s\n", P.decoder.c_str());
  exit(1);
}

// ----------------------------------------------------------------------------- costs
static double H_poisson(double lam) {
  if (lam <= 0) return 0;
  double H = 0;
  double lp = -lam;  // log pmf(0)
  for (int x = 0; x < 2000; ++x) {
    if (x > 0) lp += std::log(lam) - std::log((double)x);
    double p = std::exp(lp);
    if (p > 0) H -= p * lp / std::log(2.0);
    if (x > lam && p < 1e-18) break;
  }
  return H;
}

// CommonSense first-round cost per counter for a Poisson(lam) difference, ported from
// doro_parameter_tuning() in bin/two_party.cpp (mu1 = 0, i.e. unidirectional).
static double binom_quantile_up(int n, double p, double q) {
  double lp = n * std::log1p(-p), cdf = 0;
  for (int c = 0; c <= n; ++c) {
    if (c > 0) lp += std::log((double)(n - c + 1) / c) + std::log(p) - std::log1p(-p);
    cdf += std::exp(lp);
    if (cdf >= q - 1e-12) return c;
  }
  return n;
}
static double cs_cost_per_counter(double lam) {
  if (lam <= 1e-12) return 0.0;
  // pmf / cdf of Poisson
  std::vector<double> pmf;
  {
    double lp = -lam;
    for (int x = 0;; ++x) {
      if (x > 0) lp += std::log(lam) - std::log((double)x);
      double p = std::exp(lp);
      if (x > lam && p < 1e-18) break;
      pmf.push_back(p);
    }
  }
  int K1 = (int)pmf.size();  // support [0, K1)
  auto cdf = [&](int x) {    // P(X < x)
    if (x <= 0) return 0.0;
    double s = 0;
    for (int i = 0; i < std::min(x, K1); ++i) s += pmf[i];
    return std::min(1.0, s);
  };
  const double dmin = 1e-9, dmax = 0.1, dfinal = 1e-6, bch_block_err = 0.01;
  int mean = (int)std::floor(lam);
  double best = INFINITY;
  for (int ub = mean + 1; ub <= K1; ++ub) {
    for (int lb = mean; lb > -1; --lb) {
      double diff_err = 1.0 - cdf(ub) + cdf(lb);
      if (diff_err < dmin || diff_err > dmax) continue;
      int interval = ub - lb, cur_lb = lb, cur_ub = ub;
      double cost = std::log2((double)(ub - lb));
      bool okk = true;
      while (diff_err > dfinal) {
        double pin = 0;
        int midpoint = 0;
        for (int mid = cur_lb; mid <= cur_ub; ++mid) {
          double c = cdf(mid + interval) - cdf(mid - interval);
          if (c > pin) { pin = c; midpoint = mid; }
        }
        double rate = INFINITY;
        for (int order = 5; order <= 15; ++order) {
          int len = (1 << order) - 1;
          int cap = (int)std::ceil(binom_quantile_up(len, diff_err, 1.0 - bch_block_err));
          if (cap < 2) cap = 2;
          double chk = (double)order * cap;
          if (len <= chk) continue;
          rate = std::min(rate, chk / (len - chk));
        }
        if (!std::isfinite(rate)) { okk = false; break; }
        cost += rate;
        interval *= 2;
        diff_err = 1.0 - pin;
        cur_lb = midpoint - interval;
        cur_ub = midpoint + interval;
      }
      if (!okk || diff_err > dfinal) continue;
      best = std::min(best, cost);
    }
  }
  return best;
}

struct Costs {
  double entropy_bits;     // sum_j H(Poisson(lambda_j))  (ideal Slepian-Wolf-style)
  double cs_global_bits;   // l * cs_cost(max lambda)      (CommonSense scheme, one q for all rows)
  double cs_adaptive_bits; // sum_j cs_cost(lambda_j)      (CommonSense scheme, q per row class)
  double lambda_bulk;
};

static Costs compute_costs(const Params &P) {
  Costs C{};
  std::vector<double> mass = row_edge_mass(P);
  std::map<long, std::pair<double, long>> classes;  // quantized lambda -> (lambda, count)
  double lam_max = 0;
  for (long j = 0; j < P.l; ++j) {
    double lam = mass[j] * P.d;
    lam_max = std::max(lam_max, lam);
    long key = std::llround(lam * 2000);
    auto &c = classes[key];
    c.first = lam;
    c.second++;
  }
  for (auto &kv : classes) {
    double lam = kv.second.first;
    long cnt = kv.second.second;
    C.entropy_bits += cnt * H_poisson(lam);
    C.cs_adaptive_bits += cnt * cs_cost_per_counter(lam);
  }
  C.cs_global_bits = P.l * cs_cost_per_counter(lam_max);
  C.lambda_bulk = lam_max;
  return C;
}

static double log2_binom(double n, double k) {
  return (std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1)) / std::log(2.0);
}

// ----------------------------------------------------------------------------- main
static Params parse_args(int argc, char **argv) {
  Params P;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto eq = a.find('=');
    if (eq == std::string::npos) { fprintf(stderr, "bad arg %s\n", a.c_str()); exit(1); }
    std::string key = a.substr(0, eq), v = a.substr(eq + 1);
    if (key == "n") P.n_common = std::stol(v);
    else if (key == "d") P.d = std::stol(v);
    else if (key == "k") P.k = std::stoi(v);
    else if (key == "l") P.l = std::stol(v);
    else if (key == "ratio") P.l = std::lround(std::stod(v) * P.d);  // l = ratio * d (set d first)
    else if (key == "matrix") P.matrix = v;
    else if (key == "W") P.W = std::stol(v);
    else if (key == "wfrac") P.wfrac = std::stod(v);
    else if (key == "offsets") P.offsets = v;
    else if (key == "decoder") P.decoder = v;
    else if (key == "trials") P.trials = std::stoi(v);
    else if (key == "threads") P.threads = std::stoi(v);
    else if (key == "seed") P.seed = std::stoull(v);
    else if (key == "bp_iters") P.bp_iters = std::stoi(v);
    else if (key == "bp_damp") P.bp_damp = std::stod(v);
    else if (key == "bp_rounds") P.bp_rounds = std::stoi(v);
    else if (key == "verbose") P.verbose = std::stoi(v);
    else { fprintf(stderr, "unknown key %s\n", key.c_str()); exit(1); }
  }
  return P;
}

int main(int argc, char **argv) {
  Params P = parse_args(argc, argv);
  if (P.matrix == "uniform") P.offsets = "-";
  long W = window_of(P);
  std::atomic<int> next{0};
  std::mutex mu;
  int n_success = 0, n_zero = 0, n_silent = 0;
  double sum_secs = 0, sum_peels = 0, sum_iters = 0, sum_undec = 0;
  std::map<std::string, int> stage_count;
  auto worker = [&]() {
    Instance I;
    while (true) {
      int t = next.fetch_add(1);
      if (t >= P.trials) break;
      uint64_t s = P.seed * 1000003ULL + (uint64_t)t * 7919ULL + 17;
      build_instance(P, splitmix64(s), I);
      auto t0 = std::chrono::steady_clock::now();
      Result R = run_decoder(I, P);
      R.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      std::lock_guard<std::mutex> g(mu);
      n_success += R.correct;
      n_zero += R.zero_residual;
      n_silent += (R.zero_residual && !R.correct);
      sum_secs += R.secs;
      sum_peels += R.peels;
      sum_iters += R.iters;
      sum_undec += R.undecided;
      stage_count[R.stage + (R.correct ? ":ok" : ":fail")]++;
      if (P.verbose)
        fprintf(stderr, "trial %d ok=%d zero=%d peels=%ld iters=%ld undec=%ld %.2fs %s\n", t,
                R.correct, R.zero_residual, R.peels, R.iters, R.undecided, R.secs, R.stage.c_str());
    }
  };
  std::vector<std::thread> th;
  int nt = std::max(1, std::min(P.threads, P.trials));
  for (int i = 0; i < nt; ++i) th.emplace_back(worker);
  for (auto &t : th) t.join();

  Costs C = compute_costs(P);
  double n = (double)(P.n_common + P.d);
  double lb_bits = log2_binom(n, (double)P.d);
  std::string stages;
  for (auto &kv : stage_count) stages += (stages.empty() ? "" : ",") + kv.first + ":" + std::to_string(kv.second);
  printf("{\"n\":%ld,\"d\":%ld,\"k\":%d,\"l\":%ld,\"ratio\":%.4f,\"matrix\":\"%s\",\"W\":%ld,"
         "\"wfrac\":%.4f,\"offsets\":\"%s\",\"decoder\":\"%s\",\"trials\":%d,\"success\":%d,"
         "\"zero_residual\":%d,\"silent_errors\":%d,\"avg_secs\":%.3f,\"avg_peels\":%.1f,"
         "\"avg_iters\":%.1f,\"avg_undecided\":%.1f,\"lambda_bulk\":%.4f,\"entropy_bits\":%.0f,"
         "\"cs_global_bits\":%.0f,\"cs_adaptive_bits\":%.0f,\"lower_bound_bits\":%.0f,"
         "\"entropy_over_lb\":%.4f,\"cs_global_over_lb\":%.4f,\"cs_adaptive_over_lb\":%.4f,"
         "\"stages\":\"%s\"}\n",
         P.n_common + P.d, P.d, P.k, P.l, (double)P.l / P.d, P.matrix.c_str(),
         P.matrix == "uniform" ? 0L : W, P.matrix == "uniform" ? 0.0 : (double)W / P.l,
         P.offsets.c_str(), P.decoder.c_str(), P.trials, n_success, n_zero, n_silent,
         sum_secs / std::max(1, P.trials), sum_peels / std::max(1, P.trials),
         sum_iters / std::max(1, P.trials), sum_undec / std::max(1, P.trials),
         C.lambda_bulk, C.entropy_bits, C.cs_global_bits, C.cs_adaptive_bits, lb_bits,
         C.entropy_bits / lb_bits, C.cs_global_bits / lb_bits, C.cs_adaptive_bits / lb_bits,
         stages.c_str());
  return 0;
}
