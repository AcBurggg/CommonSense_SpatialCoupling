// bidir_sim.cpp -- bidirectional CommonSense (neither A nor B contains the other) with uniform or
// spatially coupled (SC) counting sketches.
//
// Build:  clang++ -O3 -std=c++17 -o bidir_sim bidir_sim.cpp -lpthread
// Usage:  ./bidir_sim key=value ...   (see parse_args)
//
// Protocol: a port of the ping-pong loop in bin/two_party.cpp. Alice sends M1_A (cost: the
// mod-q + BCH tuning of doro_parameter_tuning, here for a Skellam difference); both parties then
// share the residual r = M1_B - M1_A - M(x_B + x_A), and alternately decode (Bob first), send the
// residual (cost: rANS with a moment-fitted Skellam, as in two_party.cpp), send fingerprints of
// newly decoded elements (the SMF that prevents common hallucinations), and, once in the
// collision-resolving stage, send (finger1, finger2) pairs for decoded elements that hit the other
// party's fingerprints. Same status machine, resolving round and max rounds as two_party.cpp.
//
// Per-party decoders (decoder=...):
//   mp : port of Doro::DoroDecoder (L2 MP, then L1 once resolving; k = 5 "5-in-3" repair; peel cap
//        2.5 |A delta B| on the shared peel counter), i.e. CommonSense as published.
//   bp : sum-product BP over the party's own candidates. The other party's not-yet-decoded unique
//        elements are modelled as Poisson noise subtracted from each row, with the noise rate
//        moment-fitted per block of rows from the party's view of the residual. Messages persist
//        across rounds (warm start); the decoded set is re-decided every round, so earlier
//        mistakes can be undone. Elements blocked by the other party's fingerprints are clamped to 0
//        during collision avoidance, exactly like the MP decoder.
//
// Success = both decoded sets exactly right. A "silent error" is zero residual with a wrong set.

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
#include <queue>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

// ----------------------------------------------------------------------------- RNG / hashing
static inline uint64_t splitmix64(uint64_t &x) {
  uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}
static inline uint64_t mix(uint64_t x) { return splitmix64(x); }
struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed) {}
  uint64_t next() { return splitmix64(s); }
  uint64_t below(uint64_t n) { return (uint64_t)(((__uint128_t)next() * n) >> 64); }
};

// ----------------------------------------------------------------------------- params
struct Params {
  long n_common = 100000, dA = 1000, dB = 1000;
  int k = 5;
  long l = 17000;
  std::string matrix = "uniform";
  long W = 0;
  double wfrac = 0.0;
  std::string offsets = "uni";
  std::string decoder = "mp";
  int trials = 8, threads = 8;
  uint64_t seed = 1;
  int max_rounds = 10;        // "max comm rounds" in threshold3.json
  int resolving_round = 3;    // two_party.cpp default
  double failure_rate = 1e-7; // two_party.cpp default (divided by |B \ A| like the original)
  int bp_iters = 400, bp_min_iters = 10, bp_stable = 3;
  double bp_tol = 0.05;       // stop when max |change in LLR| < tol for bp_stable iterations
  double bp_damp = 0.5;
  double commit_llr = 3.0;    // BP: newly commit only when LLR > this; keep while LLR > 0
  long block_rows = 0;        // BP noise-fit / adaptive-coding block; 0 = auto (W/4 or l/64)
  int verbose = 0;
};

// ----------------------------------------------------------------------------- matrix
static long window_of(const Params &P) {
  if (P.matrix == "uniform") return P.l;
  long W = P.W > 0 ? P.W : (long)std::llround(P.wfrac * P.l);
  return std::max<long>(std::min<long>(W, P.l), 4L * P.k);
}

static void element_rows(const Params &P, long W, uint64_t key, int32_t *out) {
  Rng rng(key * 0x9E3779B97F4A7C15ULL ^ 0xD1B54A32D192ED03ULL);
  const int k = P.k;
  long base = 0, range = P.l;
  bool tle = false;
  if (P.matrix != "uniform") {
    base = (long)rng.below(P.l - W + 1);
    range = W;
    tle = P.offsets == "tle";
  }
  for (int t = 0; t < k; ++t) {
    while (true) {
      int32_t r = (int32_t)(base + ((tle && t == 0) ? 0 : (long)rng.below(range)));
      bool dup = false;
      for (int u = 0; u < t; ++u) dup |= out[u] == r;
      if (!dup) { out[t] = r; break; }
    }
  }
}

// expected edges per element landing in each row (sum over the k edges)
static std::vector<double> row_edge_mass(const Params &P) {
  std::vector<double> mass(P.l, 0.0);
  if (P.matrix == "uniform") {
    for (auto &m : mass) m = (double)P.k / P.l;
    return mass;
  }
  long W = window_of(P), npos = P.l - W + 1;
  for (long j = 0; j < P.l; ++j) {
    long lo = std::max(0L, j - (npos - 1)), hi = std::min(W - 1, j);
    double frac = lo <= hi ? (double)(hi - lo + 1) / ((double)npos * W) : 0.0;  // one uniform edge
    double m = 0;
    if (P.offsets == "tle") m = (P.k - 1) * frac + ((j < npos) ? 1.0 / npos : 0.0);
    else m = P.k * frac;
    mass[j] = m;
  }
  return mass;
}

// ----------------------------------------------------------------------------- instance
enum : uint8_t { COMMON = 0, A_ONLY = 1, B_ONLY = 2 };
struct Instance {
  long N, l;
  int k;
  std::vector<uint64_t> key;
  std::vector<uint8_t> type;
  std::vector<int32_t> cols;  // N*k
  std::vector<int32_t> y;     // M1_B - M1_A
};

static void build_instance(const Params &P, uint64_t seed, Instance &I) {
  I.N = P.n_common + P.dA + P.dB;
  I.l = P.l;
  I.k = P.k;
  long W = window_of(P);
  Rng rng(seed);
  I.key.resize(I.N);
  for (auto &x : I.key) x = rng.next();
  I.type.assign(I.N, COMMON);
  for (long e = 0; e < P.dA; ++e) I.type[e] = A_ONLY;
  for (long e = P.dA; e < P.dA + P.dB; ++e) I.type[e] = B_ONLY;
  I.cols.resize((size_t)I.N * I.k);
  for (long e = 0; e < I.N; ++e) element_rows(P, W, I.key[e], &I.cols[(size_t)e * I.k]);
  I.y.assign(I.l, 0);
  for (long e = 0; e < I.N; ++e) {
    int s = I.type[e] == B_ONLY ? 1 : I.type[e] == A_ONLY ? -1 : 0;
    if (s)
      for (int t = 0; t < I.k; ++t) I.y[I.cols[(size_t)e * I.k + t]] += s;
  }
}

// ----------------------------------------------------------------------------- Skellam / coding costs
// pmf of Poi(a) - Poi(b) for |z| <= zmax
static std::map<int, double> skellam_pmf(double a, double b, double eps = 1e-18) {
  std::map<int, double> pmf;
  if (a < 1e-12 && b < 1e-12) { pmf[0] = 1.0; return pmf; }
  int zmax = (int)std::ceil(a + b + 12 * std::sqrt(a + b) + 25);
  auto lpois = [](int x, double m) {
    if (m < 1e-300) return x == 0 ? 0.0 : -INFINITY;
    return -m + x * std::log(m) - std::lgamma(x + 1.0);
  };
  for (int z = -zmax; z <= zmax; ++z) {
    double p = 0;
    for (int yv = std::max(0, -z); yv <= zmax; ++yv) {
      double lp = lpois(z + yv, a) + lpois(yv, b);
      if (lp > -800) p += std::exp(lp);
    }
    if (p >= eps) pmf[z] = p;
  }
  return pmf;
}

// rANS cost in bits of `data` under the frequency table `pmf`, mirroring RansWrapper:
// 15-bit cumulative quantization, symbols with zero quantized frequency (or missing) escape with
// (32-bit index + symbol bytes); "default mode" when only one symbol survives.
static double rans_bits(const std::map<int, double> &pmf, const std::vector<int32_t> &data,
                        int symbol_bytes) {
  const double scale = 32768.0;
  double total = 0;
  for (auto &kv : pmf) total += kv.second;
  std::unordered_map<int, uint32_t> freq;
  double cum = 0;
  uint32_t prev = 0;
  int nonzero = 0;
  for (auto &kv : pmf) {
    cum += kv.second;
    uint32_t c = (uint32_t)(cum * scale / total);
    freq[kv.first] = c - prev;
    nonzero += (c - prev) > 0;
    prev = c;
  }
  double bits = 1;  // default-mode marker
  const double escape = 32 + 8.0 * symbol_bytes;
  if (nonzero == 1) {
    bits += 8.0 * symbol_bytes;
    for (int v : data) {
      auto it = freq.find(v);
      if (it == freq.end() || it->second == 0) bits += escape;
    }
    return bits;
  }
  bits += 32;  // flush
  for (int v : data) {
    auto it = freq.find(v);
    if (it == freq.end() || it->second == 0) bits += escape;
    else bits += -std::log2(it->second / scale);
  }
  return bits;
}

static void mean_var(const std::vector<int32_t> &v, size_t lo, size_t hi, double &m, double &s) {
  double n = (double)(hi - lo), sum = 0, sq = 0;
  for (size_t i = lo; i < hi; ++i) { sum += v[i]; sq += (double)v[i] * v[i]; }
  m = sum / n;
  s = n > 1 ? (sq - n * m * m) / (n - 1) : 0;
}

// two_party.cpp: residual message = rANS under moment-fitted Skellam + 64 bits of parameters
static double residual_bits_global(const std::vector<int32_t> &r) {
  double m, s;
  mean_var(r, 0, r.size(), m, s);
  double mu1 = std::max(0.0, (m + s) / 2), mu2 = std::max(0.0, (s - m) / 2);
  return rans_bits(skellam_pmf(mu1, mu2), r, 2) + 64;
}

// Block-adaptive variant: one Skellam fit per block of rows (1-bit "all zero" flag per block,
// 16 bits of quantized parameters otherwise), one rANS stream.
static double residual_bits_blocks(const std::vector<int32_t> &r, long block) {
  double bits = 32;
  for (size_t lo = 0; lo < r.size(); lo += block) {
    size_t hi = std::min(r.size(), lo + (size_t)block);
    bool zero = true;
    for (size_t i = lo; i < hi && zero; ++i) zero = r[i] == 0;
    bits += 1;
    if (zero) continue;
    double m, s;
    mean_var(r, lo, hi, m, s);
    double mu1 = std::max(0.0, (m + s) / 2), mu2 = std::max(0.0, (s - m) / 2);
    std::vector<int32_t> part(r.begin() + lo, r.begin() + hi);
    bits += 16 + rans_bits(skellam_pmf(mu1, mu2), part, 2) - 33;  // 2 x 8-bit log-quantized rates; stream overhead counted once
  }
  return bits;
}

// First message, ported from doro_parameter_tuning (bin/two_party.cpp) for the Skellam
// difference Poi(mu_b) - Poi(mu_a); returns bits per counter.
static double binom_quantile_up(int n, double p, double q) {
  double lp = n * std::log1p(-p), cdf = 0;
  for (int c = 0; c <= n; ++c) {
    if (c > 0) lp += std::log((double)(n - c + 1) / c) + std::log(p) - std::log1p(-p);
    cdf += std::exp(lp);
    if (cdf >= q - 1e-12) return c;
  }
  return n;
}
static double first_message_bits_per_counter(double mu_b, double mu_a) {
  auto pmf = skellam_pmf(mu_b, mu_a);
  int k2 = pmf.begin()->first - 1, k1 = pmf.rbegin()->first + 1;
  auto cdf = [&](int x) {  // P(X < x)
    double s = 0;
    for (auto &kv : pmf) {
      if (kv.first >= x) break;
      s += kv.second;
    }
    return std::min(1.0, s);
  };
  int mean = (int)std::floor(mu_b - mu_a);
  const double dmin = 1e-9, dmax = 0.1, dfinal = 1e-6, bch_block_err = 0.01;
  double best = INFINITY;
  for (int ub = mean + 1; ub <= k1; ++ub) {
    for (int lb = mean; lb > k2; --lb) {
      double diff_err = 1.0 - cdf(ub) + cdf(lb);
      if (diff_err < dmin || diff_err > dmax) continue;
      int interval = ub - lb, cur_lb = lb, cur_ub = ub;
      double cost = std::log2((double)(ub - lb));
      bool ok = true;
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
          int cap = std::max(2, (int)std::ceil(binom_quantile_up(len, diff_err, 1.0 - bch_block_err)));
          double chk = (double)order * cap;
          if (len <= chk) continue;
          rate = std::min(rate, chk / (len - chk));
        }
        if (!std::isfinite(rate)) { ok = false; break; }
        cost += rate;
        interval *= 2;
        diff_err = 1.0 - pin;
        cur_lb = midpoint - interval;
        cur_ub = midpoint + interval;
      }
      if (ok && diff_err <= dfinal) best = std::min(best, cost);
    }
  }
  return best;
}

static double H_skellam(double a, double b) {
  double H = 0;
  for (auto &kv : skellam_pmf(a, b, 1e-15)) H -= kv.second * std::log2(kv.second);
  return H;
}

// ----------------------------------------------------------------------------- SMF signature lengths
// Ported from signature_length / loss_func / finger_l_size / convex_argmax (lib/probability.hpp).
static double bin_entropy(double p) {
  if (p <= 0 || p >= 1) return 0;
  return -p * std::log2(p) - (1 - p) * std::log2(1 - p);
}
static int log2_ceil(double x) { return (int)std::ceil(std::log2(x)); }
static int finger_l_size(double alpha, double beta, double eps, double b) {
  if (alpha < 1e-30 || alpha > 1e30) return 0;
  return (int)std::ceil(-std::log2(eps) + std::log2(alpha) + std::log2(beta) + std::log2(b));
}
static double loss_func(double alpha, double beta, double b, double eps) {
  double p1 = 1 - std::exp(-beta);
  double e1 = bin_entropy(p1);
  double p2 = 1 - std::exp(-alpha * beta);
  double e2 = bin_entropy(p2 * p1);
  double recon = beta * (alpha * p1 + p2) * (finger_l_size(alpha, beta, eps, b) + log2_ceil(b / beta));
  return (e1 + e2 + recon) / beta;
}
template <typename F>
static double convex_argmax(F f, double lb, double ub, double epsilon) {
  double mid = (lb + ub) / 2, fmid = f(mid), flb = f(lb), fub = f(ub);
  while (fmid < flb && ub - lb > epsilon) { ub = mid; mid = (lb + ub) / 2; fmid = f(mid); }
  while (fmid < fub && ub - lb > epsilon) { lb = mid; mid = (lb + ub) / 2; fmid = f(mid); }
  while (ub - lb > epsilon) {
    bool right = ub - mid > mid - lb;
    double &farther = right ? ub : lb, &closer = right ? lb : ub;
    double nm = (mid + farther) / 2, fn = f(nm);
    if (fn > fmid) { closer = mid; mid = nm; fmid = fn; }
    else farther = nm;
  }
  return mid;
}
static std::pair<int, int> signature_length(double amb, double bma, double failure_rate) {
  if (bma < 1 || amb < 1) return {1, 0};
  bma *= 1.25;
  double alpha = amb / bma;
  auto lam = [&](double beta) { return -loss_func(alpha, 1.0 / beta, bma, failure_rate); };
  double p = convex_argmax(lam, 2.0, 50.0, 0.01);
  int s = (int)std::ceil(bma * p);
  return {s, finger_l_size(alpha, 1.0 / p, failure_rate, s)};
}

// ----------------------------------------------------------------------------- party (shared parts)
struct Shared {
  std::vector<int32_t> r;  // residual counters (the shared "doro" code)
  long peels = 0;          // DoroCode::num_peels, shared by both decoders
};

struct Party {
  const Instance *I;
  Shared *S;
  int sign;                         // +1 Bob (values 0/1), -1 Alice (values -1/0)
  long d_own;                       // |own unique| (prior)
  std::vector<int32_t> elems;       // local -> global
  std::vector<int64_t> row_start;   // CSR: row -> edge slots (local*k + t)
  std::vector<int64_t> row_edges;
  std::vector<int8_t> val;          // current decoded value per local (0 or sign)
  // SMF state
  uint32_t finger_s = 1;
  int finger_l = 0;
  uint64_t fseed = 0, rseed = 0;
  std::vector<uint8_t> other_fingers;            // fingerprints_ (received)
  std::vector<std::vector<int32_t>> my_fingers;  // my_fingerprints_: finger -> locals
  std::vector<uint8_t> finger_sent;              // fingers already sent
  std::vector<int32_t> new_elements;
  std::vector<uint8_t> is_new;
  std::vector<std::pair<uint32_t, uint64_t>> unresolved;
  std::vector<int32_t> unresolved_ids;  // local ids
  std::vector<uint8_t> blocked_forever;  // reverted common hallucinations
  bool resolving = false;

  int k() const { return I->k; }
  int32_t row_of(long local, int t) const { return I->cols[(size_t)elems[local] * I->k + t]; }
  uint32_t finger(long local) const {
    return (uint32_t)(((__uint128_t)mix(I->key[elems[local]] ^ fseed) * finger_s) >> 64);
  }
  uint64_t resolve_hash(long local) const {
    uint64_t h = mix(I->key[elems[local]] ^ rseed);
    return finger_l >= 64 ? h : (h & ((1ULL << finger_l) - 1));
  }
  bool blocked(long local) const {
    return blocked_forever[local] || (!resolving && finger_s > 1 && other_fingers[finger(local)]);
  }

  void init(const Instance &inst, Shared &sh, int sgn, long d, uint32_t fs, int fl, uint64_t fsd,
            uint64_t rsd) {
    I = &inst; S = &sh; sign = sgn; d_own = d;
    finger_s = fs; finger_l = fl; fseed = fsd; rseed = rsd;
    elems.clear();
    for (long e = 0; e < inst.N; ++e) {
      bool mine = sgn > 0 ? inst.type[e] != A_ONLY : inst.type[e] != B_ONLY;
      if (mine) elems.push_back((int32_t)e);
    }
    const long n = (long)elems.size();
    row_start.assign(inst.l + 1, 0);
    for (long v = 0; v < n; ++v)
      for (int t = 0; t < inst.k; ++t) row_start[row_of(v, t) + 1]++;
    for (long j = 0; j < inst.l; ++j) row_start[j + 1] += row_start[j];
    row_edges.resize(row_start.back());
    std::vector<int64_t> fill(row_start.begin(), row_start.end() - 1);
    for (long v = 0; v < n; ++v)
      for (int t = 0; t < inst.k; ++t) row_edges[fill[row_of(v, t)]++] = (int64_t)v * inst.k + t;
    val.assign(n, 0);
    other_fingers.assign(finger_s, 0);
    my_fingers.assign(finger_s, {});
    finger_sent.assign(finger_s, 0);
    is_new.assign(n, 0);
    blocked_forever.assign(n, 0);
    new_elements.clear();
    unresolved.clear();
    unresolved_ids.clear();
    resolving = false;
  }

  // DoroCode::peel(element, delta) on the shared code + result bookkeeping
  void peel(long local, int delta) {
    S->peels++;
    for (int t = 0; t < k(); ++t) S->r[row_of(local, t)] -= delta;
    val[local] = (int8_t)(val[local] + delta);
  }
  void mark_new(long local) {
    if (!is_new[local]) { is_new[local] = 1; new_elements.push_back((int32_t)local); }
  }
  void clear_new() {
    for (int32_t v : new_elements) is_new[v] = 0;
    new_elements.clear();
  }
  bool has_new_elements() const { return !new_elements.empty(); }

  // update_fingerprints(): returns the delta bit vector
  std::vector<int32_t> update_fingerprints() {
    std::vector<int32_t> delta(finger_s, 0);
    for (long v = 0; v < (long)val.size(); ++v) {
      if (!val[v]) continue;
      uint32_t f = finger(v);
      if (!finger_sent[f]) { delta[f] = 1; finger_sent[f] = 1; }
      auto &lst = my_fingers[f];
      if (std::find(lst.begin(), lst.end(), (int32_t)v) == lst.end()) lst.push_back((int32_t)v);
    }
    return delta;
  }
  void load_fingerprints(const std::vector<int32_t> &delta) {
    for (uint32_t f = 0; f < finger_s; ++f)
      if (delta[f]) other_fingers[f] = 1;
  }
  void add_unresolved_elements() {
    for (int32_t v : new_elements)
      if (val[v] && finger_s > 1 && other_fingers[finger(v)]) {
        unresolved.push_back({finger(v), resolve_hash(v)});
        unresolved_ids.push_back(v);
      }
  }
  // this = other_decoder in two_party.cpp; `dec` holds the unresolved elements
  int resolve_collision(Party &dec) {
    int collisions = 0;
    for (size_t i = 0; i < dec.unresolved.size(); ++i) {
      auto [f1, f2] = dec.unresolved[i];
      for (int32_t mine : my_fingers[f1]) {
        if (resolve_hash(mine) != f2) continue;
        ++collisions;
        int32_t theirs = dec.unresolved_ids[i];
        if (val[mine]) { int c = val[mine]; peel(mine, -c); }
        if (dec.val[theirs]) { int c = dec.val[theirs]; dec.peel(theirs, -c); }
        blocked_forever[mine] = 1;
        dec.blocked_forever[theirs] = 1;
        on_revert(mine);
        dec.on_revert(theirs);
        break;
      }
    }
    dec.unresolved.clear();
    dec.unresolved_ids.clear();
    return collisions;
  }
  virtual void on_revert(long) {}
  virtual void enter_resolving() { resolving = true; }
  virtual long decode() = 0;
  virtual ~Party() {}
};

// ----------------------------------------------------------------------------- MP (DoroDecoder port)
struct MPParty : Party {
  double lb, ub;
  bool l2 = true;
  long max_peels = 0;
  std::vector<int32_t> str;
  struct Entry {
    int32_t prio;
    uint64_t key;
    int32_t idx;
    bool operator<(const Entry &o) const { return prio != o.prio ? prio < o.prio : key < o.key; }
  };
  std::priority_queue<Entry> pq;
  std::vector<int32_t> acc, touched;
  std::vector<uint8_t> is_touched;
  std::set<std::pair<int32_t, int32_t>> pairs_5in3;

  void setup(long max_peels_) {
    lb = sign > 0 ? 0 : -1;
    ub = sign > 0 ? 1 : 0;
    max_peels = max_peels_;
    str.assign(elems.size(), 0);
    acc.assign(elems.size(), 0);
    is_touched.assign(elems.size(), 0);
  }
  int32_t sense(long v) const {
    int32_t s = 0;
    for (int t = 0; t < k(); ++t) s += S->r[row_of(v, t)];
    return s;
  }
  int32_t get_median(long v) const {
    int32_t sig[64];
    for (int t = 0; t < k(); ++t) sig[t] = S->r[row_of(v, t)];
    int32_t median = 0;
    int32_t *mid = sig + k() / 2;
    std::nth_element(sig, mid, sig + k());
    if (k() % 2 == 0) {
      int32_t *mid2 = mid - 1;
      std::nth_element(sig, mid2, mid);
      if ((long)(*mid) * (*mid2) < 0) median = 0;
      else if (std::abs(*mid) > std::abs(*mid2)) median = *mid2;
    } else median = *mid;
    return median;
  }
  int32_t sense_l1(long v) const {
    int32_t med = get_median(v);
    double o = 0, nw = 0;
    for (int t = 0; t < k(); ++t) {
      int32_t x = S->r[row_of(v, t)];
      o += std::abs(x);
      nw += std::abs(x - med);
    }
    return (int32_t)(o - nw);
  }
  int32_t strength_to_delta(long v, int32_t strength) const {
    double normalized = l2 ? (double)strength / k() : (double)get_median(v);
    double cur = val[v], nv = std::max(lb, std::min(ub, cur + normalized)), delta = nv - cur;
    if (std::fabs(delta) < 0.50001) return 0;
    return (int32_t)std::lround(delta);
  }
  void add_pq(long v, int32_t strength, int32_t delta) {
    str[v] = strength;
    bool fp_ok = finger_s <= 1 || !other_fingers[finger(v)];
    if ((resolving || fp_ok) && delta != 0)
      pq.push({std::abs(strength), I->key[elems[v]], (int32_t)v});
  }
  bool pq_empty() {
    while (!pq.empty()) {
      if (pq.top().prio == std::abs(str[pq.top().idx])) return false;
      pq.pop();
    }
    return true;
  }
  void scan() {
    for (long v = 0; v < (long)elems.size(); ++v) {
      int32_t s = l2 ? sense(v) : sense_l1(v);
      add_pq(v, s, strength_to_delta(v, s));
    }
  }
  void notify(int32_t row, int32_t change) {
    for (int64_t p = row_start[row]; p < row_start[row + 1]; ++p) {
      int32_t nb = (int32_t)(row_edges[p] / k());
      if (!is_touched[nb]) { is_touched[nb] = 1; touched.push_back(nb); acc[nb] = 0; }
      acc[nb] += change;
    }
  }
  void update_neighbors(long element, bool append_53) {
    for (int32_t nb : touched) {
      if (append_53 && k() == 5 && nb != element && std::abs(acc[nb]) >= 3)
        pairs_5in3.insert({std::min<int32_t>(element, nb), std::max<int32_t>(element, nb)});
      int32_t ns = l2 ? str[nb] + acc[nb] : sense_l1(nb);
      add_pq(nb, ns, strength_to_delta(nb, ns));
      is_touched[nb] = 0;
    }
    touched.clear();
  }
  bool peel_until_empty(long hard_cap) {
    while (!pq_empty()) {
      Entry t = pq.top();
      int32_t delta = strength_to_delta(t.idx, str[t.idx]);
      if (delta == 0) { pq.pop(); continue; }
      mark_new(t.idx);
      peel(t.idx, delta);
      for (int u = 0; u < k(); ++u) notify(row_of(t.idx, u), -delta);
      update_neighbors(t.idx, true);
      if (S->peels > hard_cap) return false;
    }
    return true;
  }
  void resolve_5in3() {
    for (auto [e1, e2] : pairs_5in3) {
      int v1 = val[e1], v2 = val[e2];
      if (!v1 && !v2) continue;
      int vv = v1 ? v1 : v2;
      int32_t cur = v1 ? e1 : e2, oth = v1 ? e2 : e1;
      int32_t strength = (sense(cur) - sense(oth)) * (-vv);
      if (strength >= 3) {
        peel(cur, -vv);
        peel(oth, vv);
        mark_new(oth);
      }
    }
  }
  long decode() override {
    long before = S->peels;
    while (!pq.empty()) pq.pop();
    clear_new();
    scan();
    long hard_cap = before + 60L * (long)elems.size() / 10 + 100000;
    for (int rnd = 0; rnd < 100; ++rnd) {
      if (!peel_until_empty(hard_cap)) break;
      if (resolving && l2) { l2 = false; scan(); }
      if (max_peels > 0 && S->peels >= max_peels) break;
      if (pq_empty()) break;
    }
    if (resolving && !pairs_5in3.empty()) resolve_5in3();
    unresolved.clear();
    unresolved_ids.clear();
    if (resolving) add_unresolved_elements();
    while (!pq.empty()) pq.pop();
    return S->peels - before;
  }
};

// ----------------------------------------------------------------------------- BP party
struct BPParty : Party {
  const Params *P = nullptr;
  std::vector<float> q2c, c2v;
  std::vector<double> tot;
  std::vector<double> noise;  // per row Poisson noise rate (other party's undecoded uniques)
  long block = 1024;
  bool started = false;
  long iters_total = 0;

  void setup(const Params &prm, long block_rows) {
    P = &prm;
    block = block_rows;
    size_t ne = elems.size() * (size_t)k();
    double pi = std::min(0.5, std::max(1e-9, (double)d_own / elems.size()));
    q2c.assign(ne, (float)pi);
    c2v.assign(ne, 0.f);
    tot.assign(elems.size(), 0.0);
    noise.assign(I->l, 0.0);
  }

  // w_j = sign * (r_j + sum of my own committed contributions) = own positives - other's noise
  void view(std::vector<int32_t> &w) const {
    w = S->r;
    for (long v = 0; v < (long)val.size(); ++v)
      if (val[v])
        for (int t = 0; t < k(); ++t) w[row_of(v, t)] += val[v];
    if (sign < 0)
      for (auto &x : w) x = -x;
  }

  long decode() override {
    std::vector<int32_t> w;
    view(w);
    // per-block Skellam moment fit: E[w] = muS - muN, Var[w] = muS + muN
    for (size_t lo = 0; lo < w.size(); lo += block) {
      size_t hi = std::min(w.size(), lo + (size_t)block);
      double m, s;
      mean_var(w, lo, hi, m, s);
      double muN = std::max(1e-4, (s - m) / 2);
      for (size_t j = lo; j < hi; ++j) noise[j] = muN;
    }
    const long n = (long)elems.size();
    const int K = k();
    double pi = std::min(0.5, std::max(1e-9, (double)d_own / n));
    const double prior = std::log(pi / (1 - pi)), CL = 30.0;
    // blocked elements (other party's fingerprints during avoidance, reverted collisions) keep
    // their current state, as in DoroDecoder where they are never pushed to the queue
    std::vector<uint8_t> frozen(n, 0);
    for (long v = 0; v < n; ++v) frozen[v] = blocked(v) ? (val[v] ? 2 : 1) : 0;
    for (long v = 0; v < n; ++v)
      if (frozen[v])
        for (int t = 0; t < K; ++t) q2c[(size_t)v * K + t] = frozen[v] == 2 ? 1.f : 0.f;
    std::vector<double> pre, suf, pois;
    std::vector<uint8_t> hard(n, 0);
    int stable = 0, it = 0;
    for (; it < P->bp_iters; ++it) {
      for (long j = 0; j < I->l; ++j) {
        int64_t s0 = row_start[j], s1 = row_start[j + 1];
        long deg = (long)(s1 - s0);
        if (deg == 0) continue;
        int wj = w[j];
        double mu = noise[j];
        int nmax = (int)std::ceil(mu + 6 * std::sqrt(mu) + 4);
        nmax = std::max(nmax, -wj + 4);
        if (wj > deg) wj = (int)deg;  // cannot happen without errors; keep the DP bounded
        // pre over v in [-nmax, wj]  (index v + nmax), suf over u in [0, wj + nmax]
        int R = wj + nmax + 1;
        if (R <= 0) { for (int64_t p = s0; p < s1; ++p) c2v[row_edges[p]] = 0.f; continue; }
        pois.assign(nmax + 1, 0.0);
        for (int x = 0; x <= nmax; ++x)
          pois[x] = std::exp(-mu + (x ? x * std::log(std::max(mu, 1e-300)) : 0.0) - std::lgamma(x + 1.0));
        pre.assign((size_t)(deg + 1) * R, 0.0);
        suf.assign((size_t)(deg + 1) * R, 0.0);
        for (int x = 0; x <= nmax; ++x) {
          int v = -x;  // value -x at index v + nmax
          if (v <= wj) pre[v + nmax] = pois[x];
        }
        for (long t = 0; t < deg; ++t) {
          double q = q2c[row_edges[s0 + t]];
          const double *c = &pre[t * R];
          double *nx = &pre[(t + 1) * R];
          double mx = 0;
          for (int i = 0; i < R; ++i) {
            nx[i] = c[i] * (1 - q) + (i ? c[i - 1] * q : 0.0);
            mx = std::max(mx, nx[i]);
          }
          if (mx > 0) for (int i = 0; i < R; ++i) nx[i] /= mx;
        }
        suf[deg * R] = 1.0;
        for (long t = deg - 1; t >= 0; --t) {
          double q = q2c[row_edges[s0 + t]];
          const double *c = &suf[(t + 1) * R];
          double *nx = &suf[t * R];
          double mx = 0;
          for (int i = 0; i < R; ++i) {
            nx[i] = c[i] * (1 - q) + (i ? c[i - 1] * q : 0.0);
            mx = std::max(mx, nx[i]);
          }
          if (mx > 0) for (int i = 0; i < R; ++i) nx[i] /= mx;
        }
        for (long t = 0; t < deg; ++t) {
          const double *A = &pre[t * R], *B = &suf[(t + 1) * R];
          // P(T_-t = target) = sum_v A[v] B[target - v], v in [-nmax, target], target - v < R
          auto at = [&](int target) {
            double s = 0;
            for (int v = -nmax; v <= target; ++v) {
              int u = target - v;
              if (u >= R) continue;
              s += A[v + nmax] * B[u];
            }
            return s;
          };
          double m0 = at(wj), m1 = at(wj - 1);
          double L;
          if (m1 <= 0 && m0 <= 0) L = 0;
          else if (m1 <= 0) L = -CL;
          else if (m0 <= 0) L = CL;
          else L = std::max(-CL, std::min(CL, std::log(m1) - std::log(m0)));
          c2v[row_edges[s0 + t]] = (float)L;
        }
      }
      double max_change = 0;
      for (long v = 0; v < n; ++v) {
        double s = prior;
        for (int t = 0; t < K; ++t) s += c2v[(size_t)v * K + t];
        max_change = std::max(max_change, std::fabs(std::max(-CL, std::min(CL, s)) - std::max(-CL, std::min(CL, tot[v]))));
        tot[v] = s;
        hard[v] = frozen[v] ? frozen[v] == 2 : s > 0;
        if (frozen[v]) continue;
        for (int t = 0; t < K; ++t) {
          double p = 1.0 / (1.0 + std::exp(-(s - c2v[(size_t)v * K + t])));
          float &q = q2c[(size_t)v * K + t];
          q = (float)(P->bp_damp * q + (1 - P->bp_damp) * p);
        }
      }
      // converged messages, not merely stable hard decisions: with coupling, confidence grows as
      // a wave over many iterations while the hard decisions can stay flat
      if (it > 0 && max_change < P->bp_tol) ++stable; else stable = 0;
      if (it + 1 >= P->bp_min_iters && stable >= P->bp_stable) { ++it; break; }
    }
    iters_total += it;
    // commit the new decision set
    clear_new();
    long changes = 0;
    for (long v = 0; v < n; ++v) {
      // hysteresis: a new commitment needs confidence, an existing one only a positive LLR
      bool keep = frozen[v] ? frozen[v] == 2 : (val[v] ? tot[v] > 0 : tot[v] > P->commit_llr);
      int want = keep ? sign : 0;
      if (want != val[v]) {
        if (want) mark_new(v);
        peel(v, want - val[v]);
        ++changes;
      }
    }
    unresolved.clear();
    unresolved_ids.clear();
    if (resolving) add_unresolved_elements();
    return changes;
  }
};

// ----------------------------------------------------------------------------- protocol
struct Result {
  bool success = false, zero = false;
  int rounds = 0;
  double first_bits = 0, residual_bits = 0, residual_block_bits = 0, residual_best_bits = 0;
  double finger_bits = 0, resolve_bits = 0;
  long bp_iters = 0;
  double secs = 0;
};

static double finger_bits(const std::vector<int32_t> &delta, uint32_t finger_s) {
  double ones = 0;
  for (int32_t b : delta) ones += b;
  if (ones < 1 || finger_s <= 1) return 0;
  double load = ones / finger_s;
  std::map<int, double> pmf{{0, 1.0 - load}, {1, load}};
  return rans_bits(pmf, delta, 1) + 32;
}

static Result run_protocol(const Params &P, const Instance &I, long block_rows) {
  Result R;
  Shared S;
  S.r = I.y;
  long W = window_of(P);
  long bulk = P.matrix == "uniform" ? P.l : P.l - W + 1;
  double mu_b = (double)P.dB * P.k / bulk, mu_a = (double)P.dA * P.k / bulk;
  R.first_bits = P.l * first_message_bits_per_counter(mu_b, mu_a);

  double fail = P.failure_rate / std::max(1L, P.dB);
  auto [fs, fl] = signature_length((double)P.dA, (double)P.dB, fail);
  uint64_t fseed = mix(P.seed ^ 0xF1F1), rseed = mix(P.seed ^ 0x2E2E);

  MPParty mpA, mpB;
  BPParty bpA, bpB;
  Party *alis, *bela;
  if (P.decoder == "mp") {
    mpA.init(I, S, -1, P.dA, fs, fl, fseed, rseed);
    mpB.init(I, S, +1, P.dB, fs, fl, fseed, rseed);
    long max_peels = (long)(2.5 * (P.dA + P.dB));
    mpA.setup(max_peels);
    mpB.setup(max_peels);
    if (P.dA == 0) mpA.lb = 0;
    if (P.dB == 0) mpB.ub = 0;
    alis = &mpA; bela = &mpB;
  } else {
    bpA.init(I, S, -1, P.dA, fs, fl, fseed, rseed);
    bpB.init(I, S, +1, P.dB, fs, fl, fseed, rseed);
    bpA.setup(P, block_rows);
    bpB.setup(P, block_rows);
    alis = &bpA; bela = &bpB;
  }
  auto residual_zero = [&]() {
    for (int32_t v : S.r)
      if (v) return false;
    return true;
  };
  enum { AVOID, RESOLVE, FINISHED } status = AVOID;
  if (P.resolving_round < 0) { alis->enter_resolving(); bela->enter_resolving(); status = RESOLVE; }
  int comm_rounds = 0, actual = 0, round_in_resolve = 0;
  bool party_is_bob = false;
  std::vector<int32_t> last_sent;
  while (!residual_zero() && comm_rounds < P.max_rounds && status != FINISHED) {
    actual = std::max(actual, comm_rounds + 1);
    party_is_bob = !party_is_bob;
    Party &dec = party_is_bob ? *bela : *alis;
    Party &oth = party_is_bob ? *alis : *bela;
    last_sent = S.r;  // the residual this party received
    dec.decode();
    if (residual_zero() && dec.unresolved.empty()) break;
    // CommonSense codes the whole residual with one fitted Skellam. Alternatives, any of which a
    // sender may pick per message (2 bits to say which): block-adaptive coding of the residual,
    // or of its change since the residual the sender received (both sides hold that one; not
    // available for the first residual message, since Alice never saw y).
    {
      double g = residual_bits_global(S.r);
      double best = g;
      for (long b : {block_rows, 4 * block_rows, 1024L}) {
        best = std::min(best, residual_bits_blocks(S.r, b));
        if (comm_rounds > 0) {
          std::vector<int32_t> diff(S.r.size());
          for (size_t j = 0; j < diff.size(); ++j) diff[j] = S.r[j] - last_sent[j];
          best = std::min(best, residual_bits_blocks(diff, b));
        }
      }
      R.residual_bits += g;
      R.residual_block_bits += residual_bits_blocks(S.r, block_rows);
      R.residual_best_bits += best + 2;
    }
    auto delta = dec.update_fingerprints();
    R.finger_bits += finger_bits(delta, fs);
    oth.load_fingerprints(delta);
    if (!dec.unresolved.empty()) {
      double cost = dec.unresolved.size() * (double)(log2_ceil(fs) + fl);
      int extra = 1;
      size_t nu = dec.unresolved.size();
      int coll = oth.resolve_collision(dec);
      if (coll > 0) { cost += nu; extra = 2; }
      R.resolve_bits += cost;
      actual = std::max(actual, comm_rounds + extra + 1);
    }
    if (!dec.has_new_elements() || comm_rounds == P.resolving_round) {
      if (status == AVOID) {
        status = RESOLVE;
        alis->enter_resolving();
        bela->enter_resolving();
      }
      if (!dec.has_new_elements() && round_in_resolve >= 2) status = FINISHED;
      ++round_in_resolve;
    }
    ++comm_rounds;
    if (P.verbose > 1) {
      long nz = 0;
      for (int32_t v : S.r) nz += v != 0;
      fprintf(stderr, "  round %d (%s): nonzero residual rows %ld, status %d\n", comm_rounds,
              party_is_bob ? "Bob" : "Alice", nz, (int)status);
      if (P.verbose > 2) {
        // fraction of each party's unique elements decoded, by decile of their first row
        for (Party *pp : {alis, bela}) {
          std::vector<long> got(10, 0), all(10, 0), fp(10, 0);
          uint8_t mine = pp->sign > 0 ? B_ONLY : A_ONLY;
          for (long v = 0; v < (long)pp->val.size(); ++v) {
            int32_t lo = (int32_t)P.l;
            for (int t = 0; t < I.k; ++t) lo = std::min(lo, pp->row_of(v, t));
            int b = (int)((long)lo * 10 / P.l);
            if (I.type[pp->elems[v]] == mine) { all[b]++; got[b] += pp->val[v] != 0; }
            else fp[b] += pp->val[v] != 0;
          }
          fprintf(stderr, "    %s decoded by decile:", pp->sign > 0 ? "Bob  " : "Alice");
          for (int b = 0; b < 10; ++b) fprintf(stderr, " %.2f", all[b] ? (double)got[b] / all[b] : 0.0);
          long f = 0;
          for (long x : fp) f += x;
          fprintf(stderr, "   false=%ld\n", f);
        }
      }
    }
  }
  R.rounds = actual;
  R.zero = residual_zero();
  bool ok = true;
  for (long v = 0; v < (long)alis->val.size() && ok; ++v)
    ok = (alis->val[v] != 0) == (I.type[alis->elems[v]] == A_ONLY);
  for (long v = 0; v < (long)bela->val.size() && ok; ++v)
    ok = (bela->val[v] != 0) == (I.type[bela->elems[v]] == B_ONLY);
  R.success = ok;
  R.bp_iters = bpA.iters_total + bpB.iters_total;
  return R;
}

// ----------------------------------------------------------------------------- main
static double log2_binom(double n, double k) {
  return (std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1)) / std::log(2.0);
}

static Params parse_args(int argc, char **argv) {
  Params P;
  double ratio = -1;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto eq = a.find('=');
    if (eq == std::string::npos) { fprintf(stderr, "bad arg %s\n", a.c_str()); exit(1); }
    std::string key = a.substr(0, eq), v = a.substr(eq + 1);
    if (key == "n") P.n_common = std::stol(v);
    else if (key == "dA") P.dA = std::stol(v);
    else if (key == "dB") P.dB = std::stol(v);
    else if (key == "k") P.k = std::stoi(v);
    else if (key == "l") P.l = std::stol(v);
    else if (key == "ratio") ratio = std::stod(v);  // l = ratio * (dA + dB)
    else if (key == "matrix") P.matrix = v;
    else if (key == "W") P.W = std::stol(v);
    else if (key == "wfrac") P.wfrac = std::stod(v);
    else if (key == "offsets") P.offsets = v;
    else if (key == "decoder") P.decoder = v;
    else if (key == "trials") P.trials = std::stoi(v);
    else if (key == "threads") P.threads = std::stoi(v);
    else if (key == "seed") P.seed = std::stoull(v);
    else if (key == "max_rounds") P.max_rounds = std::stoi(v);
    else if (key == "resolving_round") P.resolving_round = std::stoi(v);
    else if (key == "bp_iters") P.bp_iters = std::stoi(v);
    else if (key == "bp_damp") P.bp_damp = std::stod(v);
    else if (key == "bp_tol") P.bp_tol = std::stod(v);
    else if (key == "commit_llr") P.commit_llr = std::stod(v);
    else if (key == "block") P.block_rows = std::stol(v);
    else if (key == "verbose") P.verbose = std::stoi(v);
    else { fprintf(stderr, "unknown key %s\n", key.c_str()); exit(1); }
  }
  if (ratio > 0) P.l = std::lround(ratio * (P.dA + P.dB));
  if (P.matrix == "uniform") P.offsets = "-";
  return P;
}

int main(int argc, char **argv) {
  Params P = parse_args(argc, argv);
  long W = window_of(P);
  long block_rows = P.block_rows > 0 ? P.block_rows
                                     : (P.matrix == "uniform" ? std::max(256L, P.l / 64) : std::max(64L, W / 4));
  std::atomic<int> next{0};
  std::mutex mu;
  int n_ok = 0, n_zero = 0, n_silent = 0;
  double s_total = 0, s_total_blk = 0, s_rounds = 0, s_secs = 0, s_iters = 0;
  double s_first = 0, s_res = 0, s_resb = 0, s_resbest = 0, s_fing = 0, s_resolve = 0, s_total_best = 0;
  auto worker = [&]() {
    Instance I;
    while (true) {
      int t = next.fetch_add(1);
      if (t >= P.trials) break;
      uint64_t s = P.seed * 1000003ULL + (uint64_t)t * 7919ULL + 17;
      build_instance(P, splitmix64(s), I);
      auto t0 = std::chrono::steady_clock::now();
      Result R = run_protocol(P, I, block_rows);
      R.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      std::lock_guard<std::mutex> g(mu);
      n_ok += R.success;
      n_zero += R.zero;
      n_silent += R.zero && !R.success;
      if (R.success) {
        double tot = R.first_bits + R.residual_bits + R.finger_bits + R.resolve_bits;
        double totb = R.first_bits + R.residual_block_bits + R.finger_bits + R.resolve_bits;
        double totbest = R.first_bits + R.residual_best_bits + R.finger_bits + R.resolve_bits;
        s_total += tot; s_total_blk += totb; s_total_best += totbest; s_resbest += R.residual_best_bits;
        s_first += R.first_bits; s_res += R.residual_bits; s_resb += R.residual_block_bits;
        s_fing += R.finger_bits; s_resolve += R.resolve_bits; s_rounds += R.rounds;
      }
      s_secs += R.secs;
      s_iters += R.bp_iters;
      if (P.verbose)
        fprintf(stderr, "trial %d ok=%d zero=%d rounds=%d first=%.0f res=%.0f resblk=%.0f fing=%.0f resolve=%.0f %.1fs\n",
                t, R.success, R.zero, R.rounds, R.first_bits, R.residual_bits, R.residual_block_bits,
                R.finger_bits, R.resolve_bits, R.secs);
    }
  };
  std::vector<std::thread> th;
  for (int i = 0; i < std::max(1, std::min(P.threads, P.trials)); ++i) th.emplace_back(worker);
  for (auto &t : th) t.join();

  double lb = log2_binom((double)(P.n_common + P.dA), (double)P.dA) +
              log2_binom((double)(P.n_common + P.dB), (double)P.dB);
  // ideal first message: per-row entropy of the Skellam difference
  std::vector<double> mass = row_edge_mass(P);
  double first_ideal = 0;
  std::map<long, std::pair<double, long>> cls;
  for (double m : mass) { auto &c = cls[std::llround(m * 1e7)]; c.first = m; c.second++; }
  for (auto &kv : cls) first_ideal += kv.second.second * H_skellam(kv.second.first * P.dB, kv.second.first * P.dA);
  double ok = std::max(1, n_ok);
  printf("{\"n_common\":%ld,\"dA\":%ld,\"dB\":%ld,\"k\":%d,\"l\":%ld,\"ratio\":%.4f,\"matrix\":\"%s\",\"W\":%ld,"
         "\"wfrac\":%.4f,\"offsets\":\"%s\",\"decoder\":\"%s\",\"trials\":%d,\"success\":%d,\"zero_residual\":%d,"
         "\"silent_errors\":%d,\"avg_rounds\":%.2f,\"avg_total_bits\":%.0f,\"avg_total_KB\":%.2f,"
         "\"avg_total_blocks_KB\":%.2f,\"avg_first_KB\":%.2f,\"avg_residual_KB\":%.2f,"
         "\"avg_residual_blocks_KB\":%.2f,\"avg_finger_KB\":%.2f,\"avg_resolve_KB\":%.2f,"
         "\"avg_total_best_KB\":%.2f,\"avg_residual_best_KB\":%.2f,\"first_ideal_KB\":%.2f,\"lower_bound_KB\":%.2f,\"total_over_lb\":%.3f,\"total_blocks_over_lb\":%.3f,\"total_best_over_lb\":%.3f,"
         "\"block_rows\":%ld,\"avg_bp_iters\":%.0f,\"avg_secs\":%.2f}\n",
         P.n_common, P.dA, P.dB, P.k, P.l, (double)P.l / (P.dA + P.dB), P.matrix.c_str(),
         P.matrix == "uniform" ? 0L : W, P.matrix == "uniform" ? 0.0 : (double)W / P.l, P.offsets.c_str(),
         P.decoder.c_str(), P.trials, n_ok, n_zero, n_silent, s_rounds / ok, s_total / ok, s_total / ok / 8000,
         s_total_blk / ok / 8000, s_first / ok / 8000, s_res / ok / 8000, s_resb / ok / 8000, s_fing / ok / 8000,
         s_resolve / ok / 8000, s_total_best / ok / 8000, s_resbest / ok / 8000, first_ideal / 8000, lb / 8000, (s_total / ok) / lb, (s_total_blk / ok) / lb, (s_total_best / ok) / lb,
         block_rows, s_iters / std::max(1, P.trials), s_secs / std::max(1, P.trials));
  return 0;
}
