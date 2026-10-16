// End-to-end check of DoroCode::enable_spatial_coupling + ScBpDecoder on the real headers.
// Alice encodes A, Bob encodes B (A subset of B), Bob decodes B \ A from the counter difference.
//
// Build (needs wyhash headers and libONIAK's include dirs; Eigen is only used for type aliases):
//   clang++ -std=c++2b -O2 -I../libONIAK -I<wyhash parent> -I<eigen parent> test_sc_decoder.cpp
// Usage: ./test_sc_decoder [n_common d ratio window_fraction tle trials]

#include "../lib/doro.hpp"
#include "../lib/sc_decoder.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <unordered_set>

using Code = Doro::DoroCode<uint64_t, int16_t>;

int main(int argc, char** argv) {
  long n_common = argc > 1 ? std::atol(argv[1]) : 100000;
  long d = argc > 2 ? std::atol(argv[2]) : 1000;
  double ratio = argc > 3 ? std::atof(argv[3]) : 5.5;
  double wfrac = argc > 4 ? std::atof(argv[4]) : 0.1;  // 0 = uniform matrix
  bool tle = argc > 5 ? std::atoi(argv[5]) != 0 : true;
  int trials = argc > 6 ? std::atoi(argv[6]) : 5;
  const int k = 8;
  const int l = static_cast<int>(ratio * d);
  int ok = 0, silent = 0;
  for (int trial = 0; trial < trials; ++trial) {
    std::mt19937 rng(1234 + trial);
    std::mt19937_64 keys(99 + trial);
    std::unordered_set<uint64_t> setA, setB, truth;
    while (static_cast<long>(setB.size()) < n_common + d) {
      uint64_t x = keys();
      if (!setB.insert(x).second) continue;
      if (static_cast<long>(setA.size()) < n_common) setA.insert(x);
      else truth.insert(x);
    }
    Code bob(l, k, /*is_cbf*/ true, rng, 0, 0, /*track_value*/ false);
    if (wfrac > 0) bob.enable_spatial_coupling(static_cast<int>(wfrac * l), tle, rng());
    Code alice = bob.copy();
    Code::SparseVector a_vals, b_vals;
    for (auto x : setA) a_vals[x] = 1;
    for (auto x : setB) b_vals[x] = 1;
    alice.encode(a_vals);
    bob.encode(b_vals);
    // window sanity check: every column's counters lie inside one window
    if (wfrac > 0) {
      for (auto x : truth) {
        auto h = bob.hash_all(x);
        int lo = l, hi = -1;
        for (auto [idx, sign] : h) { lo = std::min(lo, idx); hi = std::max(hi, idx); }
        if (hi - lo >= bob.sc_window() || static_cast<int>(h.size()) != k) { std::puts("WINDOW VIOLATION"); return 1; }
      }
    }
    for (int j = 0; j < l; ++j) bob.code()[j] -= alice.code()[j];  // M 1_{B \ A}
    auto t0 = std::chrono::steady_clock::now();
    Doro::ScBpDecoder<Code> decoder(bob, setB);
    bool zero = decoder.decode();
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    bool correct = zero && bob.empty() && decoder.result().size() == truth.size();
    if (correct)
      for (auto& [x, v] : decoder.result()) correct &= truth.contains(x) && v == 1;
    ok += correct;
    silent += zero && !correct;
    std::printf("trial %d: zero_residual=%d correct=%d undecided=%ld bp_iters=%ld %.2fs\n", trial, zero,
                correct, decoder.num_undecided(), decoder.bp_iterations(), secs);
  }
  std::printf("n=%ld d=%ld l/d=%.2f wfrac=%.3f tle=%d: %d/%d correct, %d silent errors\n", n_common + d, d, ratio,
              wfrac, tle, ok, trials, silent);
  return ok == trials ? 0 : 2;
}
