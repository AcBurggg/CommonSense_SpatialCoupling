#pragma once
// Exact peeling + damped belief-propagation decoder for unidirectional CommonSense (A subset of B).
//
// Intended to be paired with DoroCode::enable_spatial_coupling(): the coupled matrix only pays off
// with a message-passing decoder (peeling/BP propagate a decoding wave in from the low-density
// boundary rows), whereas the greedy MP decoder in doro_decoder.hpp gets worse on it because
// localized columns are more coherent. Works on the uniform matrix too.
//
// Decoding, given the difference counters y = M 1_{B\A} (code.code()) and Bob's set B:
//   1. Exact peeling. A row with residual 0 forces all its undecided candidates to 0; a row whose
//      residual equals its number of undecided candidates forces them all to 1. Every decision is
//      implied by y, so peeling never errs, whatever the matrix.
//   2. Sum-product BP on what peeling leaves, with exact leave-one-out Poisson-binomial check
//      messages, damping, and BP-guided decimation (fix the most confident variables, re-peel).
// The residual must reach zero; only then are the decoded elements peeled from the code and
// reported. On failure the code is left untouched so the caller can fall back (MP / IBLT).
//
// See sc_study/REPORT.md for the analysis and sc_study/sc_sim.cpp for the simulator version.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Doro {

struct ScBpConfig {
  int bp_iters = 1500;   // BP iterations per decimation round
  int max_rounds = 30;   // decimation rounds (1 = plain BP)
  double damping = 0.5;  // fraction of the old variable-to-check message kept each iteration
  double decimate_fraction = 0.02;
};

template <typename DoroCodeT>
class ScBpDecoder {
public:
  using ArrType = typename DoroCodeT::ArrayT;
  using IndexType = typename DoroCodeT::IndexT;

  ScBpDecoder(DoroCodeT& code, const std::unordered_set<IndexType>& candidates, ScBpConfig config = {})
      : code_(&code), config_(config), k_(code.k()), l_(code.size()) {
    if (!code.is_cbf())
      throw std::invalid_argument("ScBpDecoder requires a counting (CBF) code");
    elements_.assign(candidates.begin(), candidates.end());
    n_ = static_cast<long>(elements_.size());
    cols_.resize(static_cast<size_t>(n_) * k_);
    for (long e = 0; e < n_; ++e) {
      auto hashes = code.hash_all(elements_[e]);
      for (int t = 0; t < k_; ++t) cols_[static_cast<size_t>(e) * k_ + t] = hashes[t].first;
    }
    row_start_.assign(l_ + 1, 0);
    for (int r : cols_) ++row_start_[r + 1];
    for (long j = 0; j < l_; ++j) row_start_[j + 1] += row_start_[j];
    row_elems_.resize(cols_.size());
    std::vector<int64_t> fill(row_start_.begin(), row_start_.end() - 1);
    for (long e = 0; e < n_; ++e)
      for (int t = 0; t < k_; ++t) row_elems_[fill[cols_[static_cast<size_t>(e) * k_ + t]]++] = static_cast<int32_t>(e);
  }

  // Returns true iff a zero residual was reached; the code is then peeled to all zeros.
  bool decode() {
    reset_peeling();
    // y = M 1_{B\A} must satisfy 0 <= y_j <= (#candidates in row j); otherwise B \ A is not
    // representable over B (e.g. A is not a subset of B, or the first-round BCH layer failed).
    for (long j = 0; j < l_; ++j)
      if (r_[j] < 0 || r_[j] > u_[j]) return false;
    long total = 0;
    for (long j = 0; j < l_; ++j) total += r_[j];
    if (total % k_ != 0) return false;
    expected_ones_ = total / k_;  // every column has exactly k ones, so sum(y) = k |B \ A|
    run_peeling();
    for (int round = 0; round < config_.max_rounds && n_unknown_ > 0; ++round) {
      if (run_bp()) break;
      decimate();
      run_peeling();
      for (long j = 0; j < l_; ++j)
        if (r_[j] < 0 || r_[j] > u_[j]) return false;  // a decimation decision was wrong
    }
    if (n_unknown_ > 0) return false;
    for (long j = 0; j < l_; ++j)
      if (r_[j] != 0) return false;
    for (long e = 0; e < n_; ++e) {
      if (state_[e] == 1) {
        code_->peel(elements_[e], 1);
        result_[elements_[e]] = 1;
      }
    }
    return true;
  }

  std::unordered_map<IndexType, ArrType>& result() { return result_; }
  long num_undecided() const { return n_unknown_; }
  long bp_iterations() const { return bp_iterations_; }

private:
  // ------------------------------------------------------------------ exact peeling
  void reset_peeling() {
    state_.assign(n_, -1);
    r_.assign(l_, 0);
    for (long j = 0; j < l_; ++j) r_[j] = code_->code()[j];
    u_.assign(l_, 0);
    for (long j = 0; j < l_; ++j) u_[j] = static_cast<int32_t>(row_start_[j + 1] - row_start_[j]);
    in_queue_.assign(l_, 0);
    queue_.clear();
    n_unknown_ = n_;
    n_ones_ = 0;
    result_.clear();
  }

  void maybe_push(int32_t j) {
    if (!in_queue_[j] && u_[j] > 0 && (r_[j] == 0 || r_[j] == u_[j])) {
      in_queue_[j] = 1;
      queue_.push_back(j);
    }
  }

  void set_value(long e, int8_t v) {
    state_[e] = v;
    --n_unknown_;
    if (v) ++n_ones_;
    for (int t = 0; t < k_; ++t) {
      int32_t j = cols_[static_cast<size_t>(e) * k_ + t];
      --u_[j];
      r_[j] -= v;
      maybe_push(j);
    }
  }

  void run_peeling() {
    for (long j = 0; j < l_; ++j) maybe_push(static_cast<int32_t>(j));
    while (!queue_.empty()) {
      int32_t j = queue_.back();
      queue_.pop_back();
      in_queue_[j] = 0;
      if (u_[j] == 0) continue;
      int8_t v;
      if (r_[j] == 0) v = 0;
      else if (r_[j] == u_[j]) v = 1;
      else continue;
      for (int64_t p = row_start_[j]; p < row_start_[j + 1]; ++p)
        if (state_[row_elems_[p]] < 0) set_value(row_elems_[p], v);
    }
  }

  // ------------------------------------------------------------------ belief propagation
  // Returns true if the hard decisions satisfy every residual (then they are committed).
  bool run_bp() {
    vars_.clear();
    for (long e = 0; e < n_; ++e)
      if (state_[e] < 0) vars_.push_back(static_cast<int32_t>(e));
    const long nv = static_cast<long>(vars_.size());
    std::vector<int32_t> rows;
    std::vector<int32_t> row_id(l_, -1);
    for (long j = 0; j < l_; ++j)
      if (u_[j] > 0) {
        row_id[j] = static_cast<int32_t>(rows.size());
        rows.push_back(static_cast<int32_t>(j));
      }
    std::vector<int64_t> rstart(rows.size() + 1, 0);
    for (long v = 0; v < nv; ++v)
      for (int t = 0; t < k_; ++t) ++rstart[row_id[col(vars_[v], t)] + 1];
    for (size_t a = 0; a < rows.size(); ++a) rstart[a + 1] += rstart[a];
    std::vector<int64_t> redge(rstart.back());
    {
      std::vector<int64_t> fill(rstart.begin(), rstart.end() - 1);
      for (long v = 0; v < nv; ++v)
        for (int t = 0; t < k_; ++t) redge[fill[row_id[col(vars_[v], t)]]++] = static_cast<int64_t>(v) * k_ + t;
    }
    const double pi = std::clamp(static_cast<double>(expected_ones_ - n_ones_) / std::max(1L, nv), 1e-9, 0.999);
    const double prior = std::log(pi / (1 - pi));
    const double clamp_llr = 30.0;
    std::vector<float> q2c(static_cast<size_t>(nv) * k_, static_cast<float>(pi));
    std::vector<float> c2v(static_cast<size_t>(nv) * k_, 0.f);
    total_llr_.assign(nv, prior);
    std::vector<double> pre, suf;
    for (int it = 0; it < config_.bp_iters; ++it) {
      ++bp_iterations_;
      for (size_t a = 0; a < rows.size(); ++a) {
        const int32_t rj = r_[rows[a]];
        const int64_t s0 = rstart[a];
        const long deg = static_cast<long>(rstart[a + 1] - s0);
        const int S = rj + 1;
        pre.assign(static_cast<size_t>(deg + 1) * S, 0.0);
        suf.assign(static_cast<size_t>(deg + 1) * S, 0.0);
        pre[0] = 1.0;
        for (long t = 0; t < deg; ++t) convolve_step(&pre[t * S], &pre[(t + 1) * S], S, q2c[redge[s0 + t]]);
        suf[deg * S] = 1.0;
        for (long t = deg - 1; t >= 0; --t) convolve_step(&suf[(t + 1) * S], &suf[t * S], S, q2c[redge[s0 + t]]);
        for (long t = 0; t < deg; ++t) {
          const double* A = &pre[t * S];
          const double* B = &suf[(t + 1) * S];
          double m0 = 0, m1 = 0;  // P(others sum to rj), P(others sum to rj - 1)
          for (int s = 0; s <= rj; ++s) m0 += A[s] * B[rj - s];
          for (int s = 0; s < rj; ++s) m1 += A[s] * B[rj - 1 - s];
          double llr;
          if (m1 <= 0 && m0 <= 0) llr = 0;
          else if (m1 <= 0) llr = -clamp_llr;
          else if (m0 <= 0) llr = clamp_llr;
          else llr = std::clamp(std::log(m1) - std::log(m0), -clamp_llr, clamp_llr);
          c2v[redge[s0 + t]] = static_cast<float>(llr);
        }
      }
      for (long v = 0; v < nv; ++v) {
        double s = prior;
        for (int t = 0; t < k_; ++t) s += c2v[static_cast<size_t>(v) * k_ + t];
        total_llr_[v] = s;
        for (int t = 0; t < k_; ++t) {
          double p = 1.0 / (1.0 + std::exp(-(s - c2v[static_cast<size_t>(v) * k_ + t])));
          float& q = q2c[static_cast<size_t>(v) * k_ + t];
          q = static_cast<float>(config_.damping * q + (1 - config_.damping) * p);
        }
      }
      bool satisfied = true;
      for (size_t a = 0; a < rows.size() && satisfied; ++a) {
        int32_t count = 0;
        for (int64_t p = rstart[a]; p < rstart[a + 1]; ++p) count += total_llr_[redge[p] / k_] > 0;
        satisfied = count == r_[rows[a]];
      }
      if (satisfied) {
        for (long v = 0; v < nv; ++v) set_value(vars_[v], total_llr_[v] > 0 ? 1 : 0);
        return true;
      }
    }
    return false;
  }

  void decimate() {
    const long nv = static_cast<long>(vars_.size());
    if (nv == 0) return;
    std::vector<std::pair<double, int32_t>> confidence;
    confidence.reserve(nv);
    for (long v = 0; v < nv; ++v) confidence.push_back({-std::fabs(total_llr_[v]), static_cast<int32_t>(v)});
    long nfix = std::max(1L, static_cast<long>(nv * config_.decimate_fraction));
    std::partial_sort(confidence.begin(), confidence.begin() + nfix, confidence.end());
    for (long i = 0; i < nfix; ++i) {
      int32_t e = vars_[confidence[i].second];
      if (state_[e] < 0) set_value(e, total_llr_[confidence[i].second] > 0 ? 1 : 0);
    }
  }

  static void convolve_step(const double* cur, double* nxt, int S, double q) {
    double mx = 0;
    for (int s = 0; s < S; ++s) {
      nxt[s] = cur[s] * (1 - q) + (s ? cur[s - 1] * q : 0.0);
      mx = std::max(mx, nxt[s]);
    }
    if (mx > 0)
      for (int s = 0; s < S; ++s) nxt[s] /= mx;  // scale cancels in the per-edge ratio
  }

  int32_t col(long e, int t) const { return cols_[static_cast<size_t>(e) * k_ + t]; }

  DoroCodeT* code_;
  ScBpConfig config_;
  int k_;
  long l_, n_ = 0;
  std::vector<IndexType> elements_;
  std::vector<int32_t> cols_;
  std::vector<int64_t> row_start_;
  std::vector<int32_t> row_elems_;
  // peeling state: -1 undecided, 0, 1
  std::vector<int8_t> state_;
  std::vector<int32_t> r_, u_;
  std::vector<uint8_t> in_queue_;
  std::vector<int32_t> queue_;
  long n_unknown_ = 0, n_ones_ = 0, expected_ones_ = 0;
  // BP scratch
  std::vector<int32_t> vars_;
  std::vector<double> total_llr_;
  long bp_iterations_ = 0;
  std::unordered_map<IndexType, ArrType> result_;
};

}  // namespace Doro
