#ifndef __DORO_HPP__
#define __DORO_HPP__

#include "oniakHash/ohash.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <inttypes.h>
#include <iostream>
#include <numeric>
#include <random>
#include <ranges>
#include <unordered_map>
#include <memory>
#include <vector>

namespace Doro {

constexpr int kDoroNumHashFuncs = 30;
constexpr int kDoroNumHashFuncsMargin = 10;

template <typename IndexType, typename ArrType = int32_t>
class DoroCode {
public:
  using SparseVector = std::unordered_map<IndexType, ArrType>;
  using PairType = std::pair<IndexType, ArrType>;
  using IndexT = IndexType;
  using ArrayT = ArrType;

  template <typename RandomDevice>
  DoroCode(int size, int k, bool is_cbf, RandomDevice& rng, int lb = 0, int ub = 0, bool track_value = true) :
    arr_(size), is_cbf_(is_cbf), k_(k), num_peels_(0), num_correct_peels_(0),
    num_recenters_(0), lb_(lb), ub_(ub), interval_(ub - lb) {
    int num_hash_funcs = std::max(kDoroNumHashFuncs, k + kDoroNumHashFuncsMargin);
    hash_funcs_.reserve(num_hash_funcs);
    std::ranges::for_each(std::views::iota(0, num_hash_funcs),
      [&](int) {hash_funcs_.emplace_back(/*mask*/ 0, /*mod*/ 2 * size, /*seed*/ rng());});
    if (track_value) {
      values_ = std::make_unique<SparseVector>();
    }
  }

  void encode(const SparseVector& kvpairs) {
    assert(!values_ || values_->empty());
    for (auto [x, val] : kvpairs) {
      auto all_hashes = hash_all(x);
      for (auto [index, sign] : all_hashes) {
        arr_[index] += sign * val;
      }
      if (values_) values_->operator[](x) = val;
    }
  }

  // if delta is 1, all counters will be REDUCED by 1 (if positive sign)
  // or INCREASED by 1 (if negative sign)
  void peel(IndexType element, ArrType delta) {
    num_peels_ += 1;
    auto all_hashes = hash_all(element);
    for (auto [index, sign] : all_hashes) {
      arr_[index] -= sign * delta;
    }
    if (values_) {
    auto element_iter = values_->find(element);
    ArrType ori_value = (element_iter != values_->end()) ? element_iter->second : 0;
    ArrType new_value = ori_value - delta;
    values_->operator[](element)= new_value;
    if (std::abs(new_value) < std::abs(ori_value))
      num_correct_peels_ += 1;
    }
  }

  // senses the signal of an element by inner product
  // was performance bottleneck
  ArrType sense(IndexType element) const {
    ArrType signal = 0;
    auto all_hashes = hash_all(element);
    for (auto [index, sign] : all_hashes) {
      signal += sign * arr_[index];
    }
    return signal;
  }

  void print_key(IndexType element) const {
    auto all_hashes = hash_all(element);
    for (auto [index, sign] : all_hashes) {
      std::print(std::cout, "{}\\{} ", index, sign * arr_[index]);
    }
    std::cout << std::endl;
  }

  // recenter a counter so that it lies in [lb, ub)
  ArrType recenter(ArrType x) {
    if (interval_ <= 0) return x;   // never recenters
    x -= std::floor(static_cast<float>(x - lb_) / interval_) * interval_;
    assert(x >= lb_ && x < ub_);
    return x;
  }

  ArrType deviation(ArrType x) const {
    if (interval_ <= 0 || (x >= lb_ && x < ub_)) return 0;
    return std::max(lb_ - x, x - ub_ + 1);
  }

  // Returns median signal. If k_ is even, returns the smaller one in absolute value.
  ArrType sense_l1(IndexType element) const {
    std::vector<ArrType> signals;
    auto all_hashes = hash_all(element);
    for (auto [index, sign] : all_hashes) {
      signals.push_back(sign * arr_[index]);
    }
    ArrType median = get_median(element);
    double old_l1 = std::accumulate(signals.begin(), signals.end(), 0.0, [](double x, ArrType& y){return x + std::abs(y);});
    double new_l1 = std::accumulate(signals.begin(), signals.end(), 0.0, [median](double x, ArrType& y){return x + std::abs(y - median);});
    return old_l1 - new_l1;
  }

  ArrType get_median(IndexType element) const {
    std::vector<ArrType> signals;
    auto all_hashes = hash_all(element);
    for (auto [index, sign] : all_hashes) {
      signals.push_back(sign * arr_[index]);
    }
    ArrType median = 0.0;
    auto med_iter = signals.begin() + k_ / 2;
    std::nth_element(signals.begin(), med_iter, signals.end());
    if (k_ % 2 == 0) {
      auto med_iter2 = med_iter - 1;
      std::nth_element(signals.begin(), med_iter2, med_iter);
      if (*med_iter * (*med_iter2) < 0)
      median = 0;
      else if (std::abs(*med_iter) > std::abs(*med_iter2))
      median = *med_iter2;
    } else median = *med_iter;
    return median;
  }

  int nonzero_num() const {
    if (!values_) return -1;
    return std::count_if(values_->begin(), values_->end(),
                         [](PairType pair) { return pair.second != 0; });
  }

  ArrType mae() const {
    if (!values_) return -1;
    return std::accumulate(
        values_->begin(), values_->end(), 0,
        [](ArrType sum, PairType pair) { return sum + std::abs(pair.second); });
  }

  void reset() {
    arr_.clear();
    values_.reset();
  }

  int size() const {
    return arr_.size();
  }

  std::vector<ArrType>& code() { return arr_; }
  const std::vector<ArrType>& code() const { return arr_; }

  int k() const { return k_; }
  int num_peels() const { return num_peels_; }
  int& num_peels() { return num_peels_; }
  int& num_correct_peels() { return num_correct_peels_; }
  int& num_recenters() { return num_recenters_; }
  ArrType interval() const { return interval_; }
  double mid_point() const { return (lb_ + ub_ - 1.0) / 2.0; }

  void show_result() const {
    std::cout << "Number of peels: " << num_peels_ << std::endl;
    std::cout << "Number of correct peels: " << num_correct_peels_ << std::endl;
    std::cout << "Number of wrong peels: " << num_peels_ - num_correct_peels_ << std::endl;
    std::cout << "Size of nonzero indices: " << nonzero_num() << std::endl;
    std::cout << "L1 Norm of indices: " << mae() << std::endl;
  }

  SparseVector* values() { return values_.get(); }

  // Spatial coupling (Walzer's scheme, as adapted by METTLE; see sc_study/REPORT.md).
  // Each element x is hashed to a position s(x) in [0, size - window], and its k counters are
  // drawn (distinct) from the window [s(x), s(x) + window). With tle, the first counter is s(x)
  // itself (METTLE's leading edge). Both parties must enable it with identical arguments before
  // encoding. The MP decoder in doro_decoder.hpp is NOT recommended with this mode; use
  // ScBpDecoder (sc_decoder.hpp) instead.
  void enable_spatial_coupling(int window, bool tle, uint32_t seed) {
    if (window < 4 * k_ || window > size())
      throw std::invalid_argument("spatial coupling window must be in [4k, size]");
    sc_window_ = window;
    sc_tle_ = tle;
    // doubled ranges carry a sign bit for non-CBF codes, like hash_funcs_
    sc_pos_hash_ = ONIAK::WYHash(/*mask*/ 0, /*mod*/ 2 * (size() - window + 1), seed);
    sc_offset_hashes_.clear();
    for (auto& hfunc : hash_funcs_)
      sc_offset_hashes_.emplace_back(/*mask*/ 0, /*mod*/ 2 * window, hfunc.seed());
  }

  int sc_window() const { return sc_window_; }
  bool sc_tle() const { return sc_tle_; }

  // each value is an (index, sign) pair
  std::vector<std::pair<int, int>> hash_all(IndexType x) const {
    if (sc_window_ > 0) return hash_all_coupled(x);
    std::vector<std::pair<int, int>> all_hashes;
    all_hashes.reserve(k_);
    for (auto& hfunc : hash_funcs_) {
      int hash_value = hfunc.hash_in_range(x);
      int idx = (hash_value >= size()) ? hash_value - size() : hash_value;
      int sign = (is_cbf_ || hash_value < size()) ? 1 : -1;
      if (!std::ranges::contains(all_hashes, idx, [](auto x) {return x.first;}))
        all_hashes.emplace_back(idx, sign);
      // otherwise reject this sample.
      if (all_hashes.size() == static_cast<size_t>(k_)) return all_hashes;
    }
    throw std::runtime_error("Error: Run out of hash functions.");
  }

  std::vector<std::pair<int, int>> hash_all_coupled(IndexType x) const {
    std::vector<std::pair<int, int>> all_hashes;
    all_hashes.reserve(k_);
    const int num_positions = size() - sc_window_ + 1;
    int pos_value = sc_pos_hash_.hash_in_range(x);
    int position = (pos_value >= num_positions) ? pos_value - num_positions : pos_value;
    if (sc_tle_)
      all_hashes.emplace_back(position, (is_cbf_ || pos_value < num_positions) ? 1 : -1);
    for (auto& hfunc : sc_offset_hashes_) {
      if (all_hashes.size() == static_cast<size_t>(k_)) return all_hashes;
      int hash_value = hfunc.hash_in_range(x);
      int offset = (hash_value >= sc_window_) ? hash_value - sc_window_ : hash_value;
      int sign = (is_cbf_ || hash_value < sc_window_) ? 1 : -1;
      int idx = position + offset;
      if (!std::ranges::contains(all_hashes, idx, [](auto x) {return x.first;}))
        all_hashes.emplace_back(idx, sign);
      // otherwise reject this sample.
    }
    if (all_hashes.size() == static_cast<size_t>(k_)) return all_hashes;
    throw std::runtime_error("Error: Run out of hash functions.");
  }

  bool is_cbf() const { return is_cbf_; }

  bool empty() const {
    return std::all_of(arr_.begin(), arr_.end(), [](ArrType val) { return val == 0; });
  }

  DoroCode copy() const {
    std::mt19937 rng;
    DoroCode copy(size(), k_, is_cbf_, rng, lb_, ub_, values_ != nullptr);
    copy.arr_ = arr_;
    copy.hash_funcs_ = hash_funcs_;
    copy.sc_window_ = sc_window_;
    copy.sc_tle_ = sc_tle_;
    copy.sc_pos_hash_ = sc_pos_hash_;
    copy.sc_offset_hashes_ = sc_offset_hashes_;
    if (values_) {
      copy.values_ = std::make_unique<SparseVector>(*values_);
    }
    return copy;
  }

private:
  std::vector<ArrType> arr_;
  std::unique_ptr<SparseVector> values_;

  // If is_cbf_ is true, then all hashed signs are positive.
  bool is_cbf_; // cbf: counting Bloom filter
  std::vector<ONIAK::WYHash> hash_funcs_;
  int k_, num_peels_, num_correct_peels_, num_recenters_;
  // [lower bounds, and upper bounds) in quantization
  // if lb >= ub, then quantization is not in effect.
  ArrType lb_, ub_;
  float interval_; // ub - lb
  // spatial coupling; disabled when sc_window_ == 0
  int sc_window_ = 0;
  bool sc_tle_ = false;
  ONIAK::WYHash sc_pos_hash_;
  std::vector<ONIAK::WYHash> sc_offset_hashes_;
};

} // namespace Doro

#endif
