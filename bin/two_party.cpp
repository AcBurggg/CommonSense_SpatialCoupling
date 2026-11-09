#include <doro/bch_wrapper.hpp>
#include <doro/doro_decoder.hpp>
#include <doro/doro.hpp>
#include <doro/load_file.hpp>
#include <doro/probability.hpp>
#include <doro/rans_wrapper.hpp>

#include "oniakDataStructure/ohist.h"
#include "oniakDebug/odebug.h"
#include "oniakMath/overylarge.h"
#include "oniakRandom/orand.h"
#include "oniakTimer/otime.h"
#include "nlohmann/json.hpp"
#include <doro/iblt.h>
#include <doro/sc_decoder.hpp>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <format>
#include <iostream>
#include <print>
#include <random>
#include <ranges>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <boost/math/distributions/binomial.hpp>

using namespace std;
using namespace Doro;
using namespace ONIAK;
using namespace nlohmann;
using namespace boost::math;

enum class Party { Alis = 0, Bela = 1 };
enum class Status { CollisionAvoiding = 0, CollisionResolving = 1, Finished = 2 };
using CounterType = int16_t;
using IndexType = VeryLargeInt<256>;
using DoroCodeType = DoroCode<IndexType, CounterType>;

// The following constants are used for parameter tuning
constexpr double diff_coding_error_min = 1e-9;
constexpr double diff_coding_error_max = 0.1;
constexpr double diff_coding_error_final = 1e-6;
constexpr double bch_block_error_rate = 0.01;
constexpr int iblt_element_margin = 100;
constexpr double iblt_element_factor = 1.25;
constexpr int iblt_cell_size = sizeof(IndexType) * 8 + 32; // 32 key-check bits.

// get the current decoding status and error counts.
// every key in result is regarded as valid as long as value is not 0.
bool get_sizes(
  const std::unordered_map<IndexType, CounterType>& result1,
  const std::unordered_map<IndexType, CounterType>& result2,
  const std::unordered_set<IndexType>& setA_minus_B, 
  const std::unordered_set<IndexType>& setB_minus_A,
  json& config, const DoroCodeType& doro
) {
  size_t A_minus_B_remaining_size = setA_minus_B.size(), B_minus_A_remaining_size = setB_minus_A.size(),
    A_intersect_B_remaining_size = 0;
  bool debug_flag = config.contains("debug") && config.at("debug");
  for (auto [key, value] : result1) {
    if (value == 0) continue;
    if (setA_minus_B.contains(key)) --A_minus_B_remaining_size;
    else if (setB_minus_A.contains(key)) --B_minus_A_remaining_size;
    else {
      ++A_intersect_B_remaining_size;
      if (debug_flag) {
        std::print(std::cout, "{} \t A intersect B in this\t", key);
        doro.print_key(key);
      }
    }
  }
  for (auto [key, value] : result2) {
    if (value == 0) continue;
    if (setA_minus_B.contains(key)) --A_minus_B_remaining_size;
    else if (setB_minus_A.contains(key)) --B_minus_A_remaining_size;
    else {
      ++A_intersect_B_remaining_size;
      if (debug_flag) {
        std::print(std::cout, "{} \t A intersect B in other\t", key);
        doro.print_key(key);
      }
    }
  }
  config["A minus B remaining"].push_back(A_minus_B_remaining_size);
  config["B minus A remaining"].push_back(B_minus_A_remaining_size);
  config["A intersect B remaining"].push_back(A_intersect_B_remaining_size);
  std::cout << "A minus B remaining: " << A_minus_B_remaining_size << std::endl;
  std::cout << "B minus A remaining: " << B_minus_A_remaining_size << std::endl;
  std::cout << "A intersect B remaining: " << A_intersect_B_remaining_size << std::endl;

  if (debug_flag) {
    for (auto value : setA_minus_B) {
      if ((!result1.contains(value) || result1.at(value) == 0)
        && (!result2.contains(value) || result2.at(value) == 0)) {
          std::print(std::cout, "{} \t A minus B\t", value);
        doro.print_key(value);
      }
    }
    for (auto value : setB_minus_A) {
      if ((!result1.contains(value) || result1.at(value) == 0)
        && (!result2.contains(value) || result2.at(value) == 0)) {
        std::print(std::cout, "{} \t B minus A\t", value);
        doro.print_key(value);
      }
    }
  }
  return A_minus_B_remaining_size == 0 && B_minus_A_remaining_size == 0 && A_intersect_B_remaining_size == 0;
}

// automatically decide the signature lengths
std::pair<int, int> signature_length(double A_minus_B_size, double B_minus_A_size, double failure_rate) {
  // signature is not needed if difference is one-sided.
  if (B_minus_A_size < 1 || A_minus_B_size < 1) return { 1, 0 };
  B_minus_A_size *= 1.25; // account for erroneous elements in the intersection
  double alpha = A_minus_B_size / B_minus_A_size;
  auto lam = [alpha, failure_rate, B_minus_A_size](double beta) {
    return -loss_func(alpha, 1.0 / beta, B_minus_A_size, failure_rate);
    };
  double p = convex_argmax(lam, /*lb*/ 2.0, /*ub*/ 50.0, /*epsilon*/ 0.01);
  int finger_s = std::ceil(B_minus_A_size * p);
  int finger_l = finger_l_size(alpha, 1.0 / p, failure_rate, finger_s);
  return { finger_s, finger_l };
}

Skellam moment_fit_skellam(const DoroCodeType& code) {
  auto [mean, variance] = ONIAK::sample_mean_variance(code.code());
  double mu1, mu2;
  if (code.is_cbf()) {
    mu1 = (mean + variance) / 2.0;
    mu2 = (variance - mean) / 2.0;
    if (mu1 < 0) mu1 = 0;
    if (mu2 < 0) mu2 = 0;
  }
  else {
    if (variance < 0) variance = 0;
    mu1 = mu2 = variance / 2.0;
  }
  return { mu1, mu2 };
}

struct bch_parameter {
  int order, capacity, midpoint, interval;
  double rate;
};

struct doro_parameter {
  int lb, ub;
  vector<bch_parameter> bch;
};

void print_doro_parameter(const doro_parameter& para) {
  cout << "lb: " << para.lb << ", ub: " << para.ub << endl;
  for (auto [i, bch] : views::enumerate(para.bch)) {
    println("BCH layer {}: order: {}, capacity: {}, midpoint: {}, rate: {}, interval: {}",
      i, bch.order, bch.capacity, bch.midpoint, bch.rate, bch.interval);
  }
}

double cost_estimation(doro_parameter para) {
  double uniform_cost = log2(para.ub - para.lb);
  for (auto& bch : para.bch) {
    uniform_cost += bch.rate;
  }
  return uniform_cost;
}

// determine the lb, ub, bch order and bch capacity automatically
doro_parameter doro_parameter_tuning(int d, int k, int A_minus_B_size, int B_minus_A_size,
  double bch_cap = bch_block_error_rate,
  double diff_coding_error = diff_coding_error_final) {
  double mu1 = static_cast<double>(A_minus_B_size) * k / static_cast<double>(d);
  double mu2 = static_cast<double>(B_minus_A_size) * k / static_cast<double>(d);
  Skellam skellam = { mu2, mu1 };
  auto [cdf, k2, k1] = skellam.cdf_map();
  int mean = std::floor(mu2 - mu1);

  doro_parameter best_para;
  double best_cost = std::numeric_limits<float>::infinity();
  assert(cdf.at(k1) >= 1.0 - 1e-6 && (cdf.at(k1) <= 1.0 + 1e-6));

  for (int ub = mean + 1; ub <= k1; ++ub) {
    for (int lb = mean; lb > k2; --lb) {
      double diff_err = 1.0 - cdf.at(ub) + cdf.at(lb);
      if (diff_err < diff_coding_error_min || diff_err > diff_coding_error_max) continue;
      doro_parameter cur_para = { lb, ub };
      int interval = ub - lb;
      int cur_lb = lb, cur_ub = ub;

      while (diff_err > diff_coding_error) {
        int midpoint = INT_MAX; // default midpoint
        double prob_in_range = 0.0;
        // We assume the original difference is between [mid - interval, mid + interval).
        for (int mid = cur_lb; mid <= cur_ub; ++mid) {
          auto val_retriever = [mean](auto map, int key) {
            if (map.contains(key)) return map.at(key);
            else return (key > mean) ? 1.0 : 0.0;
            };
          double cur_prob_in_range = val_retriever(cdf, mid + interval) - val_retriever(cdf, mid - interval);
          if (cur_prob_in_range > prob_in_range) {
            prob_in_range = cur_prob_in_range;
            midpoint = mid;
          }
        }
        assert(prob_in_range > 1.0 - diff_err);

        double cur_layer_rate = std::numeric_limits<float>::infinity();
        bch_parameter cur_layer;
        for (int bch_order = 5; bch_order <= 15; ++bch_order) {
          int bch_code_length = (1 << bch_order) - 1;
          int bch_capacity = ceil(quantile(binomial(bch_code_length, diff_err), 1.0 - bch_block_error_rate));
          if (bch_capacity < 2) bch_capacity = 2;
          double bch_check_length = bch_order * bch_capacity;
          // Cannot find a valid bch code to correct this number of errors
          if (bch_code_length <= bch_check_length) continue;
          double bch_rate = bch_check_length / (bch_code_length - bch_check_length);

          if (bch_rate < cur_layer_rate) {
            cur_layer_rate = bch_rate;
            cur_layer = { bch_order, bch_capacity, midpoint, interval, bch_rate };
          }
        }
        if (!isfinite(cur_layer_rate)) break;
        cur_para.bch.push_back(cur_layer);

        interval *= 2;
        diff_err = 1.0 - prob_in_range;
        cur_lb = midpoint - interval;
        cur_ub = midpoint + interval;
      }
      if (diff_err > diff_coding_error) continue;
      double estimated_cost = cost_estimation(cur_para);
      if (estimated_cost < best_cost) {
        best_cost = estimated_cost;
        best_para = cur_para;
      }
    }
  }

  if (!isfinite(best_cost)) {
    cerr << "Error: Could not determine the best parameters." << endl;
    exit(1);
  }
  return best_para;
}

// since the value could be negative, taking an integer division does not work.
// first cast to int, then to unsigned, as suggested by standard.
uint8_t get_bch_data(int val, float interval) {
  int bch_data = floor(val / interval);
  return bch_data;
}

// it is decoder who has unresolved elements
// return value is (resolving_cost, extra_comm_rounds)
std::pair<int, int> resolve_collisions(DoroDecoder<DoroCodeType>& decoder, DoroDecoder<DoroCodeType>& other_decoder,
  int finger_s, int finger_l) {
  int unresolved_size = decoder.unresolved_elements().size();
  double resolving_cost = unresolved_size * (log2_ceil(finger_s) + finger_l);
  if (unresolved_size > 0) {
    int extra_comm_rounds = 1; // need an additional resolving round
    int number_collisions = other_decoder.resolve_collision(decoder);
    if (number_collisions > 0) {
      resolving_cost += unresolved_size;
      extra_comm_rounds = 2;   // resolving round needs feedback
    }
    return { resolving_cost, extra_comm_rounds };
  }
  else return { 0, 0 };
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    println("Error: please specify config json file.");
    exit(1);
  }
  std::ifstream fin(argv[1]);
  json config = json::parse(fin);

  size_t seed = 0;
  if (config.contains("seed")) seed = config.at("seed");

  auto [setA, setB, setA_minus_B, setB_minus_A, A_minus_B_size, B_minus_A_size, A_intersect_B_size] 
  = load_dataset_k256_or_k32<IndexType>(config);
  auto B_size = setB.size();

  int k = config.at("k");
  int d = config.at("d");
  // Spatial coupling and decoder choice (see sc_study/REPORT.md). Defaults keep the original
  // protocol: uniform matrix and the MP ping-pong decoder.
  int sc_window = 0;
  if (config.contains("sc window")) sc_window = config.at("sc window");
  else if (config.contains("sc window fraction"))
    sc_window = static_cast<int>(std::lround(config.at("sc window fraction").get<double>() * d));
  bool sc_tle = config.contains("sc tle") && config.at("sc tle").get<bool>();
  std::string decoder_choice = config.contains("decoder") ? config.at("decoder").get<std::string>() : "mp";
  // With coupling, the bulk counters carry the load of d - sc_window + 1 positions.
  int bulk_counters = (sc_window > 0) ? d - sc_window + 1 : d;
  std::string save_path = config.at("result filename");
  int max_comm_rounds = 20;
  if (config.contains("max comm rounds")) max_comm_rounds = config.at("max comm rounds");
  // at this round, always start resolving
  // by default, resolving starts after two complete back-and-forths
  int resolving_round = 3;
  if (config.contains("resolving round")) resolving_round = config.at("resolving round");
  bool counting = config.at("counting");
  double bch_cap = bch_block_error_rate;
  if (config.contains("bch block error rate")) bch_cap = config.at("bch block error rate");
  double diff_coding_error = diff_coding_error_final;
  if (config.contains("diff coding error")) diff_coding_error = config.at("diff coding error");

  // Skip this experiment if result already exists, used for batch experimenting.
  bool skip_if_exists = false;
  double failure_rate = 1e-7;
  if (config.contains("failure rate")) failure_rate = config.at("failure rate");
  failure_rate /= B_minus_A_size; // from now on, we 
  int max_recenter_rounds = 100;
  if (config.contains("max recenter rounds")) max_recenter_rounds = config.at("max recenter rounds");
  bool force_pursue_l1 = false;
  if (config.contains("force pursue l1")) force_pursue_l1 = config.at("force pursue l1");
  bool track_values = true;
  int max_num_peels = -1;
  if (config.contains("max num peels")) max_num_peels = config.at("max num peels");
  if (config.contains("track values"))
    track_values = config.at("track values");
  if (max(setA.size(), setB.size()) > 10000000)
    track_values = false;

  if (config.contains("skip if exists")) skip_if_exists = config.at("skip if exists");
  if (skip_if_exists && filesystem::exists(save_path)) {
    cout << "Notice: Experiment results already exist. Skipping..." << endl;
    return 0;
  }
  std::ofstream fout(save_path);
  if (!fout.is_open()) {
    cout << "Error: cannot open file " << save_path << endl;
    return 1;
  }

  auto rng = std::mt19937(seed);

  auto auto_parameter = doro_parameter_tuning(bulk_counters, k, A_minus_B_size, B_minus_A_size, bch_cap, diff_coding_error);
  cout << "Automatically selected the following parameters: " << endl;
  print_doro_parameter(auto_parameter);

  DoroCodeType doro(d, k, counting, rng, auto_parameter.lb, auto_parameter.ub, track_values);
  if (sc_window > 0) doro.enable_spatial_coupling(sc_window, sc_tle, rng());
  // copy to keep the same set of hash functions
  auto first_round_code = doro.copy();
  unordered_map<IndexType, CounterType> Bela_coef, Alis_coef;
  for (auto& i : setA) {
    Alis_coef[i] = 1;
  }
  for (auto& i : setB) {
    Bela_coef[i] = 1;  // in Bela
  }
  first_round_code.encode(Alis_coef);

  vector<BCHWrapper> bch_wrappers;
  for (const auto& bch : auto_parameter.bch) {
    bch_wrappers.emplace_back(bch.order, bch.capacity);
  }
  vector<uint8_t> bch_data;
  vector<vector<uint8_t>> parity_bits;
  float interval_size = first_round_code.interval();
  for (auto& val : first_round_code.code()) {
    bch_data.push_back(get_bch_data(val, interval_size));
    val = first_round_code.recenter(val);
  }
  for (auto [i, bch] : views::enumerate(bch_wrappers)) {
    // Only LSB is encoded in each iteration. So we reveal the LSB after that.
    parity_bits.push_back(bch.encode(bch_data, /*bit-by-bit*/ true));
    for_each(bch_data.begin(), bch_data.end(), [](uint8_t& x) { x >>= 1; });
  }
  bch_data = {};

  size_t total_bch_sizes = accumulate(parity_bits.begin(), parity_bits.end(), 0,
    [](size_t sum, const auto& vec) { return sum + vec.size(); });
  // rANS encode for uniform distribution in [lb, ub)
  auto uniform_map = uniform_pmf<CounterType>(auto_parameter.lb, auto_parameter.ub);
  RansWrapper first_round_wrapper(uniform_map);
  RansCode first_round_compressed_code = first_round_wrapper.encode(first_round_code.code());
  auto first_round_decompressed_code = first_round_wrapper.decode(first_round_compressed_code);
  DEBUG_VECTOR_EQUAL(first_round_decompressed_code, first_round_code.code());
  double first_round_cost = first_round_compressed_code.size() + total_bch_sizes;
  first_round_code.code() = std::move(first_round_decompressed_code);
  first_round_compressed_code.clear();

  doro.encode(Bela_coef);
  auto doro_value = doro.values();
  if (doro_value) {
    for (auto &item : setA) {
      if (doro_value->contains(item))
        doro_value->operator[](item) -= 1;
      else
        doro_value->operator[](item) = -1;
    }
  }
  vector<uint8_t> data_decode;
  vector<CounterType> diff_vec;
  for (auto [i, val] : views::enumerate(doro.code())) {
    CounterType diff = val - first_round_code.code()[i];
    diff = doro.recenter(diff);
    diff_vec.push_back(diff);
    CounterType alis_val = val - diff; // Guess of original code by Alis.
    data_decode.push_back(get_bch_data(alis_val, interval_size));
  }
  first_round_code.reset();
  int num_failed_blocks = 0, num_errors_found = 0;
  for (auto [i, bch] : views::enumerate(bch_wrappers)) {
    auto bch_param = auto_parameter.bch[i];
    auto data_decode_copy = data_decode;
    for_each(data_decode_copy.begin(), data_decode_copy.end(), [i](uint8_t& x) { x >>= i; });
    vector<size_t> error_pos =
      bch.get_error_positions(data_decode_copy, parity_bits[i], /*bit-by-bit*/ true, &num_failed_blocks, &num_errors_found);
    for (size_t pos : error_pos) {
      CounterType& diff = diff_vec[pos];
      diff += (diff < bch_param.midpoint) ? bch_param.interval : -bch_param.interval;
      CounterType alis_val = doro.code()[pos] - diff;
      data_decode[pos] = get_bch_data(alis_val, interval_size);
    }
  }
  config["number of failed blocks"] = num_failed_blocks;
  config["number of corrected errors"] = num_errors_found;
  assert(diff_vec.size() == doro.code().size());
  // We only look at the difference from this moment on.
  doro.code() = std::move(diff_vec);
  parity_bits = {};
  data_decode = {};

  double lambda = static_cast<double>(A_minus_B_size) * k / d;
  auto A_minus_B_map = get_pmf(lambda, d, counting);

  lambda = static_cast<double>(A_intersect_B_size) * k / d;
  auto A_intersect_B_map = get_pmf(lambda, d, counting);
  lambda = static_cast<double>(B_minus_A_size) * k / d;
  auto B_minus_A_map = get_pmf(lambda, d, counting);
  lambda = static_cast<double>(B_size) * k / d;
  auto B_map = get_pmf(lambda, d, counting);
  double B_entropy = entropy(B_map);
  double all_entropy = doro_entropy(A_intersect_B_map, A_minus_B_map, B_minus_A_map);
  // This is conditional entropy H(A | B).
  double first_round_entropy_cost = (all_entropy - B_entropy) * d;

  config["comm costs"] = json::array({ first_round_cost });
  config["doro costs"] = json::array({ first_round_cost });
  config["theoretical entropy costs"] = json::array({ first_round_entropy_cost });
  config["finger costs"] = json::array({ 0 });
  config["resolving costs"] = json::array({ 0 });
  config["time"] = json::array({});
  config["num peels"] = json::array({});
  config["A minus B remaining"] = json::array({});
  config["B minus A remaining"] = json::array({});
  config["A intersect B remaining"] = json::array({});
  config["number of recenters"] = 0;

  Party party = Party::Alis;
  auto [finger_s, finger_l] = signature_length(A_minus_B_size, B_minus_A_size, failure_rate);
  double finger_beta = static_cast<double>(B_minus_A_size) / finger_s;
  println("Use the following signature lengths: s = {}, l = {}, beta = {}.", finger_s, finger_l, finger_beta);
  assert(finger_s > 0 && finger_l >= 0 and finger_l < 64);
  WYHash finger_hash(/*mask*/ 0, /*mode*/ finger_s, /*seed*/ rng());
  unsigned int mask2 = (1 << (finger_l % 32)) - 1;
  WYHash resolving_hash(mask2, /*mode*/ 0, /*seed*/ rng(), /*seed2*/ rng(), /*mode64*/ finger_l >= 32);

  Status status = Status::CollisionAvoiding;
  DecodeConfig dconf_alis(/*verbose*/ false, /*debug*/ false, /*lb*/ -1,
                          /*ub*/ 0, max_num_peels, PursuitChoice::L2),
      dconf_bela(/*verbose*/ false, /*debug*/ false, /*lb*/ 0, /*ub*/ 1,
                 max_num_peels, PursuitChoice::L2);
  if (A_minus_B_size == 0) dconf_alis.lb = 0;
  if (B_minus_A_size == 0) dconf_bela.ub = 0;
  if (force_pursue_l1) {
    dconf_alis.pursuit_choice = PursuitChoice::L1;
    dconf_bela.pursuit_choice = PursuitChoice::L1;
  }
  DoroDecoder<DoroCodeType> decoder_alis(doro, setA, dconf_alis, max_recenter_rounds, &finger_hash, &resolving_hash),
    decoder_bela(doro, setB, dconf_bela, max_recenter_rounds, &finger_hash, &resolving_hash);
  StopWatch sw;
  int comm_rounds = 0, actual_comm_rounds = 0, round_in_resolve = 0;
  bool success = false;

  if (resolving_round < 0) {
    decoder_alis.enter_resolving();
    decoder_bela.enter_resolving();
  }

  // Peeling + BP decoder (unidirectional only). On success the protocol finishes in one round;
  // otherwise the code is untouched and the original MP rounds / IBLT fallback run as before.
  bool primary_decoder_success = false;
  if (decoder_choice == "bp") {
    if (A_minus_B_size != 0 || !counting) {
      println("Error: decoder \"bp\" currently supports unidirectional SetX (A subset of B) with counting codes.");
      exit(1);
    }
    ScBpDecoder<DoroCodeType> sc_decoder(doro, setB);
    bool decoded = sc_decoder.decode();
    primary_decoder_success = decoded;
    config["bp iterations"] = sc_decoder.bp_iterations();
    config["bp undecided"] = sc_decoder.num_undecided();
    if (decoded) {
      actual_comm_rounds = 1;
      config["time"].push_back(sw.peek());
      config["num peels"].push_back(doro.num_peels());
      success = get_sizes(sc_decoder.result(), decoder_alis.result(), setA_minus_B, setB_minus_A, config, doro);
      status = Status::Finished;
    } else {
      cout << "BP decoder did not reach zero residual. Falling back to MP." << endl;
    }
  }

  while (!doro.empty() && comm_rounds < max_comm_rounds && status != Status::Finished) {
    actual_comm_rounds = std::max(actual_comm_rounds, comm_rounds + 1);
    // who is current decoder?
    party = (party == Party::Alis) ? Party::Bela : Party::Alis;
    DoroDecoder<DoroCodeType>& decoder = (party == Party::Alis) ? decoder_alis : decoder_bela,
      & other_decoder = (party == Party::Alis) ? decoder_bela : decoder_alis;
    int num_peels = decoder.decode();

    if (doro.empty() && decoder.unresolved_elements().empty()) {
      config["time"].push_back(sw.peek());
      config["num peels"].push_back(num_peels);
      success = success || get_sizes(decoder.result(), other_decoder.result(), setA_minus_B, setB_minus_A, config, doro);
      primary_decoder_success = success;
      break;
    }

    Skellam code_distribution = moment_fit_skellam(doro);
    auto code_frequency = code_distribution.pmf_map<CounterType>();
    double doro_cost = 0;
    RansWrapper rans_wrapper(code_frequency);
    RansCode doro_compressed_code = rans_wrapper.encode(doro.code());
    auto decompressed_code = rans_wrapper.decode(doro_compressed_code);

    double skellam_entropy = entropy(frequency_count(doro.code()));
    double theoretical_entropy = skellam_entropy * d;
    doro_cost = doro_compressed_code.size() + 64;
    doro_compressed_code.clear();
    // 8 bytes to transmit two Skellam parameters in float.
    config["theoretical entropy costs"].push_back(theoretical_entropy);
    config["doro costs"].push_back(doro_cost);
    DEBUG_VECTOR_EQUAL(decompressed_code, doro.code());
    doro.code() = std::move(decompressed_code);  // just use these codes.

    auto fingerprints_delta = decoder.update_fingerprints();
    double new_extra_size = std::accumulate(fingerprints_delta.begin(), fingerprints_delta.end(), 0.0);
    float actual_load = new_extra_size / finger_s;
    unordered_map<int8_t, double> finger_pmf = { {0, 1.0 - actual_load}, {1, actual_load} };
    RansWrapper<int8_t, double> finger_rans(finger_pmf);
    double finger_cost = 0;
    if (new_extra_size >= 1 && finger_s > 1) {
      RansCode compressed_fingers = finger_rans.encode(fingerprints_delta);
      finger_cost = compressed_fingers.size() + 32; // 4 bytes to transmit actual load.
      auto decompressed_fingers = finger_rans.decode(compressed_fingers);
      other_decoder.load_fingerprints(decompressed_fingers);
      DEBUG_VECTOR_EQUAL(decompressed_fingers, fingerprints_delta);
    }

    auto [resolving_cost, extra_comm_rounds] = resolve_collisions(decoder, other_decoder, finger_s, finger_l);
    actual_comm_rounds = std::max(actual_comm_rounds, comm_rounds + extra_comm_rounds + 1);
    config["finger costs"].push_back(finger_cost);
    config["resolving costs"].push_back(resolving_cost);
    double comm_cost = doro_cost + finger_cost + resolving_cost;
    config["comm costs"].push_back(comm_cost);
    success = success || get_sizes(decoder.result(), other_decoder.result(),
      setA_minus_B, setB_minus_A, config, doro);
    config["time"].push_back(sw.peek());
    config["num peels"].push_back(doro.num_peels());

    if (!decoder.has_new_elements() || comm_rounds == resolving_round) {
      if (status == Status::CollisionAvoiding) {
        status = Status::CollisionResolving;
        decoder_alis.enter_resolving();
        decoder_bela.enter_resolving();
        config["round entering resolving"] = actual_comm_rounds;
      } 
      if (!decoder.has_new_elements() && round_in_resolve >= 2) {
        status = Status::Finished;
      }
      ++round_in_resolve;
    }
    if (status == Status::Finished || comm_rounds == max_comm_rounds - 1) {
      // final IBLT stage.
      if (!doro.empty()) {
        cout << "Doro failed. Fall back to IBLT." << endl;
        int estimated_diff = sample_mean_variance(doro.code()).second * doro.code().size() * iblt_element_factor + iblt_element_margin;
        config["estimated diff"] = estimated_diff;
        IBLT<IndexType> iblt(estimated_diff, sizeof(IndexType));
        decoder.encode_iblt(iblt);
        comm_rounds = actual_comm_rounds;
        ++actual_comm_rounds;
        int iblt_size = iblt_cell_size * iblt.hashTableSize();
        config["comm costs"].push_back(iblt_size);

        // IBLT reaches the other party. My extra elements are returned.
        vector<IndexType> my_elements = other_decoder.decode_iblt(std::move(iblt));
        int my_element_size = my_elements.size() * sizeof(IndexType) * 8;
        if (!my_elements.empty()) ++actual_comm_rounds;
        config["iblt size"] = iblt_size + my_element_size;
        decoder.add_iblt_elements(my_elements);

        // finally resolve collisions between decoder and other
        auto [resolving_cost, extra_comm_rounds] = resolve_collisions(other_decoder, decoder, finger_s, finger_l);
        actual_comm_rounds = std::max(actual_comm_rounds, comm_rounds + extra_comm_rounds + 1);

        // update reconciliation success and costs
        config["resolving costs"].push_back(resolving_cost);
        config["comm costs"].push_back(my_element_size + resolving_cost);
        success = success || get_sizes(decoder.result(), other_decoder.result(),
          setA_minus_B, setB_minus_A, config, doro);
        config["time"].push_back(sw.peek());
      }
      status = Status::Finished;
    }
    ++comm_rounds;
  }

  config["number of recenters"] = doro.num_recenters();
  config["communication rounds"] = actual_comm_rounds;
  double total_comm_cost = 0;
  for (auto cost : config.at("comm costs")) {
    total_comm_cost += cost.get<double>();
  }
  config["total communication cost"] = total_comm_cost;
  config["total theoretical entropy cost"] = std::accumulate(config.at("theoretical entropy costs").begin(),
    config.at("theoretical entropy costs").end(), 0.0, [](double sum, const json& val) { return sum + val.get<double>(); });
  config["primary decoder success"] = primary_decoder_success;
  config["success"] = success;

  FILE *file = fopen("/proc/self/status", "r");  // Linux only
  if (file != NULL) {
    char line[128];
    while (fgets(line, 128, file) != NULL) {
      if (strstr(line, "VmHWM:") != NULL) {
        config["peak memory"] = line;
      }
    }
    fclose(file);
  }
  fout << std::setw(4) << config << endl;  // indent 4
  return 0;
}