#pragma once

#include "model.hpp"
#include "util.hpp"

#include <algorithm>
#include <cstdint>
#include <format>
#include <limits>
#include <logger.hpp>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace biogeosim {

typedef uint64_t dist_base_t;

class dist_t {
public:
  dist_t() = default;

  dist_t(const std::string &);

  dist_t(const dist_t &) = default;
  dist_t(dist_t &&)      = default;

  dist_t &operator=(const dist_t &) = default;
  dist_t &operator=(dist_t &&)      = default;

  dist_t &operator|=(const dist_t &d);

  dist_t(uint64_t d, uint16_t s) : _dist{d}, _regions{s} {}

  explicit dist_t(uint16_t r) : _dist{0}, _regions{r} {}

  inline constexpr size_t popcount() const {
    return static_cast<size_t>(__builtin_popcountll(_dist));
  }

  inline constexpr size_t clz() const {
    return static_cast<size_t>(__builtin_clzll(_dist));
  }

  inline constexpr size_t clz_min() const {
    return static_cast<size_t>(__builtin_clzll(_dist | 1ull));
  }

  inline constexpr size_t log2() const {
    constexpr size_t BITS_IN_BYTE = 8;
    return sizeof(_dist) * BITS_IN_BYTE - clz();
  }

  inline constexpr uint16_t regions() const { return _regions; }

  inline constexpr uint64_t operator[](size_t i) const { return bextr(i); }

  dist_t operator^(dist_t d) const {
    return {_dist ^ d._dist, std::max(_regions, d._regions)};
  }
  dist_t operator|(dist_t d) const {
    return {_dist | d._dist, std::max(_regions, d._regions)};
  }
  dist_t operator&(dist_t d) const {
    return {_dist & d._dist, std::max(_regions, d._regions)};
  }

  explicit operator uint64_t() const { return _dist; }

  dist_t operator&(uint64_t d) const { return {_dist & d, _regions}; }
  dist_t operator~() const { return {~_dist, _regions}; }

  bool operator==(dist_t d) const {
    return d._dist == _dist && _regions == d._regions;
  }

  bool operator!=(dist_t d) const { return !(d == *this); }

  dist_t negate_bit(size_t index) const {
    return {_dist ^ (1ull << index), _regions};
  }

  dist_t operator+(uint64_t d) const { return {_dist + d, _regions}; }

  explicit operator bool() const { return static_cast<bool>(_dist); }

  size_t index(size_t max_areas) const {
    size_t skips = compute_skips(_dist, max_areas);
    return _dist - skips;
  }

  dist_t next_dist(uint32_t n) const {
    auto d = *this + 1;
    while (d.popcount() > n) { d = d + 1; }
    return d;
  }

  std::string to_str() const {
    std::ostringstream oss;
    oss << *this;
    return oss.str();
  }

  friend std::ostream &operator<<(std::ostream &os, dist_t dist) {
    for (size_t i = dist._regions; i; --i) { os << dist.bextr(i - 1); }
    return os;
  }

private:
  constexpr uint64_t bextr(size_t index) const {
    return (_dist >> index) & 1ull;
  }

  static auto compute_skips_power_of_2(size_t k, size_t n) -> size_t {
    size_t skips = 0;
    for (size_t i = n + 1; i < k; ++i) {
      skips += biogeosim::util::combinations(k - 1, i);
    }
    return skips;
  }

  static auto compute_skips(size_t i, size_t n) -> size_t {
    constexpr size_t BITS_IN_BYTE = 8;
    size_t           skips        = 0;
    while (i != 0 && n != 0) {
      size_t first_index  = sizeof(i) * BITS_IN_BYTE - __builtin_clzll(i | 1);
      skips              += compute_skips_power_of_2(first_index, n);
      n                  -= 1;
      i                  -= 1 << (first_index - 1);
    }
    skips += i;
    return skips;
  }

  uint64_t _dist;
  uint16_t _regions;
};

bool valid_dist(dist_t d, const substitution_model_t &model);
bool valid_dist(dist_t d, size_t regions);

class transition_t {
public:
  transition_t() = default;
  transition_t(double t, dist_t i, dist_t f)
      : waiting_time{t}, initial_state{i}, final_state{f} {}

  double waiting_time = std::numeric_limits<double>::infinity();
  dist_t initial_state;
  dist_t final_state;
};

transition_t sample(dist_t                                  init_dist,
                    const substitution_model_t             &model,
                    std::uniform_random_bit_generator auto &gen) {
  auto [d, e] = model.rates();

  bool singleton = init_dist.popcount() == 1;

  std::exponential_distribution<double> dis_die(d);
  std::exponential_distribution<double> exp_die(e);

  std::vector<transition_t> rolls(model.region_count());

  for (size_t i = 0; i < model.region_count(); ++i) {
    if (singleton && init_dist[i]) { continue; }
    double waiting_time = init_dist[i] ? exp_die(gen) : dis_die(gen);
    rolls[i] = transition_t{waiting_time, init_dist, init_dist.negate_bit(i)};
  }

  return *std::min_element(rolls.begin(), rolls.end(), [](auto a, auto b) {
    return a.waiting_time < b.waiting_time;
  });
}

enum class split_type_e { singleton, allopatric, sympatric, jump, invalid };

struct split_t {
  dist_t       left;
  dist_t       right;
  split_type_e type;

  std::string to_nhx_string() {
    std::ostringstream oss;
    oss << "left-split=" << left << ":"
        << "right-split=" << right << ":";
    oss << "split-type=";

    return oss.str();
  }

  std::string type_string() const {
    switch (type) {
    case split_type_e::singleton:
      return "singleton";
    case split_type_e::allopatric:
      return "allopatric";
    case split_type_e::sympatric:
      return "sympatric";
    case split_type_e::jump:
      return "jump";
    case split_type_e::invalid:
      return "invalid";
    }
    throw std::runtime_error{"Did not cover all cases"};
  }
};

std::vector<transition_t>
generate_samples(dist_t                                  init_dist,
                 double                                  brlen,
                 const substitution_model_t             &model,
                 std::uniform_random_bit_generator auto &gen) {
  std::vector<transition_t> results;
  while (true) {
    auto r  = sample(init_dist, model, gen);
    brlen  -= r.waiting_time;
    if (brlen < 0.0) { return results; }
    LOG_DEBUG("adding transition from %lb to %lb",
              static_cast<uint64_t>(r.initial_state),
              static_cast<uint64_t>(r.final_state));
    init_dist = r.final_state;
    results.push_back(r);
  }
}

split_type_e determine_split_type(dist_t init_dist,
                                  dist_t left_dist,
                                  dist_t right_dist,
                                  bool   jumps_ok);

/*
 * There are three types of splitting:
 *  - Singleton
 *  - Allopatric
 *  - Sympatric
 *  Allopatric and Sympatric are not the names used in the original Ree paper,
 *  and they shouldn't be used in user facing descriptions, as they are very
 *  misleading. Regions are large enough that both Allopatry and Sympatry can
 *  occur, but the idea maps well, so I use it internally.
 *
 *  In additoin, Metzke introduced a jump parameter, making it +J. We optionally
 *  support this option.
 */
split_t split_dist(dist_t                                  init_dist,
                   const substitution_model_t             &model,
                   std::uniform_random_bit_generator auto &gen) {
  // Singleton case
  if (!model.jumps_ok() && init_dist.popcount() == 1) {
    LOG_DEBUG("Splitting a singleton: %lb", static_cast<uint64_t>(init_dist));
    return {init_dist, init_dist, split_type_e::singleton};
  }
  auto max_dist = (1ul << init_dist.regions()) - 1;
  std::uniform_int_distribution<dist_base_t> dist_gen(1, max_dist);

  std::bernoulli_distribution sympatry_coin(
      model.cladogenesis_params().sympatry);

  std::bernoulli_distribution allopatry_coin(
      model.cladogenesis_params().allopatry);

  std::bernoulli_distribution copy_coin(model.cladogenesis_params().copy);

  std::bernoulli_distribution jump_coin(model.cladogenesis_params().jump);

  dist_t       left_dist, right_dist;
  split_type_e split_type;
  size_t sample_count = 0;

  while (true) {
    sample_count++;
    left_dist  = dist_t{dist_gen(gen), init_dist.regions()};
    right_dist = dist_t{dist_gen(gen), init_dist.regions()};
    split_type = determine_split_type(
        init_dist, left_dist, right_dist, model.jumps_ok());

    if (split_type == split_type_e::invalid) { continue; }
    if (split_type == split_type_e::sympatric) {
      if (sympatry_coin(gen)) { break; }
    } else if (split_type == split_type_e::allopatric) {
      if (allopatry_coin(gen)) { break; }
    } else if (split_type == split_type_e::singleton) {
      if (copy_coin(gen)) { break; }
    } else if (split_type == split_type_e::jump) {
      if (jump_coin(gen)) { break; }
    }
  }
  LOG_DEBUG("Splitting took %lu samples", sample_count);
  return {left_dist, right_dist, split_type};
}

} // namespace biogeosim
