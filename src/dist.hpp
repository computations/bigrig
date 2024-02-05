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

#if __BMI2__
#include <x86intrin.h>
#endif

namespace bigrig {

typedef uint64_t dist_base_t;

class dist_t {
public:
  dist_t() = default;

  dist_t(const std::string &);

  constexpr dist_t(const dist_t &) = default;
  constexpr dist_t(dist_t &&)      = default;

  constexpr dist_t &operator=(const dist_t &) = default;
  constexpr dist_t &operator=(dist_t &&)      = default;

  dist_t &operator|=(const dist_t &d);

  constexpr dist_t(uint64_t d, uint16_t s) : _dist{d}, _regions{s} {}

  constexpr explicit dist_t(uint16_t r) : _dist{0}, _regions{r} {}

  inline constexpr size_t popcount() const {
    return static_cast<size_t>(__builtin_popcountll(_dist));
  }

  inline constexpr size_t unpopcount() const { return regions() - popcount(); }

  inline constexpr bool full() const { return regions() == popcount(); }

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

  constexpr inline dist_t operator^(dist_t d) const {
    return {_dist ^ d._dist, std::max(_regions, d._regions)};
  }

  constexpr inline dist_t operator|(dist_t d) const {
    return {_dist | d._dist, std::max(_regions, d._regions)};
  }

  constexpr inline dist_t operator&(dist_t d) const {
    return {_dist & d._dist, std::max(_regions, d._regions)};
  }

  inline explicit operator uint64_t() const { return _dist; }

  constexpr inline dist_t operator&(uint64_t d) const {
    return {_dist & d, _regions};
  }
  constexpr inline dist_t operator~() const { return {~_dist, _regions}; }

  constexpr inline bool operator==(dist_t d) const {
    return d._dist == _dist && _regions == d._regions;
  }

  constexpr inline bool operator!=(dist_t d) const { return !(d == *this); }

  constexpr inline dist_t negate_bit(size_t index) const {
    return {_dist ^ (1ull << index), _regions};
  }

  constexpr inline dist_t operator+(uint64_t d) const {
    return {_dist + d, _regions};
  }

  constexpr inline explicit operator bool() const {
    return static_cast<bool>(_dist);
  }

  constexpr inline size_t index(size_t max_areas) const {
    size_t skips = compute_skips(_dist, max_areas);
    return _dist - skips;
  }

  constexpr inline size_t set_index(size_t index) const {
#if __BMI2__
    return __builtin_ctz(_pdep_u64(_dist, index));
#else
    size_t tmp_index = 0;
    while (true) {
      if (index == 0 && (*this)[tmp_index] == 1) { break; }
      if ((*this)[tmp_index] == 1) { index -= 1; }
      tmp_index++;
    }
    return tmp_index;
#endif
  }

  constexpr inline size_t unset_index(size_t index) const {
#if __BMI2__
    uint64_t mask = (1ul << regions()) - 1;
    return __builtin_ctz(_pdep_u64(~_dist & mask, index));
#else
    size_t tmp_index = 0;
    while (true) {
      if (index == 0 && (*this)[tmp_index] == 0) { break; }
      if ((*this)[tmp_index] == 0) { index -= 1; }
      tmp_index++;
    }
    return tmp_index;
#endif
  }

  constexpr inline dist_t next_dist(uint32_t n) const {
    auto d = *this + 1;
    while (d.popcount() > n) { d = d + 1; }
    return d;
  }

  std::string to_str() const;

  friend std::ostream &operator<<(std::ostream &os, dist_t dist) {
    for (size_t i = dist._regions; i; --i) { os << dist.bextr(i - 1); }
    return os;
  }

private:
  constexpr inline uint64_t bextr(size_t index) const {
    return (_dist >> index) & 1ull;
  }

  constexpr static auto compute_skips_power_of_2(size_t k, size_t n) -> size_t {
    size_t skips = 0;
    for (size_t i = n + 1; i < k; ++i) {
      skips += bigrig::util::combinations(k - 1, i);
    }
    return skips;
  }

  constexpr static auto compute_skips(size_t i, size_t n) -> size_t {
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
  return sample_analytic(init_dist, model, gen);
}

transition_t sample_rejection(dist_t                                  init_dist,
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

  auto min_ele
      = *std::min_element(rolls.begin(), rolls.end(), [](auto a, auto b) {
          return a.waiting_time < b.waiting_time;
        });
  LOG_DEBUG("waiting time: %f", min_ele.waiting_time);
  return min_ele;
}

transition_t sample_analytic(dist_t                                  init_dist,
                             const substitution_model_t             &model,
                             std::uniform_random_bit_generator auto &gen) {
  auto [dispersion_rate, extinction_rate] = model.rates();

  bool singleton = init_dist.popcount() == 1;

  double dispersion_weight = dispersion_rate * init_dist.unpopcount();
  double extinction_weight = extinction_rate * init_dist.popcount();
  double waiting_time_parameter
      = dispersion_weight + (singleton ? 0.0 : extinction_weight);

  std::exponential_distribution<double> wait_time_distribution(
      waiting_time_parameter);

  double waiting_time = wait_time_distribution(gen);

  std::bernoulli_distribution type_coin(dispersion_weight
                                        / (waiting_time_parameter));

  size_t negate_index = 0;
  if (type_coin(gen)) {
    std::uniform_int_distribution<> index_picker(0, init_dist.unpopcount());
    negate_index = init_dist.unset_index(index_picker(gen));
  } else {
    std::uniform_int_distribution<> index_picker(0, init_dist.popcount());
    negate_index = init_dist.set_index(index_picker(gen));
  }

  LOG_DEBUG("waiting time: %f", waiting_time);

  return {waiting_time, init_dist, init_dist.negate_bit(negate_index)};
}

} // namespace bigrig
