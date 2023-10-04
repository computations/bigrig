#pragma once

#include "model.hpp"
#include "util.hpp"

#include <algorithm>
#include <cstdint>
#include <format>
#include <limits>
#include <logger.hpp>
#include <random>
#include <string>
#include <vector>

namespace biogeosim {

class dist_t {
public:
  dist_t() = default;

  dist_t(const dist_t &) = default;
  dist_t(dist_t &&) = default;

  dist_t &operator=(const dist_t &) = default;
  dist_t &operator=(dist_t &&) = default;

  dist_t(uint64_t d) : _dist{d} {}

  constexpr size_t popcount() const {
    return static_cast<size_t>(__builtin_popcountll(_dist));
  }

  constexpr size_t clz() const {
    return static_cast<size_t>(__builtin_clzll(_dist));
  }

  constexpr size_t clz_min() const {
    return static_cast<size_t>(__builtin_clzll(_dist | 1ull));
  }

  constexpr size_t log2() const {
    constexpr size_t BITS_IN_BYTE = 8;
    return sizeof(_dist) * BITS_IN_BYTE - clz();
  }

  constexpr uint64_t operator[](size_t i) const { return bextr(i); }

  dist_t operator^(dist_t d) const { return _dist ^ d._dist; }
  dist_t operator|(dist_t d) const { return _dist | d._dist; }
  dist_t operator&(dist_t d) const { return _dist & d._dist; }

  explicit operator uint64_t() const { return _dist; }

  dist_t operator&(uint64_t d) const { return _dist & d; }
  dist_t operator~() const { return ~_dist; }

  bool operator==(dist_t d) const { return d._dist == _dist; }

  dist_t negate_bit(size_t index) const { return _dist ^ (1ull << index); }

  dist_t operator+(uint64_t d) const { return _dist + d; }

  size_t index(size_t max_areas) const {
    size_t skips = compute_skips(_dist, max_areas);
    return _dist - skips;
  }

  dist_t next_dist(uint32_t n) const {
    auto d = *this + 1;
    while (d.popcount() > n) {
      d = d + 1;
    }
    return d;
  }

  friend std::ostream &operator<<(std::ostream &os, dist_t dist) {
    os << std::format("{:b}", dist._dist);
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
    size_t skips = 0;
    while (i != 0 && n != 0) {
      size_t first_index = sizeof(i) * BITS_IN_BYTE - __builtin_clzll(i | 1);
      skips += compute_skips_power_of_2(first_index, n);
      n -= 1;
      i -= 1 << (first_index - 1);
    }
    skips += i;
    return skips;
  }

  uint64_t _dist;
};

bool valid_dist(dist_t d, const substitution_model_t &model);

class transition_t {
public:
  transition_t() = default;
  transition_t(double t, dist_t i, dist_t f)
      : waiting_time{t}, initial_state{i}, final_state{f} {}

  double waiting_time = std::numeric_limits<double>::infinity();
  dist_t initial_state;
  dist_t final_state;
};

transition_t sample(dist_t init_dist, const substitution_model_t &model,
                    std::uniform_random_bit_generator auto &gen) {

  auto [d, e] = model.rates();

  bool singleton = init_dist.popcount() == 1;

  std::exponential_distribution<double> exp_die(e);
  std::exponential_distribution<double> dis_die(d);

  std::vector<transition_t> rolls(model.region_count());

  for (size_t i = 0; i < model.region_count(); ++i) {
    if (init_dist[i] && singleton) {
      continue;
    }
    double waiting_time = init_dist[i] ? exp_die(gen) : dis_die(gen);
    rolls[i] = transition_t{waiting_time, init_dist, init_dist.negate_bit(i)};
  }

  return *std::min_element(rolls.begin(), rolls.end(), [](auto a, auto b) {
    return a.waiting_time < b.waiting_time;
  });
}

std::vector<transition_t>
generate_samples(double brlen, const substitution_model_t &model,
                 std::uniform_random_bit_generator auto &gen) {
  std::vector<transition_t> results;
  while (true) {
    auto r = sample(brlen, model, gen);
    brlen -= r.waiting_time;
    if (brlen < 0.0) {
      return results;
    }
    LOG_DEBUG("adding transition from %b to %b", r.initial_state,
              r.final_state);
    results.push_back(r);
  }
}

/*
 * There are three types of splitting:
 *  - Singleton
 *  - Allopatric
 *  - Sympatric
 *  Allopatric and Sympatric are not the names used in the original Ree paper,
 *  and they shouldn't be used in user facing descriptions, as they are very
 *  misleading. Regions are large enough that both Allopatry and Sympatry can
 *  occur, but the idea maps well, so I use it internally.
 */
std::pair<dist_t, dist_t>
split_dist(dist_t init_dist, const substitution_model_t &model,
           std::uniform_random_bit_generator auto &gen) {
  // Singleton case
  if (init_dist.popcount() == 1) {
    LOG_DEBUG("Splitting a singleton: %b", init_dist);
    return {init_dist, init_dist};
  }

  std::bernoulli_distribution coin(model.splitting_prob());
  std::uniform_int_distribution<size_t> die(0, init_dist.popcount() - 1);

  size_t flipped_index = die(gen);
  for (size_t i = 0; i < flipped_index + 1; ++i) {
    if (!init_dist[i]) {
      flipped_index++;
    }
  }

  dist_t left_dist = init_dist;
  dist_t right_dist = left_dist.negate_bit(flipped_index);
  /* In the allopatric case, we need to remove the index from the init dist. */
  if (coin(gen)) {
    LOG_DEBUG("%s", "Allopatric Split");
    init_dist = init_dist.negate_bit(flipped_index);
  }

  std::bernoulli_distribution left_or_right(0.5);
  if (left_or_right(gen)) {
    std::swap(left_dist, right_dist);
  }

  LOG_DEBUG("Spilt %b into: %b, %b", init_dist, left_dist, right_dist);
  return {left_dist, right_dist};
}

} // namespace biogeosim
