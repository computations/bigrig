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

  /**
   * Returns the number of occupied (I.E. full) regions.
   */
  inline constexpr size_t full_region_count() const { return popcount(); }

  /**
   * Returns the number of empty regions.
   */
  inline constexpr size_t empty_region_count() const { return unpopcount(); }

  inline constexpr bool singleton() const { return full_region_count() == 1; }

  /**
   * Returns true if all regions are occupied.
   */
  inline constexpr bool full() const { return regions() == popcount(); }

  /**
   * Returns true if all regions are empty. Hypothetically, we should never
   * return a region with this being true but it's good to have anyways.
   */
  inline constexpr bool empty() const { return _dist == 0; }

  /**
   * Returns the highest or last full region (by index) for the current dist.
   */
  inline constexpr size_t last_full_region() const { return log2(); }

  /**
   * Returns the number of regions for the current dist. This is _not_ the
   * number of full regions, but the number of possible regions.
   */
  inline constexpr uint16_t regions() const { return _regions; }

  /**
   * Check if the dist is valid, constrained to a given number of regions, given
   * by the substitution model.
   */
  inline bool valid_dist(const biogeo_model_t &model) {
    return valid_dist(model.region_count());
  }

  /**
   * Check if the dist is valid, constrained to a given number of regions.
   */
  inline bool valid_dist(size_t required_regions) {
    if (required_regions != regions()) { return false; }
    return valid_dist();
  }

  /**
   * We don't check on dist_t construction if the dist is valid, so we have a
   * function to let us know if the dist is valid. Here valid means there are
   * no regions other than the ones allowed.
   */
  inline bool valid_dist() const {
    auto mask = valid_region_mask();
    return !static_cast<uint64_t>(*this & (~mask));
  }

  inline constexpr uint64_t operator[](size_t i) const { return bextr(i); }

  /**
   * Compute the symmetric difference, I.E. the xor of the set.
   */
  constexpr inline dist_t region_symmetric_difference(dist_t d) const {
    return *this ^ d;
  }

  /**
   * Compute the size of the symmetric difference, efficiently.
   */
  constexpr inline size_t region_symmetric_difference_size(dist_t d) const {
    return region_symmetric_difference(d).popcount();
  }

  constexpr inline dist_t operator^(dist_t d) const {
    return {_dist ^ d._dist, std::max(_regions, d._regions)};
  }

  /**
   * Compute the union of the dists.
   */
  constexpr inline dist_t region_union(dist_t d) { return *this | d; }

  constexpr inline dist_t operator|(dist_t d) const {
    return {_dist | d._dist, std::max(_regions, d._regions)};
  }

  constexpr inline dist_t region_intersection(dist_t d) { return *this & d; }

  /**
   * Returns true if the difference between two dists is exactly one region.
   */
  constexpr inline bool one_region_off(dist_t d) {
    return (*this ^ d).popcount() == 1;
  }

  constexpr inline dist_t operator&(dist_t d) const {
    return {_dist & d._dist, std::max(_regions, d._regions)};
  }

  inline explicit operator uint64_t() const { return _dist; }

  constexpr inline dist_t mask(uint64_t d) const { return *this & d; }

  /**
   * And operator for the purposes of masking the region.
   */
  constexpr inline dist_t operator&(uint64_t d) const {
    return {_dist & d, _regions};
  }

  constexpr inline dist_t invert_dist() const { return ~*this; }

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

  /**
   * Computes the index of the dist, given a maximum number of areas. We should
   * think about this dist being the ith dist in a a list ordered as if the dist
   * is a binary number, and this function returns i.
   * */
  constexpr inline size_t index(size_t max_areas) const {
    size_t skips = compute_skips(_dist, max_areas);
    return _dist - skips;
  }

  /**
   * Consider an arbitrary dist, which has some number of full and empty
   * regions. If we want to, for example, pick a random empty region, then we
   * can pick a number from 0 to `empty_region_count()`. However, that number is
   * not an index, and so we need  convert that number into an index. This
   * function does that.
   *
   * Specifically, it computes the index if we want to turn an empty region into
   * a full region.
   */
  constexpr inline size_t set_index(size_t index) const {
    size_t tmp_index = 0;
    while (true) {
      if (index == 0 && bextr(tmp_index)) { break; }
      if (bextr(tmp_index)) { index -= 1; }
      tmp_index++;
    }
    return tmp_index;
  }

  constexpr inline dist_t set_by_count(size_t count) {
    auto index = set_index(count);
    return negate_bit(index);
  }

  /**
   * This is the version of `set_index` if we want to turn a full region into an
   * empty region.
   */
  constexpr inline size_t unset_index(size_t index) const {
    size_t tmp_index = 0;
    while (true) {
      if (index == 0 && !bextr(tmp_index)) { break; }
      if (!bextr(tmp_index)) { index -= 1; }
      tmp_index++;
    }
    return tmp_index;
  }

  constexpr inline dist_t unset_by_count(size_t count) {
    auto index = unset_index(count);
    return negate_bit(index);
  }

  /**
   * Computes the next valid dist, given a max number of regions.
   */
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
  /**
   * Computes the number of set bits. In this case, it is the region count.
   */
  inline constexpr size_t popcount() const {
    return static_cast<size_t>(__builtin_popcountll(_dist));
  }

  /**
   * Computes the number of unset bits, given the region count restriction.
   */
  inline constexpr size_t unpopcount() const {
    return static_cast<size_t>(
        __builtin_popcountll((~_dist) & valid_region_mask()));
  }

  /**
   * Extract a specific bit and set it to the first bit. Specifically, returns 0
   * if the ith bit is unset, and 1 if it is set.
   */
  constexpr inline uint64_t bextr(size_t index) const {
    return (_dist >> index) & 1ull;
  }

  /**
   * Computes a fast log2 that is rounded down to the nearest integer.
   */
  inline constexpr size_t log2() const {
    constexpr size_t BITS_IN_BYTE = 8;
    return sizeof(_dist) * BITS_IN_BYTE - clz();
  }

  /**
   * Count leading zeros
   */
  inline constexpr size_t clz() const {
    return static_cast<size_t>(__builtin_clzll(_dist));
  }

  /**
   * Count leading zeros, with a saftey for when _dist == 0
   */
  inline constexpr size_t clz_min() const {
    return static_cast<size_t>(__builtin_clzll(_dist | 1ull));
  }

  /**
   * Inner function for compute skips.
   */
  constexpr static size_t compute_skips_power_of_2(size_t k, size_t n) {
    size_t skips = 0;
    for (size_t i = n + 1; i < k; ++i) {
      skips += bigrig::util::combinations(k - 1, i);
    }
    return skips;
  }

  /**
   * Compute the number of skipped dists for dist i given a limit n
   */
  constexpr static size_t compute_skips(size_t i, size_t n) {
    constexpr size_t BITS_IN_BYTE = 8;
    size_t           skips        = 0;
    while (i != 0 && n != 0) {
      size_t first_index  = sizeof(i) * BITS_IN_BYTE - __builtin_clzll(i | 1ul);
      skips              += compute_skips_power_of_2(first_index, n);
      n                  -= 1;
      i                  -= 1 << (first_index - 1);
    }
    skips += i;
    return skips;
  }

  constexpr uint64_t valid_region_mask() const {
    uint64_t mask = (1 << regions()) - 1;
    return static_cast<uint64_t>(mask);
  }

  uint64_t _dist;
  uint16_t _regions;
};

/**
 * Data class to store the results of a sample. Records an initial state (as a
 * dist) a final state (as a dist) and the waiting time.
 */
class transition_t {
public:
  transition_t() = default;
  transition_t(double t, dist_t i, dist_t f)
      : waiting_time{t}, initial_state{i}, final_state{f} {}

  double waiting_time = std::numeric_limits<double>::infinity();
  dist_t initial_state;
  dist_t final_state;
};

/**
 * Wrapper function around sample for refactoring. I keep it around in case I
 * want to revert to the rejection method.
 */
transition_t sample(dist_t                                  init_dist,
                    const biogeo_model_t                   &model,
                    std::uniform_random_bit_generator auto &gen) {
  return sample_analytic(init_dist, model, gen);
}

/**
 * Samples a `transition_t` via rejection. That is, it imagines each dist
 * as a independent Poisson process, and rolls them all. It picks the lowest of
 * all the processes. Linear in the number regions. This function exists mostly
 * to use as a check against the results of `sample_analytic`.
 */
transition_t sample_rejection(dist_t                                  init_dist,
                              const biogeo_model_t                   &model,
                              std::uniform_random_bit_generator auto &gen) {
  auto [d, e] = model.rates();

  bool singleton = init_dist.singleton();

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

/**
 * Samples a `transition_t` by combining the independent processes, and only
 * rolling once for the waiting time. There are 2 additional rolls, one for the
 * type, and one for the region.
 *
 * Linear in the number of regions, if there is no BMI2 instruction set.
 */
transition_t sample_analytic(dist_t                                  init_dist,
                             const biogeo_model_t                   &model,
                             std::uniform_random_bit_generator auto &gen) {
  double dispersion_weight = model.dispersion_weight(init_dist);
  double total_weight      = model.total_rate_weight(init_dist);

  std::exponential_distribution<double> wait_time_distribution(total_weight);

  double waiting_time = wait_time_distribution(gen);

  std::bernoulli_distribution type_coin(dispersion_weight / (total_weight));

  auto type = type_coin(gen);

  auto index_picker
      = type
          ? std::uniform_int_distribution<>(0, init_dist.empty_region_count())
          : std::uniform_int_distribution<>(0, init_dist.full_region_count());

  auto new_dist = type ? init_dist.unset_by_count(index_picker(gen))
                       : init_dist.set_by_count(index_picker(gen));

  return {waiting_time, init_dist, new_dist};
}

} // namespace bigrig
