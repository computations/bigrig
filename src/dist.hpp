#pragma once

#include "model.hpp"
#include "util.hpp"

#include <cstdint>
#include <logger.hpp>
#include <optional>
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

  dist_t negate_bit(size_t index) const {
    return _dist ^ (_dist & (1ull << index));
  }

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

class transition_t {
public:
  transition_t() = default;
  transition_t(double t, dist_t i, dist_t f)
      : _t{t}, _initial_state{i}, _final_state{f} {}

private:
  double _t;
  dist_t _initial_state;
  dist_t _final_state;
};

std::optional<transition_t> sample(double t, dist_t init_dist,
                                   const substitution_model_t &model,
                                   std::uniform_random_bit_generator auto &gen);

} // namespace biogeosim
