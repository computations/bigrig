#pragma once

#include "util.hpp"

#include <cstdint>
#include <string>
#include <vector>

class biogeosim_dist_t {
public:
  biogeosim_dist_t() = default;

  biogeosim_dist_t(const biogeosim_dist_t &) = default;
  biogeosim_dist_t(biogeosim_dist_t &&)      = default;

  biogeosim_dist_t &operator=(const biogeosim_dist_t &) = default;
  biogeosim_dist_t &operator=(biogeosim_dist_t &&)      = default;

  biogeosim_dist_t(uint64_t d) : _dist{d} {}

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

  constexpr uint64_t operator[](size_t i) const { return (_dist >> i) & 1ull; }

  biogeosim_dist_t operator^(biogeosim_dist_t d) const {
    return _dist ^ d._dist;
  }
  biogeosim_dist_t operator|(biogeosim_dist_t d) const {
    return _dist | d._dist;
  }
  biogeosim_dist_t operator&(biogeosim_dist_t d) const {
    return _dist & d._dist;
  }

  biogeosim_dist_t operator+(uint64_t d) const { return _dist + d; }

  size_t index(size_t max_areas) {
    size_t skips = compute_skips(_dist, max_areas);
    return _dist - skips;
  }

  biogeosim_dist_t next_dist(uint32_t n) const {
    auto d = *this + 1;
    while (d.popcount() > n) { d = d + 1; }
    return d;
  }

private:
  static auto compute_skips_power_of_2(size_t k, size_t n) -> size_t {
    size_t skips = 0;
    for (size_t i = n + 1; i < k; ++i) { skips += combinations(k - 1, i); }
    return skips;
  }

  static auto compute_skips(size_t i, size_t n) -> size_t {
    constexpr size_t BITS_IN_BYTE = 8;
    size_t           skips        = 0;
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
