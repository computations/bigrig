#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <string>
#include <vector>

namespace bigrig::util {

constexpr size_t factorial_table_size = 11;

constexpr std::array<size_t, factorial_table_size> factorial_table = {
    1,
    1,
    2,
    6,
    24,
    120,
    720,
    5'040,
    40'320,
    362'880,
    3'628'800,
};

constexpr inline auto factorial(uint64_t i) -> size_t {
  if (i < factorial_table_size) { return factorial_table.at(i); }
  size_t f = factorial_table[factorial_table_size - 1];
  for (size_t k = factorial_table_size; k <= i; ++k) {
    f *= static_cast<double>(k);
  }
  return f;
}

constexpr inline auto combinations(uint64_t n, uint64_t i) -> size_t {
  return factorial(n) / (factorial(i) * factorial(n - i));
}

inline auto compute_base26(size_t i) -> std::string {
  size_t length
      = std::max(static_cast<size_t>(std::ceil(
                     std::log(static_cast<double>(i + 1)) / std::log(26.0))),
                 1UL);
  std::string ret;
  ret.reserve(length);

  for (size_t j = 0; j < length; j++) {
    ret += 'a' + (i % 26);
    i   /= 26;
  }

  return ret;
}

inline auto generate_area_names(size_t region_count)
    -> std::vector<std::string> {
  return std::ranges::iota_view{0ul, region_count}
       | std::views::transform(
             [](auto i) -> std::string { return compute_base26(i); })
       | std::ranges::to<std::vector>();
}

template <typename T> constexpr T xorshift(const T &n, int i) {
  return n ^ (n >> i);
}

constexpr auto PHYILP_EXT = ".phy";
constexpr auto NEWICK_EXT = ".nwk";
constexpr auto YAML_EXT   = ".yaml";
constexpr auto JSON_EXT   = ".json";
constexpr auto CSV_EXT    = ".csv";

constexpr size_t VECTOR_INITIAL_RESERVE_COUNT = 8;

} // namespace bigrig::util
