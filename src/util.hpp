#pragma once
#include <array>
#include <cstddef>
#include <cstdint>

namespace biogeosim::util {

constexpr size_t factorial_table_size = 11;

constexpr std::array<size_t, factorial_table_size> factorial_table = {
    1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800,
};

constexpr inline auto factorial(uint64_t i) -> size_t {
  if (i < factorial_table_size) {
    return factorial_table.at(i);
  }
  size_t f = factorial_table[factorial_table_size - 1];
  for (size_t k = factorial_table_size; k <= i; ++k) {
    f *= static_cast<double>(k);
  }
  return f;
}

constexpr inline auto combinations(uint64_t n, uint64_t i) -> size_t {
  return factorial(n) / (factorial(i) * factorial(n - i));
}

} // namespace biogeosim::util
