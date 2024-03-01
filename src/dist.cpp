#include "dist.hpp"

#include "split.hpp"

#include <stdexcept>

namespace bigrig {

/**
 * String constrctor. Takes a string and makes a dist. When modifying this
 * function, care should be taken to ensure that the order is correct. The
 * string is indexed in big endian fashion, but the uint64_t in little endian
 * fashion.
 */
dist_t::dist_t(const std::string &dist_string) {
  if (dist_string.size() > 64) {
    throw std::runtime_error{"Tried to make a dist with too many regions"};
  }
  _regions = static_cast<uint16_t>(dist_string.size());

  _dist = 0;

  for (const auto &c : dist_string) {
    _dist  |= (c == '1') ? 1 : 0;
    _dist <<= 1;
  }
  _dist >>= 1; // Need to shift down one since we overshot
}

dist_t &dist_t::operator|=(const dist_t &d) {
  dist_t tmp{*this};
  tmp = tmp | d;
  std::swap(tmp, *this);
  return *this;
}

std::string dist_t::to_str() const {
  std::ostringstream oss;
  oss << *this;
  return oss.str();
}

dist_t make_full_dist(size_t regions) {
  return {(1ul << regions) - 1, static_cast<uint16_t>(regions)};
}

dist_t make_singleton_dist(size_t regions) {
  return {1ul, static_cast<uint16_t>(regions)};
}

} // namespace bigrig
