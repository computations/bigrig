#include "dist.hpp"

#include "split.hpp"

#include <stdexcept>

namespace bigrig {

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

split_type_e
determine_split_type(dist_t init_dist, dist_t left_dist, dist_t right_dist) {
  if (left_dist.full_region_count() > right_dist.full_region_count()) {
    std::swap(left_dist, right_dist);
  }

  if ((left_dist | right_dist) == init_dist) {
    if (left_dist == right_dist && left_dist.full_region_count() == 1) {
      return split_type_e::singleton;
    }
    if (!left_dist.singleton() && !right_dist.singleton()) {
      return split_type_e::invalid;
    }
    if ((left_dist & right_dist).full_region_count() == 1) {
      return split_type_e::sympatric;
    }
    if ((left_dist & right_dist).full_region_count() == 0) {
      return split_type_e::allopatric;
    }
  } else if (left_dist.singleton() && right_dist == init_dist
             && (left_dist | init_dist).full_region_count()
                    == init_dist.full_region_count() + 1) {
    return split_type_e::jump;
  }
  return split_type_e::invalid;
}

} // namespace bigrig
