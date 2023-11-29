#include "dist.hpp"
#include "split.hpp"

#include <stdexcept>

namespace biogeosim {

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

uint64_t valid_region_mask(size_t region_count) {
  uint64_t mask = 0;
  for (size_t i = 0; i < region_count; ++i) { mask |= 1ul << i; }
  return mask;
}

bool valid_dist(dist_t d, const substitution_model_t &model) {
  return valid_dist(d, model.region_count());
}

bool valid_dist(dist_t d, size_t regions) {
  auto mask = valid_region_mask(regions);
  return !static_cast<uint64_t>(d & (~mask));
}

split_type_e determine_split_type(dist_t init_dist,
                                  dist_t left_dist,
                                  dist_t right_dist,
                                  bool   jumps_ok) {
  // Jump case
  if (jumps_ok && (((left_dist | right_dist) ^ init_dist).popcount() == 1)
      && !(left_dist & right_dist)
      && (left_dist.popcount() == 1 || right_dist.popcount() == 1)) {
    return split_type_e::jump;
  }
  // other cases
  if ((left_dist | right_dist) == init_dist) {
    if (left_dist == right_dist && left_dist.popcount() == 1) {
      return split_type_e::singleton;
    }
    if (left_dist.popcount() != 1 && right_dist.popcount() != 1) {
      return split_type_e::invalid;
    }
    if ((left_dist & right_dist).popcount() == 1) {
      return split_type_e::sympatric;
    }
    if ((left_dist & right_dist).popcount() == 0) {
      return split_type_e::allopatric;
    }
  }
  return split_type_e::invalid;
}

std::string split_t::to_nhx_string() const {
  std::ostringstream oss;
  oss << "left-split=" << left << ":"
      << "right-split=" << right << ":";
  oss << "split-type=";

  return oss.str();
}

std::string split_t::type_string() const {
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


} // namespace biogeosim
