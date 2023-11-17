#include "dist.hpp"

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
} // namespace biogeosim
