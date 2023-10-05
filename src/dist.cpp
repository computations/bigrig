#include "dist.hpp"

namespace biogeosim {

uint64_t valid_region_mask(size_t region_count) {
  uint64_t mask = 0;
  for (size_t i = 0; i < region_count; ++i) {
    mask |= 1ul << i;
  }
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
