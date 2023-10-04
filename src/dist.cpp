#include "dist.hpp"

namespace biogeosim {
bool valid_dist(dist_t d, const substitution_model_t &model) {
  auto mask = model.valid_mask();
  return !static_cast<uint64_t>(d & (~mask));
}
} // namespace biogeosim
