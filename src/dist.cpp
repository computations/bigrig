#include "dist.hpp"
#include <random>

namespace biogeosim {
std::optional<transition_t>
sample(double t, dist_t init_dist, const substitution_model_t &model,
       std::uniform_random_bit_generator auto &gen) {

  auto [e, d] = model.rates();

  double sum = model.compute_denominator(init_dist.popcount());
  bool singleton = init_dist.popcount() == 1;

  std::uniform_real_distribution<double> die;
  auto roll = die(gen) * sum;
  LOG_DEBUG("roll: %d", roll);

  /* The order is important, we MUST traverse the same way each time, or else
   * the distribution gets messed up. Additionally, if popcount == 1, then we
   * can't loose a region, since that would constitute an extinction. */
  [[assume(model.region_count() != 0)]];
  for (size_t i = 0; i < model.region_count(); ++i) {
    /* If we have a single region, and _this_ is the active region, we skip
     * this loop iteration */
    if (singleton && init_dist[i]) {
      continue;
    }

    double rate = init_dist[i] ? e : d;
    roll -= rate;

    if (roll < 0.0) {
      LOG_DEBUG("Roll resulted in region %d flipping", i);
      auto next_state = init_dist.negate_bit(i);
      std::exponential_distribution<double> exp(rate);
      double waiting_time = exp(gen);
      return std::optional<transition_t>{waiting_time, init_dist, next_state};
    }
  }
  return {};
}
}; // namespace biogeosim
