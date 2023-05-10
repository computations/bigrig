#pragma once

#include "dist.hpp"
#include "model.hpp"
#include <logger.hpp>
#include <memory>
#include <random>
#include <vector>

class biogeosim_node_t {
public:
  biogeosim_node_t() = default;

  void add_child(const std::shared_ptr<biogeosim_node_t> &n) {
    _children.push_back(n);
  }

  template <typename T>
    requires(std::uniform_random_bit_generator<T>)
  void sample(biogeosim_dist_t                      initial_distribution,
              const biogeosim_substitution_model_t &model,
              T                                    &gen) {
    auto [e, d] = model.rates();

    double sum = model.compute_denominator(initial_distribution.popcount());
    bool   singleton = initial_distribution.popcount() == 1;

    std::uniform_real_distribution<double> die;
    auto                                   roll = die(gen) * sum;
    LOG_DEBUG("roll: %d", roll);

    /* The order is important, we MUST traverse the same way each time, or else
     * the distribution gets messed up. Additionally, if popcount == 1, then we
     * can't loose a region, since that would constitute an extinction. */

    for (size_t i = 0; i < model.region_count(); ++i) {
      /* If we have a single region, and _this_ is the active region, we skip
       * this loop iteration */
      if (singleton && initial_distribution[i]) { continue; }

      roll -= initial_distribution[i] ? e : d;
      if (roll < 0.0) {
        LOG_DEBUG("Roll resulted in region %d flipping", i);
        _state = initial_distribution.negate_bit(i);
        break;
      }
    }
  }

private:
  std::vector<std::shared_ptr<biogeosim_node_t>> _children;
  double                                         _brlen;
  biogeosim_dist_t                               _state;
};
