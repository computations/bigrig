#pragma once
#include "dist.hpp"
#include "model.hpp"
#include "period.hpp"

namespace bigrig {
/**
 * Wrapper function around sample for refactoring. I keep it around in case I
 * want to revert to the rejection method.
 */
transition_t spread(dist_t                                  init_dist,
                    const biogeo_model_t                   &model,
                    std::uniform_random_bit_generator auto &gen,
                    operation_mode_e mode = operation_mode_e::FAST) {
  if (mode == operation_mode_e::FAST) {
    return spread_analytic(init_dist, model, gen);
  } else if (mode == operation_mode_e::SIM) {
    return spread_rejection(init_dist, model, gen);
  }
  throw std::runtime_error{"Run mode not recognized"};
}

/**
 * Samples a `transition_t` via rejection. That is, it imagines each dist
 * as a independent Poisson process, and rolls them all. It picks the lowest of
 * all the processes. Linear in the number regions. This function exists mostly
 * to use as a check against the results of `sample_analytic`.
 */
transition_t spread_rejection(dist_t                                  init_dist,
                              const biogeo_model_t                   &model,
                              std::uniform_random_bit_generator auto &gen) {
  auto [d, e] = model.rates();

  bool       singleton    = init_dist.singleton();
  const auto region_count = init_dist.regions();

  std::exponential_distribution<double> exp_die(e);

  std::vector<transition_t> rolls(region_count);

  for (size_t i = 0; i < region_count; ++i) {
    if (singleton && init_dist[i]) { continue; }
    double waiting_time
        = init_dist[i]
            ? exp_die(gen)
            : std::exponential_distribution<double>(
                  model.dispersion_weight_for_index(init_dist, i))(gen);
    rolls[i] = transition_t{waiting_time, init_dist, init_dist.flip_region(i)};
  }

  auto min_ele
      = *std::min_element(rolls.begin(), rolls.end(), [](auto a, auto b) {
          return a.waiting_time < b.waiting_time;
        });
  LOG_DEBUG("waiting time: {}", min_ele.waiting_time);
  return min_ele;
}

inline transition_t
spread_flip_region(dist_t                                  init_dist,
                   const biogeo_model_t                   &model,
                   std::uniform_random_bit_generator auto &gen) {
  auto [d, e]         = model.rates();
  double total_weight = model.total_rate_weight(init_dist);

  std::uniform_real_distribution<double> region_dist(0, total_weight);
  double                                 region_roll = region_dist(gen);

  for (size_t i = 0; i < init_dist.regions(); ++i) {
    if (init_dist[i]) {
      if (!init_dist.singleton() || model.extinction_allowed()) {
        region_roll -= e;
      }
    } else {
      region_roll -= model.dispersion_weight_for_index(init_dist, i);
    }
    if (region_roll <= 0) {
      auto new_dist = init_dist.flip_region(i);
      return {std::numeric_limits<double>::infinity(), init_dist, new_dist};
    }
  }
  throw std::runtime_error{"Failed to spread correctly"};
}

/**
 * Samples a `transition_t` by combining the independent processes, and only
 * rolling once for the waiting time. There are 2 additional rolls, one for the
 * type, and one for the region.
 *
 * Linear in the number of regions, due to the `[un]set_by_count` call
 */
transition_t spread_analytic(dist_t                                  init_dist,
                             const biogeo_model_t                   &model,
                             std::uniform_random_bit_generator auto &gen) {
  double total_weight = model.total_rate_weight(init_dist);
  std::exponential_distribution<double> wait_time_distribution(total_weight);

  double waiting_time = wait_time_distribution(gen);

  auto transition         = spread_flip_region(init_dist, model, gen);
  transition.waiting_time = waiting_time;
  return transition;
}

/**
 * Generate transitions for a branch.
 *
 * Simulates transitions by simulating a single transition and subtracting the
 * waiting time from the branch length. If the branch length is still positive
 * or zero, then we record the sample and continue. The process repeats until
 * the branch length is negative.
 */
std::vector<transition_t>
simulate_transitions(dist_t                                  init_dist,
                     const period_list_t                    &periods,
                     std::uniform_random_bit_generator auto &gen,
                     operation_mode_e                        mode) {
  std::vector<transition_t> results;
  results.reserve(util::VECTOR_INITIAL_RESERVE_COUNT);
  double remainder = 0;
  for (const auto &current_period : periods) {
    double brlen = current_period.length();
    while (true) {
      auto r          = spread(init_dist, current_period.model(), gen, mode);
      r.period_index  = current_period.index();
      r.waiting_time += remainder;
      remainder       = 0;
      auto tmp_brlen  = brlen - r.waiting_time;

      if (tmp_brlen < 0.0) {
        remainder = brlen;
        break;
      }

      brlen     = tmp_brlen;
      init_dist = r.final_state;
      results.push_back(r);
    }
  }
  return results;
}
} // namespace bigrig
