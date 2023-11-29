#include "dist.hpp"

#include <random>

namespace biogeosim {
enum class split_type_e { singleton, allopatric, sympatric, jump, invalid };

struct split_t {
  dist_t       left;
  dist_t       right;
  split_type_e type;

  std::string to_nhx_string() const;
  std::string type_string() const;
};

std::vector<transition_t>
generate_samples(dist_t                                  init_dist,
                 double                                  brlen,
                 const substitution_model_t             &model,
                 std::uniform_random_bit_generator auto &gen) {
  std::vector<transition_t> results;
  while (true) {
    auto r  = sample(init_dist, model, gen);
    brlen  -= r.waiting_time;
    if (brlen < 0.0) { return results; }
    LOG_DEBUG("adding transition from %lb to %lb",
              static_cast<uint64_t>(r.initial_state),
              static_cast<uint64_t>(r.final_state));
    init_dist = r.final_state;
    results.push_back(r);
  }
}

split_type_e roll_split_type(dist_t                                  init_dist,
                             const substitution_model_t             &model,
                             std::uniform_random_bit_generator auto &gen) {
  double allo_weight = model.allopatry_weight(init_dist);
  double sym_weight  = model.sympatry_weight(init_dist);
  double jump_weight = model.jump_weight(init_dist);

  double sum = allo_weight + sym_weight + jump_weight;

  allo_weight /= sum;
  sym_weight  /= sum;
  jump_weight /= sum;

  double roll = std::uniform_real_distribution<double>()(gen);

  const std::array<const std::pair<double, split_type_e>, 3> options{
      std::pair{allo_weight, split_type_e::allopatric},
      std::pair{sym_weight, split_type_e::sympatric},
      std::pair{jump_weight, split_type_e::jump}};

  for (auto o : options) {
    if (roll <= o.first) { return o.second; }
    roll -= o.first;
  }
  return split_type_e::invalid;
}

/*
 * There are three types of splitting:
 *  - Singleton
 *  - Allopatric
 *  - Sympatric
 *  Allopatric and Sympatric are not the names used in the original Ree paper,
 *  and they shouldn't be used in user facing descriptions, as they are very
 *  misleading. Regions are large enough that both Allopatry and Sympatry can
 *  occur, but the idea maps well, so I use it internally.
 *
 *  In additoin, Matzke introduced a jump parameter, making it +J. We optionally
 *  support this option.
 */
split_t split_dist(dist_t                                  init_dist,
                   const substitution_model_t             &model,
                   std::uniform_random_bit_generator auto &gen) {
  // Singleton case
  if (!model.jumps_ok() && init_dist.popcount() == 1) {
    LOG_DEBUG("Splitting a singleton: %lb", static_cast<uint64_t>(init_dist));
    return {init_dist, init_dist, split_type_e::singleton};
  }

  auto   type      = roll_split_type(init_dist, model, gen);
  size_t max_index = type == split_type_e::jump ? init_dist.unpopcount()
                                                : init_dist.popcount();
  size_t flipped_index
      = std::uniform_int_distribution<size_t>(0, max_index - 1)(gen);


  if (type == split_type_e::jump) {
    flipped_index = init_dist.unset_index(flipped_index);
  } else {
    flipped_index = init_dist.set_index(flipped_index);
  }

  dist_t left_dist{init_dist};
  dist_t right_dist;

  if (type == split_type_e::allopatric) {
    left_dist = init_dist.negate_bit(flipped_index);
  }
  right_dist = {1ul << flipped_index, init_dist.regions()};

  std::bernoulli_distribution coin(0.5);
  if (coin(gen)) { std::swap(left_dist, right_dist); }

  return {left_dist, right_dist, type};
}

split_type_e determine_split_type(dist_t init_dist,
                                  dist_t left_dist,
                                  dist_t right_dist,
                                  bool   jumps_ok);

split_t
split_dist_rejection_method(dist_t                                  init_dist,
                            const substitution_model_t             &model,
                            std::uniform_random_bit_generator auto &gen) {
  // Singleton case
  if (!model.jumps_ok() && init_dist.popcount() == 1) {
    LOG_DEBUG("Splitting a singleton: %lb", static_cast<uint64_t>(init_dist));
    return {init_dist, init_dist, split_type_e::singleton};
  }
  auto max_dist = (1ul << init_dist.regions()) - 1;
  std::uniform_int_distribution<dist_base_t> dist_gen(1, max_dist);

  std::bernoulli_distribution sympatry_coin(
      model.cladogenesis_params().sympatry);

  std::bernoulli_distribution allopatry_coin(
      model.cladogenesis_params().allopatry);

  std::bernoulli_distribution copy_coin(model.cladogenesis_params().copy);

  std::bernoulli_distribution jump_coin(model.cladogenesis_params().jump);

  dist_t       left_dist, right_dist;
  split_type_e split_type;
  size_t       sample_count = 0;

  while (true) {
    sample_count++;
    left_dist  = dist_t{dist_gen(gen), init_dist.regions()};
    right_dist = dist_t{dist_gen(gen), init_dist.regions()};
    split_type = determine_split_type(
        init_dist, left_dist, right_dist, model.jumps_ok());

    if (split_type == split_type_e::invalid) { continue; }
    if (split_type == split_type_e::sympatric) {
      if (sympatry_coin(gen)) { break; }
    } else if (split_type == split_type_e::allopatric) {
      if (allopatry_coin(gen)) { break; }
    } else if (split_type == split_type_e::singleton) {
      if (copy_coin(gen)) { break; }
    } else if (split_type == split_type_e::jump) {
      if (jump_coin(gen)) { break; }
    }
  }
  LOG_DEBUG("Splitting took %lu samples", sample_count);
  return {left_dist, right_dist, split_type};
}
} // namespace biogeosim