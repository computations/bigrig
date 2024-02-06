#include "dist.hpp"

#include <random>

constexpr size_t VECTOR_INITIAL_RESERVE = 10;

namespace bigrig {
enum class split_type_e { singleton, allopatric, sympatric, jump, invalid };

std::string type_string(const split_type_e &st);

/**
 * Struct to contain the information of a split.
 */
struct split_t {
  dist_t       left;
  dist_t       right;
  dist_t       top;
  split_type_e type;

  std::string to_nhx_string() const;
  std::string to_type_string() const;
};

std::vector<transition_t>
generate_samples(dist_t                                  init_dist,
                 double                                  brlen,
                 const substitution_model_t             &model,
                 std::uniform_random_bit_generator auto &gen) {
  std::vector<transition_t> results;
  results.reserve(VECTOR_INITIAL_RESERVE);
  while (true) {
    auto r  = sample(init_dist, model, gen);
    brlen  -= r.waiting_time;
    if (brlen < 0.0) { return results; }
    LOG_DEBUG("adding transition from %s to %s",
              r.initial_state.to_str().c_str(),
              r.final_state.to_str().c_str());
    init_dist = r.final_state;
    results.push_back(r);
  }
  return results;
}

split_type_e roll_split_type(dist_t                                  init_dist,
                             const substitution_model_t             &model,
                             std::uniform_random_bit_generator auto &gen) {
  if (init_dist.singleton()) {
    double copy_weight = model.copy_weight(init_dist);
    double jump_weight = model.jump_weight(init_dist);

    double sum = copy_weight + jump_weight;

    jump_weight /= sum;

    std::bernoulli_distribution jump_coin(jump_weight);
    if (jump_coin(gen)) { return split_type_e::jump; }
    return split_type_e::singleton;
  }
  double allo_weight = model.allopatry_weight(init_dist);
  double sym_weight  = model.sympatry_weight(init_dist);
  double jump_weight = model.jump_weight(init_dist);

  double sum = allo_weight + sym_weight + jump_weight;

  double roll = std::uniform_real_distribution<double>()(gen) * sum;

  const std::array<const std::pair<double, split_type_e>, 3> options{
      std::pair{allo_weight, split_type_e::allopatric},
      std::pair{sym_weight, split_type_e::sympatric},
      std::pair{jump_weight, split_type_e::jump}};

  for (auto o : options) {
    if (roll <= o.first) { return o.second; }
    roll -= o.first;
  }

  LOG_ERROR("Rolled an invalid split. roll: %f, allo_weight: %f, sym_weight: "
            "%f, jump_weight: %f",
            roll,
            allo_weight,
            sym_weight,
            jump_weight);
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
  if (!model.jumps_ok() && init_dist.singleton()) {
    LOG_DEBUG("Splitting a singleton: %s", init_dist.to_str().c_str());
    return {init_dist, init_dist, init_dist, split_type_e::singleton};
  }

  auto type = roll_split_type(init_dist, model, gen);

  if (type == split_type_e::singleton) {
    return {init_dist, init_dist, init_dist, split_type_e::singleton};
  }

  size_t max_index = (type == split_type_e::jump)
                       ? init_dist.empty_region_count()
                       : init_dist.full_region_count();
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

  return {left_dist, right_dist, init_dist, type};
}

split_type_e
determine_split_type(dist_t init_dist, dist_t left_dist, dist_t right_dist);

split_t
split_dist_rejection_method(dist_t                                  init_dist,
                            const substitution_model_t             &model,
                            std::uniform_random_bit_generator auto &gen) {
  // Singleton case
  if (!model.jumps_ok() && init_dist.singleton()) {
    LOG_DEBUG("Splitting a singleton: %s", init_dist.to_str().c_str());
    return {init_dist, init_dist, init_dist, split_type_e::singleton};
  }
  auto max_dist = (1ul << init_dist.regions()) - 1;
  std::uniform_int_distribution<dist_base_t> dist_gen(1, max_dist);

  auto model_params = model.cladogenesis_params().normalize();

  std::bernoulli_distribution sympatry_coin(model_params.sympatry);

  std::bernoulli_distribution allopatry_coin(model_params.allopatry);

  std::bernoulli_distribution copy_coin(model_params.copy);

  std::bernoulli_distribution jump_coin(model_params.jump);

  dist_t       left_dist, right_dist;
  split_type_e split_type;
  size_t       sample_count = 0;

  while (true) {
    sample_count++;
    left_dist  = dist_t{dist_gen(gen), init_dist.regions()};
    right_dist = dist_t{dist_gen(gen), init_dist.regions()};
    split_type = determine_split_type(init_dist, left_dist, right_dist);

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
  return {left_dist, right_dist, init_dist, split_type};
}
} // namespace bigrig
