#include "dist.hpp"

constexpr size_t VECTOR_INITIAL_RESERVE = 8;

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
                     double                                  brlen,
                     const biogeo_model_t                   &model,
                     std::uniform_random_bit_generator auto &gen,
                     operation_mode_e                        mode) {
  std::vector<transition_t> results;
  results.reserve(VECTOR_INITIAL_RESERVE);
  while (true) {
    auto r  = spread(init_dist, model, gen, mode);
    brlen  -= r.waiting_time;
    if (brlen < 0.0) { return results; }
    init_dist = r.final_state;
    results.push_back(r);
  }
  return results;
}

/**
 * Simulates a split type based on some model parameters.
 *
 * When using the fast method for computing splits, we simulate a split type
 * first, and then compute the split based on that type. The specific
 * probability for each type depends on the dist to be split. See the model
 * documentation for more information.
 *
 * If the weights for the model are something like s = 0.0, v = 0.0, j = 1.0,
 * and init_dist.full() == true, then this function does some _wacky_ stuff.
 * I've decided to not fix it, as it is a pretty _weird_ parameter set, that
 * makes no sense. I basically need to check that the weight for
 */
split_type_e roll_split_type(dist_t                                  init_dist,
                             const biogeo_model_t                   &model,
                             std::uniform_random_bit_generator auto &gen) {
  if (init_dist.singleton()) {
    double total_weight = model.total_singleton_weight(init_dist);
    double jump_weight  = model.jump_weight(init_dist);

    jump_weight /= total_weight;

    std::bernoulli_distribution jump_coin(jump_weight);
    if (jump_coin(gen)) { return split_type_e::jump; }
    return split_type_e::singleton;
  }
  double allo_weight = model.allopatry_weight(init_dist);
  double sym_weight  = model.sympatry_weight(init_dist);
  double jump_weight = model.jump_weight(init_dist);

  double total_weight = model.total_nonsingleton_weight(init_dist);

  /*
   * There is a function in the standard lib that will do this, but I measured
   * it to be slower than this... So we are just going to stick with this
   * method.
   */
  double roll = std::uniform_real_distribution<double>(0, total_weight)(gen);

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
split_t split_dist(dist_t                                  init_dist,
                   const biogeo_model_t                   &model,
                   std::uniform_random_bit_generator auto &gen,
                   operation_mode_e mode = operation_mode_e::FAST) {
  if (mode == operation_mode_e::FAST) {
    return split_dist_fast(init_dist, model, gen);
  } else if (mode == operation_mode_e::SIM) {
    return split_dist_rejection_method(init_dist, model, gen);
  }
  throw std::runtime_error{"Could not recognize operation mode"};
}

/*
 * There are three types of splitting:
 *  - Singleton
 *  - Allopatric
 *  - Sympatric
 * Allopatric and Sympatric are not the names used in the original Ree paper,
 * and they shouldn't be used in user facing descriptions, as they are very
 * misleading. In the Ree paper, they use the terms "case 1" and "case 2".
 * However, these cases are very directly inspired by the processes of
 * allopatric and sympatric speciation. Additionally, Matzke uses _essentially_
 * these names for the parameters to his cladogenesis model, so what are you
 * going to do.
 *
 * In addition, Matzke introduced a jump parameter, making it +J. We optionally
 * support this option.
 */
split_t split_dist_fast(dist_t                                  init_dist,
                        const biogeo_model_t                   &model,
                        std::uniform_random_bit_generator auto &gen) {
  // Special check for the singleton case
  if (!model.jumps_ok() && init_dist.singleton()) {
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
    left_dist = init_dist.flip_region(flipped_index);
  }
  right_dist = {1ul << flipped_index, init_dist.regions()};

  std::bernoulli_distribution coin(0.5);
  if (coin(gen)) { std::swap(left_dist, right_dist); }

  return {left_dist, right_dist, init_dist, type};
}

split_type_e
determine_split_type(dist_t init_dist, dist_t left_dist, dist_t right_dist);

/**
 * Split a dist via a rejection method.
 *
 * In this type, we generate 2 _completely_ random dists, and then check to see
 * which kind of split this is. If it is a valid type, we return the generated
 * split with probability equal to the corresponding normalized parameter.
 *
 * This function does _not_ support duplicity with allopatric and copy splits.
 * This is a good argument against duplicity. However, I could probably change
 * this function so that it _does_ count duplicity.
 */
split_t
split_dist_rejection_method(dist_t                                  init_dist,
                            const biogeo_model_t                   &model,
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
