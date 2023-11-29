#include "dist.hpp"

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

split_type_e determine_split_type(dist_t init_dist,
                                  dist_t left_dist,
                                  dist_t right_dist,
                                  bool   jumps_ok);

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
 *  In additoin, Metzke introduced a jump parameter, making it +J. We optionally
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
