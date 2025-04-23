#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <model.hpp>
#include <node.hpp>
#include <period.hpp>
#include <rng.hpp>

TEST_CASE("node constructors") {
  SECTION("default constructor") { bigrig::node_t n1; }
}

TEST_CASE("simulate tree") {
  constexpr size_t REGIONS = 1;
  constexpr size_t iters   = 1'000;
  pcg64_fast       gen(Catch::getSeed());

  auto rate_params = GENERATE(bigrig::rate_params_t{0.0, 0.0});
  auto clado_params
      = GENERATE(bigrig::cladogenesis_params_t{0.0, 0.0, 1.0, 0.0});
  auto tree_params
      = GENERATE(bigrig::tree_params_t{1.0}, bigrig::tree_params_t{2.0});
  auto ext_allowed = GENERATE(false);
  auto duration    = GENERATE(1.0, 2.0);

  INFO("cladogenesis: " << tree_params.cladogenesis);
  INFO("duration: " << duration);

  bigrig::period_params_t period_params;

  period_params.start      = 0.0;
  period_params.rates      = rate_params;
  period_params.clado      = clado_params;
  period_params.tree       = tree_params;
  period_params.extinction = ext_allowed;

  bigrig::period_list_t periods{{period_params}};

  bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1, REGIONS});

  const double lambda   = period_params.tree->cladogenesis;
  const double t        = duration;
  const double brlen_mu = (2 / lambda) * (std::exp(lambda * t) - 1);
  const double leaf_mu  = 2 * std::exp(lambda * t);

  double brlen_sum    = 0;
  double brlen_sum_sq = 0;
  double leaf_sum     = 0;
  double leaf_sum_sq  = 0;

  /* We need to simulate 2 nodes, because the math assumes a zero length root
   * branch. However, our functions do not generate an immediate speciation
   * event, so we need to fake it by running the simulation twice
   */
  for (size_t i = 0; i < iters; ++i) {
    bigrig::node_t n1;
    n1.simulate_tree(init_dist, duration, periods, gen);

    bigrig::node_t n2;
    n2.simulate_tree(init_dist, duration, periods, gen);

    double cur_brlen_sum  = n1.brlen_sum() + n2.brlen_sum();
    brlen_sum            += cur_brlen_sum;
    brlen_sum_sq         += cur_brlen_sum * cur_brlen_sum;

    auto leaf_count  = n1.leaf_count() + n2.leaf_count();
    leaf_sum        += leaf_count;
    leaf_sum_sq     += leaf_count * leaf_count;
  }

  double brlen_mean = brlen_sum / iters;
  double brlen_std
      = (brlen_sum_sq - brlen_sum * brlen_sum / iters) / (iters - 1);
  double brlen_t_statistic = (brlen_mean - brlen_mu) / (std::sqrt(brlen_std));

  INFO("brlen mean: " << brlen_mean);
  INFO("brlen mu: " << brlen_mu);
  INFO("brlen std: " << brlen_std);
  INFO("brlen t: " << brlen_t_statistic);

  CHECK_THAT(brlen_t_statistic, Catch::Matchers::WithinAbs(0, 4));
  CHECK_THAT(brlen_mean - brlen_mu, Catch::Matchers::WithinAbs(0, 4));

  double leaf_mean = leaf_sum / iters;
  double leaf_std  = (leaf_sum_sq - leaf_sum * leaf_sum / iters) / (iters - 1);
  double leaf_t_statistic = (leaf_mean - leaf_mu) / std::sqrt(leaf_std);

  INFO("leaf mean: " << leaf_sum / iters);
  INFO("leaf mu: " << leaf_mu);
  INFO("leaf std: " << leaf_std);
  INFO("leaf t: " << leaf_t_statistic);

  CHECK_THAT(leaf_t_statistic, Catch::Matchers::WithinAbs(0, 4));
  CHECK_THAT(leaf_mean - leaf_mu, Catch::Matchers::WithinAbs(0, 4));
}
