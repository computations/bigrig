#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <model.hpp>
#include <node.hpp>
#include <period.hpp>
#include <rng.hpp>

constexpr double compute_std(double sum, double sum_sq, size_t iters) {
  return (sum_sq - sum * sum / iters) / (iters - 1);
}

constexpr double compute_t_stat(double mean, double mu, double std) {
  return (mean - mu) / std::sqrt(std);
}

constexpr bigrig::period_params_t
make_period_params(bigrig::rate_params_t         rate_params,
                   bigrig::cladogenesis_params_t clado_params,
                   bigrig::tree_params_t         tree_params,
                   bool                          ext_allowed) {
  bigrig::period_params_t period_params;

  period_params.start      = 0.0;
  period_params.rates      = rate_params;
  period_params.clado      = clado_params;
  period_params.tree       = tree_params;
  period_params.extinction = ext_allowed;
  return period_params;
}

bigrig::period_list_t
make_single_period(bigrig::rate_params_t         rate_params,
                   bigrig::cladogenesis_params_t clado_params,
                   bigrig::tree_params_t         tree_params,
                   bool                          ext_allowed) {
  return {{make_period_params(
      rate_params, clado_params, tree_params, ext_allowed)}};
}

TEST_CASE("node constructors") {
  SECTION("default constructor") { bigrig::node_t n1; }
}

TEST_CASE("simulate tree") {
  pcg64_fast       gen(Catch::getSeed());
  constexpr size_t iters = 10'000;

  SECTION("pure birth") {
    constexpr size_t REGIONS = 1;

    auto rate_params = GENERATE(bigrig::rate_params_t{0.0, 0.0});
    auto clado_params
        = GENERATE(bigrig::cladogenesis_params_t{0.0, 0.0, 1.0, 0.0});
    auto tree_params
        = GENERATE(bigrig::tree_params_t{1.0}, bigrig::tree_params_t{2.0});
    auto ext_allowed = GENERATE(false);
    auto duration    = GENERATE(1.0, 2.0);

    INFO("cladogenesis: " << tree_params.cladogenesis);
    INFO("duration: " << duration);

    bigrig::period_list_t periods = make_single_period(
        rate_params, clado_params, tree_params, ext_allowed);

    bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1, REGIONS});

    const double lambda   = tree_params.cladogenesis;
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

    double brlen_mean        = brlen_sum / iters;
    double brlen_std         = compute_std(brlen_sum, brlen_sum_sq, iters);
    double brlen_t_statistic = compute_t_stat(brlen_mean, brlen_mu, brlen_std);

    INFO("brlen mean: " << brlen_mean);
    INFO("brlen mu: " << brlen_mu);
    INFO("brlen std: " << brlen_std);
    INFO("brlen t: " << brlen_t_statistic);

    CHECK_THAT(brlen_t_statistic, Catch::Matchers::WithinAbs(0, 4));
    CHECK_THAT(brlen_mean - brlen_mu, Catch::Matchers::WithinAbs(0, 4));

    double leaf_mean        = leaf_sum / iters;
    double leaf_std         = compute_std(leaf_sum, leaf_sum_sq, iters);
    double leaf_t_statistic = compute_t_stat(leaf_mean, leaf_mu, leaf_std);

    INFO("leaf mean: " << leaf_sum / iters);
    INFO("leaf mu: " << leaf_mu);
    INFO("leaf std: " << leaf_std);
    INFO("leaf t: " << leaf_t_statistic);

    CHECK_THAT(leaf_t_statistic, Catch::Matchers::WithinAbs(0, 4));
    CHECK_THAT(leaf_mean - leaf_mu, Catch::Matchers::WithinAbs(0, 4));
  }

  SECTION("birth death") {
    constexpr size_t REGIONS = 1;

    auto rate_params = GENERATE(bigrig::rate_params_t{0.0, 1.0},
                                bigrig::rate_params_t{0.0, 2.0});
    auto clado_params
        = GENERATE(bigrig::cladogenesis_params_t{0.0, 0.0, 1.0, 0.0});
    auto tree_params
        = GENERATE(bigrig::tree_params_t{1.0}, bigrig::tree_params_t{2.0});
    auto ext_allowed = GENERATE(true);
    auto duration    = GENERATE(0.5, 1.0, 2.0);

    bigrig::period_list_t periods = make_single_period(
        rate_params, clado_params, tree_params, ext_allowed);

    INFO("cladogenesis: " << tree_params.cladogenesis);
    INFO("extinction: " << rate_params.ext);
    INFO("duration: " << duration);

    const double net_lambda = tree_params.cladogenesis - rate_params.ext;
    const double rate_ratio = tree_params.cladogenesis / rate_params.ext;
    const double t          = duration;
    const double leaf_mu    = 2 * std::exp(net_lambda * t);
    const double rate_time  = net_lambda * t;
    const double f_rho      = (rate_ratio * std::exp(rate_time) - 1)
                       / ((rate_ratio - 1) * std::exp(rate_time));
    const double brlen_mu
        = 2 * (std::exp(net_lambda * t) / rate_params.ext) * std::log(f_rho);

    INFO("f_rho: " << f_rho);

    bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1, REGIONS});

    double leaf_sum    = 0;
    double leaf_sum_sq = 0;

    double brlen_sum    = 0;
    double brlen_sum_sq = 0;

    /* We need to simulate 2 nodes, because the math assumes a zero length root
     * branch. However, our functions do not generate an immediate speciation
     * event, so we need to fake it by running the simulation twice
     */
    for (size_t i = 0; i < iters; ++i) {
      bigrig::node_t n1;
      n1.simulate_tree(init_dist, duration, periods, gen);
      n1.assign_abs_time_root();

      bigrig::node_t n2;
      n2.simulate_tree(init_dist, duration, periods, gen);
      n2.assign_abs_time_root();

      double cur_brlen_sum
          = n1.reconstructed_brlen_sum() + n2.reconstructed_brlen_sum();
      brlen_sum    += cur_brlen_sum;
      brlen_sum_sq += cur_brlen_sum * cur_brlen_sum;

      auto leaf_count
          = n1.reconstructed_leaf_count() + n2.reconstructed_leaf_count();
      leaf_sum    += leaf_count;
      leaf_sum_sq += leaf_count * leaf_count;
    }

    INFO("leaf_sum: " << leaf_sum);
    INFO("leaf_sum_sq: " << leaf_sum_sq);
    INFO("brlen_sum: " << brlen_sum);
    INFO("brlen_sum_sq: " << brlen_sum_sq);

    double leaf_mean        = leaf_sum / iters;
    double leaf_std         = compute_std(leaf_sum, leaf_sum_sq, iters);
    double leaf_t_statistic = compute_t_stat(leaf_mean, leaf_mu, leaf_std);

    CHECK_THAT(leaf_t_statistic, Catch::Matchers::WithinAbs(0, 4));

    /* The epxected values for branch lengths breaks down at rate ratio == 1,
    * because we would have to evaluate a 0/0 expression. So, we skip that part,
    * and simply check all the other cases
    */
    if (rate_ratio != 1) {
      double brlen_mean = brlen_sum / iters;
      double brlen_std  = compute_std(brlen_sum, brlen_sum_sq, iters);
      double brlen_t_statistic
          = compute_t_stat(brlen_mean, brlen_mu, brlen_std);

      INFO("brlen_mean: " << brlen_mean);
      INFO("brlen_mu: " << brlen_mu);
      INFO("brlen_std: " << brlen_std);

      CHECK_THAT(brlen_t_statistic, Catch::Matchers::WithinAbs(0, 4));
    }
  }
}
