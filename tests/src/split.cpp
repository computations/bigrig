#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <dist.hpp>
#include <model.hpp>
#include <pcg_random.hpp>
#include <split.hpp>

TEST_CASE("splitting", "[sample]") {
  constexpr size_t regions = 4;
  pcg64_fast       gen(Catch::getSeed());

  bigrig::substitution_model_t model;
  model.set_params(1.0, 1.0)
      .set_cladogenesis_params(1.0, 1.0, 1.0, 0.0)
      .set_region_count(4)
      .set_two_region_duplicity(true);

  SECTION("singleton") {
    bigrig::dist_t init_dist = {0b1000, regions};
    auto           sp        = bigrig::split_dist(init_dist, model, gen);
    CHECK(sp.left == sp.right);
    REQUIRE(sp.type == bigrig::split_type_e::singleton);
  }

  SECTION("allopatry") {
    bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1110, regions},
                                        bigrig::dist_t{0b1100, regions},
                                        bigrig::dist_t{0b1011, regions},
                                        bigrig::dist_t{0b1111, regions});
    model.set_cladogenesis_params(0.0, 0.0, 1.0, 0.0);
    auto sp = bigrig::split_dist(init_dist, model, gen);

    INFO("init: " << init_dist);
    INFO("d1: " << sp.left);
    INFO("d2: " << sp.right);
    CHECK(sp.left);
    CHECK(sp.right);
    REQUIRE(sp.type == bigrig::split_type_e::allopatric);
    CHECK(sp.left != sp.right);
    CHECK((sp.left.popcount() == 1 || sp.right.popcount() == 1));
    CHECK((sp.left | sp.right) == init_dist);
    CHECK((sp.left & sp.right).popcount() == 0);
  }

  SECTION("sympatry") {
    bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1110, regions},
                                        bigrig::dist_t{0b1100, regions},
                                        bigrig::dist_t{0b1011, regions},
                                        bigrig::dist_t{0b1111, regions});
    model.set_cladogenesis_params(0.0, 1.0, 0.0, 0.0);
    auto sp = bigrig::split_dist(init_dist, model, gen);

    INFO("init: " << init_dist);
    INFO("d1: " << sp.left);
    INFO("d2: " << sp.right);
    CHECK(sp.left);
    CHECK(sp.right);
    REQUIRE(sp.type == bigrig::split_type_e::sympatric);
    CHECK(sp.left != sp.right);
    CHECK((sp.left.popcount() == 1 || sp.right.popcount() == 1));
    CHECK((sp.left | sp.right) == init_dist);
    CHECK((sp.left & sp.right).popcount()
          == std::min(sp.left.popcount(), sp.right.popcount()));
  }

  SECTION("jump") {
    bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1110, regions},
                                        bigrig::dist_t{0b1100, regions},
                                        bigrig::dist_t{0b1011, regions});
    model.set_cladogenesis_params(0.0, 0.0, 0.0, 1.0);
    auto sp = bigrig::split_dist(init_dist, model, gen);

    INFO("init: " << init_dist);
    INFO("d1: " << sp.left);
    INFO("d2: " << sp.right);
    CHECK(sp.left);
    CHECK(sp.right);
    REQUIRE(sp.type == bigrig::split_type_e::jump);
    CHECK(sp.left != sp.right);
    CHECK((sp.left.popcount() == 1 || sp.right.popcount() == 1));
    CHECK((sp.left | sp.right) != init_dist);
    CHECK((sp.left & sp.right).popcount() == 0);
  }

  SECTION("benchmark") {
    SECTION("generate") {
      bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1110, regions},
                                          bigrig::dist_t{0b1010, regions},
                                          bigrig::dist_t{0b1111, regions},
                                          bigrig::dist_t{0b11'1100, 6},
                                          bigrig::dist_t{0b11'0000, 6},
                                          bigrig::dist_t{0b11'1111, 6});
      model.set_cladogenesis_params(1.0, 1.0, 1.0, 0.0);
      INFO("init dist:" << init_dist);

      BENCHMARK("split_dist: " + init_dist.to_str()) {
        return bigrig::split_dist(init_dist, model, gen);
      };
    }

    SECTION("singletons") {
      bigrig::dist_t init_dist{0b1000, regions};
      BENCHMARK("split_dist singleton") {
        return bigrig::split_dist(init_dist, model, gen);
      };
    }
  }
}

/*
 * Test the fast, optimized version against a "dumb" but accurate version. The
 * parameter for abs_tol is a compromise between getting results and accuracy.
 * The `GENERATE` expressions will try the following code with each of the
 * listed parameters. So, this code | init_dist X params | = 25 times.
 * In addition, there are 4 checks each.
 */
TEST_CASE("split regression") {
  constexpr size_t regions = 4;
  pcg64_fast       gen(Catch::getSeed());

#if D_RIGOROUS
  /* These values for iters and abs_tol ensure with 99.999% confidence that
   * error is less than 0.0001 */
  constexpr size_t iters   = 487'791'396;
  constexpr double abs_tol = 1.0e-4;
#else
  /* These values for iters and abs_tol ensure with 99.999% confidence that
   * error is less than 0.01 */
  constexpr size_t iters   = 48'780;
  constexpr double abs_tol = 1.0e-2;
#endif

  auto keys = {bigrig::split_type_e::jump,
               bigrig::split_type_e::sympatric,
               bigrig::split_type_e::allopatric,
               bigrig::split_type_e::singleton,
               bigrig::split_type_e::invalid};

  bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1000, regions},
                                      bigrig::dist_t{0b1110, regions},
                                      bigrig::dist_t{0b1100, regions},
                                      bigrig::dist_t{0b1011, regions},
                                      bigrig::dist_t{0b1111, regions});
  bigrig::cladogenesis_params_t params
      = GENERATE(bigrig::cladogenesis_params_t{1.0, 1.0, 1.0, 0.0},
                 bigrig::cladogenesis_params_t{2.0, 1.0, 1.0, 0.0},
                 bigrig::cladogenesis_params_t{2.0, 1.0, 2.0, 0.0},
                 bigrig::cladogenesis_params_t{1.0, 1.0, 1.0, 1.0},
                 bigrig::cladogenesis_params_t{1.0, 1.0, 2.0, 1.0});

  bigrig::substitution_model_t model;
  model.set_params(1.0, 1.0).set_region_count(4).set_two_region_duplicity(true);
  INFO("init dist:" << init_dist);
  INFO("model params: " << params.to_debug_string());

  model.set_cladogenesis_params(params);
  std::unordered_map<bigrig::split_type_e, size_t> regression_split_type_counts;
  std::unordered_map<bigrig::split_type_e, size_t> split_type_counts;
  for (size_t i = 0; i < iters; ++i) {
    auto rej_res = bigrig::split_dist_rejection_method(init_dist, model, gen);
    regression_split_type_counts[rej_res.type] += 1;

    auto new_res = bigrig::split_dist(init_dist, model, gen);
    split_type_counts[new_res.type] += 1;
  }
  for (const auto &key : keys) {
    double reg_val
        = static_cast<double>(regression_split_type_counts[key]) / iters;
    double new_val = static_cast<double>(split_type_counts[key]) / iters;
    INFO("clado type: " << bigrig::type_string(key));
    CHECK_THAT(reg_val, Catch::Matchers::WithinAbs(new_val, abs_tol));
  }
}
