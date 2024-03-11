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
  constexpr size_t REGIONS = 4;
  pcg64_fast       gen(Catch::getSeed());

  bigrig::biogeo_model_t model;
  model.set_params(1.0, 1.0)
      .set_cladogenesis_params(1.0, 1.0, 1.0, 0.0)
      .set_region_count(REGIONS)
      .set_two_region_duplicity(false);

  SECTION("singleton") {
    bigrig::dist_t init_dist = {0b1000, REGIONS};
    auto           sp        = bigrig::split_dist(init_dist, model, gen);
    CHECK(sp.left == sp.right);
    REQUIRE(sp.type == bigrig::split_type_e::singleton);
  }

  SECTION("allopatry") {
    bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1110, REGIONS},
                                        bigrig::dist_t{0b1100, REGIONS},
                                        bigrig::dist_t{0b1011, REGIONS},
                                        bigrig::dist_t{0b1111, REGIONS});
    model.set_cladogenesis_params(0.0, 0.0, 1.0, 0.0);
    auto sp = bigrig::split_dist(init_dist, model, gen);

    INFO("init: " << init_dist);
    INFO("d1: " << sp.left);
    INFO("d2: " << sp.right);
    CHECK(sp.left);
    CHECK(sp.right);
    REQUIRE(sp.type == bigrig::split_type_e::allopatric);
    CHECK(sp.left != sp.right);
    CHECK((sp.left.full_region_count() == 1
           || sp.right.full_region_count() == 1));
    CHECK((sp.left | sp.right) == init_dist);
    CHECK((sp.left & sp.right).full_region_count() == 0);
  }

  SECTION("sympatry") {
    bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1110, REGIONS},
                                        bigrig::dist_t{0b1100, REGIONS},
                                        bigrig::dist_t{0b1011, REGIONS},
                                        bigrig::dist_t{0b1111, REGIONS});
    model.set_cladogenesis_params(0.0, 1.0, 0.0, 0.0);
    auto sp = bigrig::split_dist(init_dist, model, gen);

    INFO("init: " << init_dist);
    INFO("d1: " << sp.left);
    INFO("d2: " << sp.right);
    CHECK(sp.left);
    CHECK(sp.right);
    REQUIRE(sp.type == bigrig::split_type_e::sympatric);
    CHECK(sp.left != sp.right);
    CHECK((sp.left.full_region_count() == 1
           || sp.right.full_region_count() == 1));
    CHECK((sp.left | sp.right) == init_dist);
    CHECK(
        (sp.left & sp.right).full_region_count()
        == std::min(sp.left.full_region_count(), sp.right.full_region_count()));
  }

  SECTION("jump") {
    bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1110, REGIONS},
                                        bigrig::dist_t{0b1100, REGIONS},
                                        bigrig::dist_t{0b1011, REGIONS});
    model.set_cladogenesis_params(0.0, 0.0, 0.0, 1.0);
    auto sp = bigrig::split_dist(init_dist, model, gen);

    INFO("init: " << init_dist);
    INFO("d1: " << sp.left);
    INFO("d2: " << sp.right);
    CHECK(sp.left);
    CHECK(sp.right);
    REQUIRE(sp.type == bigrig::split_type_e::jump);
    CHECK(sp.left != sp.right);
    CHECK((sp.left.full_region_count() == 1
           || sp.right.full_region_count() == 1));
    CHECK((sp.left | sp.right) != init_dist);
    CHECK((sp.left & sp.right).full_region_count() == 0);
  }

  SECTION("jump with full range") {
    bigrig::dist_t init_dist{0b1111, REGIONS};

    model.set_cladogenesis_params(0.0, 1.0, 0.0, 1.0);
    auto sp = bigrig::split_dist(init_dist, model, gen);

    INFO("init: " << init_dist);
    INFO("d1: " << sp.left);
    INFO("d2: " << sp.right);
    REQUIRE(sp.type == bigrig::split_type_e::sympatric);
  }

  SECTION("benchmark") {
    SECTION("non-singletons") {
      bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1110, REGIONS},
                                          bigrig::dist_t{0b1010, REGIONS},
                                          bigrig::dist_t{0b1111, REGIONS},
                                          bigrig::dist_t{0b11'1100, 6},
                                          bigrig::dist_t{0b11'0000, 6},
                                          bigrig::dist_t{0b11'1111, 6},
                                          bigrig::dist_t{0b1111'1111, 8});
      SECTION("without jumps") {
        model.set_cladogenesis_params(1.0, 1.0, 1.0, 0.0);
        INFO("init dist:" << init_dist);

        BENCHMARK("split_dist: " + init_dist.to_str()) {
          return bigrig::split_dist(init_dist, model, gen);
        };
      }

      SECTION("with jumps") {
        model.set_cladogenesis_params(1.0, 1.0, 1.0, 1.0);
        INFO("init dist:" << init_dist);

        BENCHMARK("split_dist: " + init_dist.to_str()) {
          return bigrig::split_dist(init_dist, model, gen);
        };
      }
    }

    SECTION("singletons") {
      bigrig::dist_t init_dist{0b1000, REGIONS};
      BENCHMARK("split_dist singleton") {
        return bigrig::split_dist(init_dist, model, gen);
      };
    }
  }
}

TEST_CASE("split stats") {
  pcg64_fast gen(Catch::getSeed());

#if D_RIGOROUS
  /* 99.999% confidence that error is less than 0.0001 */
  constexpr size_t iters   = 1'886'084'219;
  constexpr double abs_tol = 1.0e-4;
#else
  /* 99.999% confidence that error is less than 0.01 */
  constexpr size_t iters   = 188'609;
  constexpr double abs_tol = 1.0e-2;
#endif

  bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1000, 4},
                                      bigrig::dist_t{0b1110, 4},
                                      bigrig::dist_t{0b1100, 4},
                                      bigrig::dist_t{0b1011, 4},
                                      bigrig::dist_t{0b1111, 4},
                                      bigrig::dist_t{0b1111, 4},
                                      bigrig::dist_t{0b10000000, 8},
                                      bigrig::dist_t{0b11000111, 8},
                                      bigrig::dist_t{0b10101010, 8},
                                      bigrig::dist_t{0b11111111, 8});

  bigrig::cladogenesis_params_t params
      = GENERATE(bigrig::cladogenesis_params_t{1.0, 1.0, 1.0, 0.0},
                 bigrig::cladogenesis_params_t{2.0, 1.0, 1.0, 0.0},
                 bigrig::cladogenesis_params_t{2.0, 1.0, 2.0, 0.0},
                 bigrig::cladogenesis_params_t{1.0, 1.0, 1.0, 1.0},
                 bigrig::cladogenesis_params_t{1.0, 1.0, 2.0, 1.0},
                 bigrig::cladogenesis_params_t{1.0, 1.0, 1.0, 2.0},
                 bigrig::cladogenesis_params_t{1.0, 1.2, 1.1, 1.0},
                 bigrig::cladogenesis_params_t{2.0, 1.0, 1.0, 1.0});

  std::unordered_map<bigrig::split_type_e, size_t> split_type_counts;

  bigrig::biogeo_model_t model;
  model.set_params(1.0, 1.0).set_cladogenesis_params(params);

  for (size_t i = 0; i < iters; ++i) {
    auto res                     = bigrig::split_dist(init_dist, model, gen);
    split_type_counts[res.type] += 1;
  }

  if (!model.jumps_ok()) {
    CHECK(split_type_counts[bigrig::split_type_e::jump] == 0);
  }

  if (init_dist.singleton()) {
    auto denom = model.total_singleton_weight(init_dist);

    CHECK_THAT(
        static_cast<double>(split_type_counts[bigrig::split_type_e::jump])
            / iters,
        Catch::Matchers::WithinAbs(model.jump_weight(init_dist) / denom,
                                   abs_tol));

    CHECK_THAT(
        static_cast<double>(split_type_counts[bigrig::split_type_e::singleton])
            / iters,
        Catch::Matchers::WithinAbs(model.copy_weight(init_dist) / denom,
                                   abs_tol));

    CHECK(split_type_counts[bigrig::split_type_e::sympatric] == 0);
    CHECK(split_type_counts[bigrig::split_type_e::allopatric] == 0);
  } else {
    auto denom = model.total_nonsingleton_weight(init_dist);

    CHECK_THAT(
        static_cast<double>(split_type_counts[bigrig::split_type_e::jump])
            / iters,
        Catch::Matchers::WithinAbs(model.jump_weight(init_dist) / denom,
                                   abs_tol));

    CHECK_THAT(
        static_cast<double>(split_type_counts[bigrig::split_type_e::allopatric])
            / iters,
        Catch::Matchers::WithinAbs(model.allopatry_weight(init_dist) / denom,
                                   abs_tol));

    CHECK_THAT(
        static_cast<double>(split_type_counts[bigrig::split_type_e::sympatric])
            / iters,
        Catch::Matchers::WithinAbs(model.sympatry_weight(init_dist) / denom,
                                   abs_tol));

    CHECK(split_type_counts[bigrig::split_type_e::singleton] == 0);
  }

  /* Compute some stats */
}

/*
 * Test the fast, optimized version against a "dumb" but accurate version. The
 * parameter for abs_tol is a compromise between getting results and accuracy.
 * The `GENERATE` expressions will try the following code with each of the
 * listed parameters. So, this code | init_dist X params | = 25 times.
 * In addition, there are 4 checks each.
 */
TEST_CASE("split regression") {
  // constexpr size_t REGIONS = 4;
  pcg64_fast gen(Catch::getSeed());

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

  bigrig::dist_t                init_dist = GENERATE(bigrig::dist_t{0b1000, 4},
                                      bigrig::dist_t{0b1110, 4},
                                      bigrig::dist_t{0b1100, 4},
                                      bigrig::dist_t{0b1011, 4},
                                      bigrig::dist_t{0b1111, 4},
                                      bigrig::dist_t{0b11'0011, 6});
  bigrig::cladogenesis_params_t params
      = GENERATE(bigrig::cladogenesis_params_t{1.0, 1.0, 1.0, 0.0},
                 bigrig::cladogenesis_params_t{2.0, 1.0, 1.0, 0.0},
                 bigrig::cladogenesis_params_t{2.0, 1.0, 2.0, 0.0},
                 bigrig::cladogenesis_params_t{1.0, 1.0, 1.0, 1.0},
                 bigrig::cladogenesis_params_t{1.0, 1.0, 2.0, 1.0});

  bigrig::biogeo_model_t model;
  model.set_params(1.0, 1.0).set_two_region_duplicity(false);
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
    double rejection_value
        = static_cast<double>(regression_split_type_counts[key]) / iters;
    double new_val = static_cast<double>(split_type_counts[key]) / iters;
    INFO("clado type: " << bigrig::type_string(key));
    CHECK_THAT(rejection_value, Catch::Matchers::WithinAbs(new_val, abs_tol));
  }
}
