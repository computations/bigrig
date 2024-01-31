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

  biogeosim::substitution_model_t model;
  model.set_params(1.0, 1.0)
      .set_cladogenesis_params(1.0, 1.0, 1.0, 0.0)
      .set_region_count(4)
      .set_two_region_duplicity(true);

  SECTION("singleton") {
    biogeosim::dist_t init_dist = {0b1000, regions};
    auto              sp        = biogeosim::split_dist(init_dist, model, gen);
    CHECK(sp.left == sp.right);
    REQUIRE(sp.type == biogeosim::split_type_e::singleton);
  }

  SECTION("allopatry") {
    biogeosim::dist_t init_dist = GENERATE(biogeosim::dist_t{0b1110, regions},
                                           biogeosim::dist_t{0b1100, regions},
                                           biogeosim::dist_t{0b1011, regions},
                                           biogeosim::dist_t{0b1111, regions});
    model.set_cladogenesis_params(0.0, 0.0, 1.0, 0.0);
    auto sp = biogeosim::split_dist(init_dist, model, gen);

    INFO("init: " << init_dist);
    INFO("d1: " << sp.left);
    INFO("d2: " << sp.right);
    CHECK(sp.left);
    CHECK(sp.right);
    REQUIRE(sp.type == biogeosim::split_type_e::allopatric);
    CHECK(sp.left != sp.right);
    CHECK((sp.left.popcount() == 1 || sp.right.popcount() == 1));
    CHECK((sp.left | sp.right) == init_dist);
    CHECK((sp.left & sp.right).popcount() == 0);
  }

  SECTION("sympatry") {
    biogeosim::dist_t init_dist = GENERATE(biogeosim::dist_t{0b1110, regions},
                                           biogeosim::dist_t{0b1100, regions},
                                           biogeosim::dist_t{0b1011, regions},
                                           biogeosim::dist_t{0b1111, regions});
    model.set_cladogenesis_params(0.0, 1.0, 0.0, 0.0);
    auto sp = biogeosim::split_dist(init_dist, model, gen);

    INFO("init: " << init_dist);
    INFO("d1: " << sp.left);
    INFO("d2: " << sp.right);
    CHECK(sp.left);
    CHECK(sp.right);
    REQUIRE(sp.type == biogeosim::split_type_e::sympatric);
    CHECK(sp.left != sp.right);
    CHECK((sp.left.popcount() == 1 || sp.right.popcount() == 1));
    CHECK((sp.left | sp.right) == init_dist);
    CHECK((sp.left & sp.right).popcount()
          == std::min(sp.left.popcount(), sp.right.popcount()));
  }

  SECTION("jump") {
    biogeosim::dist_t init_dist = GENERATE(biogeosim::dist_t{0b1110, regions},
                                           biogeosim::dist_t{0b1100, regions},
                                           biogeosim::dist_t{0b1011, regions});
    model.set_cladogenesis_params(0.0, 0.0, 0.0, 1.0);
    auto sp = biogeosim::split_dist(init_dist, model, gen);

    INFO("init: " << init_dist);
    INFO("d1: " << sp.left);
    INFO("d2: " << sp.right);
    CHECK(sp.left);
    CHECK(sp.right);
    REQUIRE(sp.type == biogeosim::split_type_e::jump);
    CHECK(sp.left != sp.right);
    CHECK((sp.left.popcount() == 1 || sp.right.popcount() == 1));
    CHECK((sp.left | sp.right) != init_dist);
    CHECK((sp.left & sp.right).popcount() == 0);
  }

  SECTION("benchmark") {
    SECTION("generate") {
      biogeosim::dist_t init_dist = GENERATE(biogeosim::dist_t{0b1110, regions},
                                             biogeosim::dist_t{0b1010, regions},
                                             biogeosim::dist_t{0b1111, regions},
                                             biogeosim::dist_t{0b11'1100, 6},
                                             biogeosim::dist_t{0b11'0000, 6},
                                             biogeosim::dist_t{0b11'1111, 6});
      model.set_cladogenesis_params(1.0, 1.0, 1.0, 0.0);
      INFO("init dist:" << init_dist);

      BENCHMARK("split_dist: " + init_dist.to_str()) {
        return biogeosim::split_dist(init_dist, model, gen);
      };
    }

    SECTION("singletons") {
      biogeosim::dist_t init_dist{0b1000, regions};
      BENCHMARK("split_dist singleton") {
        return biogeosim::split_dist(init_dist, model, gen);
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
TEST_CASE("regression") {
  constexpr size_t regions = 4;
  pcg64_fast       gen(Catch::getSeed());

#if D_RIGOROUS
  /* These values for iters and abs_tol ensure with 99.999% confidence that
   * error is less than 0.0001 */
  constexpr size_t iters   = 487'791'396;
  constexpr double abs_tol = 1.0e-4;
#else
  /* These values for iters and abs_tol ensure with 99.999% confidence that error
   * is less than 0.01 */
  constexpr size_t iters   = 48'780;
  constexpr double abs_tol = 1.0e-2;
#endif

  auto keys = {biogeosim::split_type_e::jump,
               biogeosim::split_type_e::sympatric,
               biogeosim::split_type_e::allopatric,
               biogeosim::split_type_e::singleton,
               biogeosim::split_type_e::invalid};

  biogeosim::dist_t init_dist = GENERATE(biogeosim::dist_t{0b1000, regions},
                                         biogeosim::dist_t{0b1110, regions},
                                         biogeosim::dist_t{0b1100, regions},
                                         biogeosim::dist_t{0b1011, regions},
                                         biogeosim::dist_t{0b1111, regions});
  biogeosim::cladogenesis_params_t params
      = GENERATE(biogeosim::cladogenesis_params_t{1.0, 1.0, 1.0, 0.0},
                 biogeosim::cladogenesis_params_t{2.0, 1.0, 1.0, 0.0},
                 biogeosim::cladogenesis_params_t{2.0, 1.0, 2.0, 0.0},
                 biogeosim::cladogenesis_params_t{1.0, 1.0, 1.0, 1.0},
                 biogeosim::cladogenesis_params_t{1.0, 1.0, 2.0, 1.0});

  biogeosim::substitution_model_t model;
  model.set_params(1.0, 1.0).set_region_count(4).set_two_region_duplicity(true);
  INFO("init dist:" << init_dist);
  INFO("model params: " << params.to_debug_string());

  model.set_cladogenesis_params(params);
  std::unordered_map<biogeosim::split_type_e, size_t>
      regression_split_type_counts;
  std::unordered_map<biogeosim::split_type_e, size_t> split_type_counts;
  for (size_t i = 0; i < iters; ++i) {
    auto rej_res
        = biogeosim::split_dist_rejection_method(init_dist, model, gen);
    regression_split_type_counts[rej_res.type] += 1;

    auto new_res = biogeosim::split_dist(init_dist, model, gen);
    split_type_counts[new_res.type] += 1;
  }
  for (const auto &key : keys) {
    double reg_val
        = static_cast<double>(regression_split_type_counts[key]) / iters;
    double new_val = static_cast<double>(split_type_counts[key]) / iters;
    INFO("clado type: " << biogeosim::type_string(key));
    CHECK_THAT(reg_val, Catch::Matchers::WithinAbs(new_val, abs_tol));
  }
}
