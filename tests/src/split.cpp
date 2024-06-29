#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
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
    model.set_cladogenesis_params(1.0, 0.0, 0.0, 0.0);
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

/* This test uses a G-test to test against the expected weights of the model.
 */
TEST_CASE("split g-test") {
  pcg64_fast gen(Catch::getSeed());

#if D_RIGOROUS
  /* These values are computed to ensure a probability of
   *  - Type I error = 0.00001
   *  - Type II error = 0.00001
   *  - And deviation is less than = 0.0001
   */
  constexpr size_t iters = 2'019'696'124;
  constexpr double q     = 18.420680743952584;
#else
  /* These values are computed to ensure a probability of
   *  - Type I error = 0.00001
   *  - Type II error = 0.00001
   *  - And deviation is less than = 0.01
   */
  constexpr size_t iters = 201'970;
  constexpr double q     = 18.420680743952584;
#endif

  auto keys = {bigrig::split_type_e::jump,
               bigrig::split_type_e::sympatric,
               bigrig::split_type_e::allopatric,
               bigrig::split_type_e::singleton,
               bigrig::split_type_e::invalid};

  bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1000, 4},
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
  std::unordered_map<bigrig::split_type_e, size_t> split_type_counts;
  std::unordered_map<bigrig::split_type_e, size_t> expected_split_type_counts;

  auto expected_weights = model.normalized_cladogenesis_params(init_dist);

  expected_split_type_counts[bigrig::split_type_e::jump]
      = expected_weights.jump * iters;
  expected_split_type_counts[bigrig::split_type_e::sympatric]
      = expected_weights.sympatry * iters;
  expected_split_type_counts[bigrig::split_type_e::allopatric]
      = expected_weights.allopatry * iters;
  expected_split_type_counts[bigrig::split_type_e::singleton]
      = expected_weights.copy * iters;
  expected_split_type_counts[bigrig::split_type_e::invalid] = 0;

  for (size_t i = 0; i < iters; ++i) {
    auto new_res = bigrig::split_dist(init_dist, model, gen);
    split_type_counts[new_res.type] += 1;
  }

  double g = 0.0;
  for (const auto &key : keys) {
    if (expected_split_type_counts[key] == 0) { continue; }
    double g_i = static_cast<double>(split_type_counts[key])
               / static_cast<double>(expected_split_type_counts[key]);
    g_i  = std::log(g_i);
    g   += g_i * split_type_counts[key];
  }
  INFO("split_type_counts[jump]: "
       << split_type_counts[bigrig::split_type_e::jump]);
  INFO("regression_split_type_counts[jump]: "
       << expected_split_type_counts[bigrig::split_type_e::jump]);

  INFO("split_type_counts[sympatry]: "
       << split_type_counts[bigrig::split_type_e::sympatric]);
  INFO("regression_split_type_counts[sympatry]: "
       << expected_split_type_counts[bigrig::split_type_e::sympatric]);

  INFO("split_type_counts[allopatric]: "
       << split_type_counts[bigrig::split_type_e::allopatric]);
  INFO("regression_split_type_counts[allopatric]: "
       << expected_split_type_counts[bigrig::split_type_e::allopatric]);

  INFO("split_type_counts[singleton]: "
       << split_type_counts[bigrig::split_type_e::singleton]);
  INFO("regression_split_type_counts[singleton]: "
       << expected_split_type_counts[bigrig::split_type_e::singleton]);

  CHECK(g < q);
}

TEST_CASE("split regression g-test") {
  pcg64_fast gen(Catch::getSeed());

  /* These values are computed to ensure a probability of
   *  - Type I error = 0.00001
   *  - Type II error = 0.00001
   *  - And deviation is less than = 0.01
   */
  constexpr size_t iters = 201'970;
  constexpr double q     = 18.420680743952584;

  auto keys = {bigrig::split_type_e::jump,
               bigrig::split_type_e::sympatric,
               bigrig::split_type_e::allopatric,
               bigrig::split_type_e::singleton,
               bigrig::split_type_e::invalid};

  bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b1000, 4},
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

  double g = 0.0;
  for (const auto &key : keys) {
    if (regression_split_type_counts[key] == 0) { continue; }
    double g_i = static_cast<double>(split_type_counts[key])
               / static_cast<double>(regression_split_type_counts[key]);
    g_i  = std::log(g_i);
    g   += g_i * split_type_counts[key];
  }
  INFO("split_type_counts[jump]: "
       << split_type_counts[bigrig::split_type_e::jump]);
  INFO("regression_split_type_counts[jump]: "
       << regression_split_type_counts[bigrig::split_type_e::jump]);

  INFO("split_type_counts[sympatry]: "
       << split_type_counts[bigrig::split_type_e::sympatric]);
  INFO("regression_split_type_counts[sympatry]: "
       << regression_split_type_counts[bigrig::split_type_e::sympatric]);

  INFO("split_type_counts[allopatric]: "
       << split_type_counts[bigrig::split_type_e::allopatric]);
  INFO("regression_split_type_counts[allopatric]: "
       << regression_split_type_counts[bigrig::split_type_e::allopatric]);

  INFO("split_type_counts[singleton]: "
       << split_type_counts[bigrig::split_type_e::singleton]);
  INFO("regression_split_type_counts[singleton]: "
       << regression_split_type_counts[bigrig::split_type_e::singleton]);

  CHECK(g < q);
}

TEST_CASE("split index chi2 test", "[sample]") {
  pcg64_fast gen(Catch::getSeed());

  /* chi square lut. df -> 99-%ile */
  constexpr double chi2_lut[]{
      6.6348966010212145, 9.21034037197618,   11.344866730144373,
      13.276704135987622, 15.08627246938899,  16.811893829770927,
      18.475306906582357, 20.090235029663233, 21.665994333461924,
      23.209251158954356, 24.724970311318277, 26.216967305535853,
      27.68824961045705,  29.141237740672796, 30.57791416689249,
      31.999926908815176, 33.40866360500461,  34.805305734705065,
      36.19086912927004,  37.56623478662507,  38.93217268351607,
      40.289360437593864, 41.638398118858476, 42.97982013935165,
      44.31410489621915,  45.64168266628317,  46.962942124751436,
      48.27823577031548,  49.58788447289881,  50.89218131151707,
      52.19139483319193,  53.48577183623535,  54.77553976011035,
      56.06090874778906,  57.3420734338592,   58.61921450168706,
      59.89250004508689,  61.1620867636897,   62.4281210161849,
      63.690739751564465, 64.9500713352112,   66.20623628399322,
      67.45934792232582,  68.7095129693454,   69.95683206583814,
      71.20140024831149,  72.44330737654823,  73.68263852010573,
      74.91947430847816,  76.1538912490127,   77.38596201613736,
      78.6157557150025,   79.84333812225145,  81.0687719062971,
      82.29211682919967,  83.51342993198946,  84.73276570506393,
      85.95017624510335,  87.16571139978757,  88.37941890144937,
      89.59134449068712,  90.80153203083871,  92.01002361413214};

  constexpr size_t trials = 1e4;

  bigrig::dist_t init_dist = GENERATE(
      bigrig::dist_t{0b01'1100'1100, 10},
      bigrig::dist_t{0b11'0100'0001, 10},
      bigrig::dist_t{0b00'1111'1011, 10},
      bigrig::dist_t{0b11'0101'0100, 10},
      bigrig::dist_t{
          0b001'1110'1000'0111'1000'1111'1100'1111'1101'1101'1110'1101'1111'0110'1001'0101,
          63},
      bigrig::dist_t{
          0b111'1111'1000'0011'0110'1011'1111'0010'1001'0000'0101'0011'1100'1110'0000'0100,
          63});

  bigrig::biogeo_model_t model;
  model.set_params(1.0, 1.0).set_two_region_duplicity(false);

  std::vector<size_t> index_counts;
  index_counts.resize(init_dist.regions());

  SECTION("sympatry") {
    model.set_cladogenesis_params(0.0, 1.0, 0.0, 0.0);
    for (size_t i = 0; i < trials; ++i) {
      auto   split       = bigrig::split_dist(init_dist, model, gen);
      auto   tmp_dist    = (split.left & split.right);
      size_t split_index = tmp_dist.last_full_region() - 1;
      index_counts[split_index] += 1;
    }

    double chi2 = 0;
    double expected_count
        = trials / static_cast<double>(init_dist.full_region_count());
    size_t df = init_dist.full_region_count() - 1;
    for (auto c : index_counts) {
      if (c == 0) { continue; }
      double num  = c - expected_count;
      num *= num;
      chi2       += num / expected_count;
    }
  }

  SECTION("allopatry") {
    model.set_cladogenesis_params(1.0, 0.0, 0.0, 0.0);
    for (size_t i = 0; i < trials; ++i) {
      auto split = bigrig::split_dist(init_dist, model, gen);

      auto   tmp_dist    = split.left.singleton() ? split.left : split.right;
      size_t split_index = tmp_dist.last_full_region() - 1;
      index_counts[split_index] += 1;
    }
    double chi2 = 0;
    double expected_count
        = trials / static_cast<double>(init_dist.full_region_count());
    size_t df = init_dist.full_region_count() - 1;
    for (auto c : index_counts) {
      if (c == 0) { continue; }
      double num  = c - expected_count;
      num *= num;
      chi2       += num / expected_count;
    }
    CHECK(chi2 < chi2_lut[df]);
  }

  SECTION("jump") {
    model.set_cladogenesis_params(0.0, 0.0, 0.0, 1.0);
    for (size_t i = 0; i < trials; ++i) {
      auto split = bigrig::split_dist(init_dist, model, gen);

      auto   tmp_dist    = split.left.singleton() ? split.left : split.right;
      size_t split_index = tmp_dist.last_full_region() - 1;
      index_counts[split_index] += 1;
    }
    double chi2 = 0;
    double expected_count
        = trials / static_cast<double>(init_dist.empty_region_count());
    size_t df = init_dist.empty_region_count() - 1;
    for (auto c : index_counts) {
      if (c == 0) { continue; }
      double num  = c - expected_count;
      num *= num;
      chi2       += num / expected_count;
    }
    CHECK(chi2 < chi2_lut[df]);
  }
}
