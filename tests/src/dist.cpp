#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <dist.hpp>
#include <iostream>
#include <math.h>
#include <model.hpp>
#include <pcg_random.hpp>
#include <random>
#include <sstream>

TEST_CASE("dist operations", "[dist]") {
  constexpr size_t regions = 4;
  bigrig::dist_t   d       = {0b1001, regions};
  bigrig::dist_t   e       = {0b1101, regions};
  SECTION("bitwise ops") {
    CHECK((d ^ e) == bigrig::dist_t{0b0100, regions});
    CHECK((d.region_symmetric_difference(e)) == bigrig::dist_t{0b0100, regions});
    CHECK((d | e) == bigrig::dist_t{0b1101, regions});
    CHECK((d.region_union(e)) == bigrig::dist_t{0b1101, regions});
    CHECK((d & e) == bigrig::dist_t{0b1001, regions});
    CHECK((d.region_intersection(e)) == bigrig::dist_t{0b1001, regions});
  }

  SECTION("access operator") {
    CHECK(d[0] == 1);
    CHECK(d[1] == 0);
    CHECK(d[2] == 0);
    CHECK(d[3] == 1);

    CHECK(e[0] == 1);
    CHECK(e[1] == 0);
    CHECK(e[2] == 1);
    CHECK(e[3] == 1);
  }

  SECTION("last set bit") {
    CHECK(d.last_full_region() == 4);
    CHECK((bigrig::dist_t(0b1, 1).last_full_region()) == 1);
  }

  SECTION("valid") {
    bigrig::dist_t valid{0b11'0011, 6};
    bigrig::dist_t invalid{0b11'0011, 5};
    CHECK(valid.valid_dist(6));
    CHECK(valid.valid_dist());
    CHECK(!invalid.valid_dist(5));
    CHECK(!invalid.valid_dist());
  }

  SECTION("addition") {
    auto f = d + 1;
    CHECK(f == bigrig::dist_t{0b1010, d.regions()});
  }

  SECTION("negate bit") {
    for (size_t i = 0; i < regions; ++i) {
      auto tmp = d.negate_bit(i);
      CHECK(tmp != d);
      CHECK((tmp ^ d).full_region_count() == 1);
      CHECK(tmp.region_symmetric_difference(d).full_region_count() == 1);
      CHECK(tmp.region_symmetric_difference_size(d) == 1);
      CHECK(tmp.one_region_off(d));
    }
  }

  SECTION("string constructor") {
    CHECK(bigrig::dist_t("1010") == bigrig::dist_t{0b1010, 4});
    CHECK(bigrig::dist_t("0000") == bigrig::dist_t{0b0000, 4});
    CHECK(bigrig::dist_t("1011111") == bigrig::dist_t{0b101'1111, 7});
  }

  SECTION("stream operator") {
    CHECK("1010" == bigrig::dist_t{0b1010, 4}.to_str());
    CHECK("0000" == bigrig::dist_t{0b0000, 4}.to_str());
    CHECK(bigrig::dist_t("1011111").to_str()
          == bigrig::dist_t{0b101'1111, 7}.to_str());
  }
}

TEST_CASE("sample", "[sample]") {
  constexpr double dis = 1.0;
  constexpr double ext = 1.0;

  uint16_t regions = GENERATE(4, 8, 16, 32);

  bigrig::biogeo_model_t model(dis, ext, regions, true);
  bigrig::dist_t               init_dist = {0b0101, regions};

  pcg64_fast gen(Catch::getSeed());

  auto t = bigrig::sample(init_dist, model, gen);
  CHECK(t.initial_state == init_dist);
  CHECK(t.initial_state != t.final_state);

  BENCHMARK("sample " + std::to_string(regions) + " regions") {
    return bigrig::sample(init_dist, model, gen);
  };
}

TEST_CASE("stats for sample", "[sample][stats]") {
  constexpr size_t regions    = 4;
  constexpr double expected_t = 4.0;

#if D_RIGOROUS
  /* 99.999% confidence that error is less than 0.001 */
  constexpr size_t iters   = 1'886'084'219;
  constexpr double abs_tol = 1.0e-4;
#else
  /* 99.999% confidence that error is less than 0.01 */
  constexpr size_t iters   = 188'609;
  constexpr double abs_tol = 1.0e-2;
#endif

  pcg64_fast gen(Catch::getSeed());

  bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b0001, regions},
                                      bigrig::dist_t{0b0010, regions},
                                      bigrig::dist_t{0b0100, regions},
                                      bigrig::dist_t{0b1000, regions},
                                      bigrig::dist_t{0b1010, regions},
                                      bigrig::dist_t{0b1110, regions},
                                      bigrig::dist_t{0b1111, regions});

  double dis = GENERATE(0.25, 0.66, 1.0, 2.0);
  double ext = GENERATE(0.25, 0.66, 1.0, 2.0);

  double average_rate = dis * init_dist.empty_region_count();

  if (!init_dist.singleton() && !init_dist.empty()) {
    average_rate += ext * init_dist.full_region_count();
  }

  double mu    = 1 / (average_rate);
  double sigma = mu * mu;

  bigrig::biogeo_model_t model(dis, ext, regions, true);

  INFO("dis: " << dis << " ext: " << ext << " dist: " << init_dist);

  double sum   = 0;
  double sumsq = 0;

  for (size_t i = 0; i < iters; ++i) {
    auto t  = bigrig::sample(init_dist, model, gen);
    sum    += t.waiting_time;
    sumsq  += t.waiting_time * t.waiting_time;
  }

  double mean = sum / iters;
  double std  = (sumsq - sum * sum / iters) / (iters - 1);
  double t    = (mean - mu) / (sqrt(sigma / iters));

  INFO("mean: " << mean << " mu: " << mu << " std: " << std
                << " sigma: " << sigma);

  CHECK_THAT(t, Catch::Matchers::WithinAbs(0, 4));
  CHECK_THAT(mean - mu, Catch::Matchers::WithinAbs(0, abs_tol));
}

TEST_CASE("sample regression") {
  constexpr size_t regions = 4;

#if D_RIGOROUS
  /* 99.999% confidence that error is less than 0.001 */
  constexpr size_t iters   = 1'886'084'219;
  constexpr double abs_tol = 1.0e-4;
#else
  /* 99.999% confidence that error is less than 0.01 */
  constexpr size_t iters   = 188'609;
  constexpr double abs_tol = 1.0e-2;
#endif

  pcg64_fast gen(Catch::getSeed());

  bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b0001, regions},
                                      bigrig::dist_t{0b0100, regions},
                                      bigrig::dist_t{0b1010, regions},
                                      bigrig::dist_t{0b1110, regions},
                                      bigrig::dist_t{0b1111, regions});

  double dis = GENERATE(0.25, 0.66, 1.0, 2.0);
  double ext = GENERATE(0.25, 0.66, 1.0, 2.0);

  bigrig::biogeo_model_t model(dis, ext, regions, true);

  double rej_total = 0;
  double ana_total = 0;

  for (size_t i = 0; i < iters; ++i) {
    auto rej_res  = bigrig::sample_rejection(init_dist, model, gen);
    rej_total    += rej_res.waiting_time;
    auto ana_res  = bigrig::sample_analytic(init_dist, model, gen);
    ana_total    += ana_res.waiting_time;
  }

  double rej_mean = rej_total / iters;
  double ana_mean = ana_total / iters;

  CHECK_THAT(rej_mean - ana_mean, Catch::Matchers::WithinAbs(0, abs_tol));
}
