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
  constexpr size_t  regions = 4;
  biogeosim::dist_t d       = {0b1001, regions};
  biogeosim::dist_t e       = {0b1101, regions};
  SECTION("bitwise ops") {
    CHECK((d ^ e) == biogeosim::dist_t{0b0100, regions});
    CHECK((d | e) == biogeosim::dist_t{0b1101, regions});
    CHECK((d & e) == biogeosim::dist_t{0b1001, regions});
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

  SECTION("log2") {
    CHECK(d.log2() == 4);
    CHECK((biogeosim::dist_t(0b1, 1).log2()) == 1);
  }

  SECTION("valid") {
    CHECK(biogeosim::valid_dist({0b11'0011, 6}, 6));
    CHECK(!biogeosim::valid_dist({0b11'0011, 6}, 5));
  }

  SECTION("addition") {
    auto f = d + 1;
    CHECK(f == biogeosim::dist_t{0b1010, d.regions()});
  }

  SECTION("negate bit") {
    for (size_t i = 0; i < regions; ++i) {
      auto tmp = d.negate_bit(i);
      CHECK(tmp != d);
      CHECK((tmp ^ d).popcount() == 1);
    }
  }

  SECTION("string constructor") {
    CHECK(biogeosim::dist_t("1010") == biogeosim::dist_t{0b1010, 4});
    CHECK(biogeosim::dist_t("0000") == biogeosim::dist_t{0b0000, 4});
    CHECK(biogeosim::dist_t("1011111") == biogeosim::dist_t{0b101'1111, 7});
  }

  SECTION("stream operator") {
    CHECK("1010" == biogeosim::dist_t{0b1010, 4}.to_str());
    CHECK("0000" == biogeosim::dist_t{0b0000, 4}.to_str());
    CHECK(biogeosim::dist_t("1011111").to_str()
          == biogeosim::dist_t{0b101'1111, 7}.to_str());
  }
}

TEST_CASE("sample", "[sample]") {
  constexpr size_t regions = 4;
  constexpr double dis     = 1.0;
  constexpr double ext     = 1.0;

  constexpr double brlen = 1.0;

  biogeosim::substitution_model_t model(dis, ext, regions, true);
  biogeosim::dist_t               init_dist = {0b0101, regions};

  pcg64_fast gen(Catch::getSeed());

  auto t = biogeosim::sample(init_dist, model, gen);
  CHECK(t.initial_state == init_dist);
  CHECK(t.initial_state != t.final_state);

  BENCHMARK("sample") { return biogeosim::sample(init_dist, model, gen); };
}

TEST_CASE("stats for sample", "[sample][stats]") {
  constexpr size_t regions = 4;
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


  pcg64 gen(Catch::getSeed());

  biogeosim::dist_t init_dist = GENERATE(biogeosim::dist_t{0b0001, regions},
                                         biogeosim::dist_t{0b0010, regions},
                                         biogeosim::dist_t{0b0100, regions},
                                         biogeosim::dist_t{0b1000, regions},
                                         biogeosim::dist_t{0b1010, regions},
                                         biogeosim::dist_t{0b1010, regions},
                                         biogeosim::dist_t{0b1110, regions},
                                         biogeosim::dist_t{0b1111, regions});

  double dis = GENERATE(0.25, 0.66, 1.0, 2.0);
  double ext = GENERATE(0.25, 0.66, 1.0, 2.0);

  double average_rate = dis * (regions - init_dist.popcount());

  if (init_dist.popcount() > 1) { average_rate += ext * init_dist.popcount(); }

  double mu    = 1 / (average_rate);
  double sigma = mu * mu;

  biogeosim::substitution_model_t model(dis, ext, regions, true);

  double brlen = GENERATE(0.5, 1.0, 1.5);
  INFO("dis: " << dis << " ext: " << ext << " brlen: " << brlen
               << " dist: " << init_dist);

  double sum   = 0;
  double sumsq = 0;

  for (size_t i = 0; i < iters; ++i) {
    auto t  = biogeosim::sample(init_dist, model, gen);
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
