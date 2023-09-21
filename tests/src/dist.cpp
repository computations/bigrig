#include "model.hpp"
#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <dist.hpp>
#include <iostream>
#include <math.h>
#include <random>

TEST_CASE("dist operations", "[dist]") {
  constexpr size_t regions = 4;
  biogeosim::dist_t d = 0b1001;
  biogeosim::dist_t e = 0b1101;
  SECTION("bitwise ops") {
    CHECK((d ^ e) == 0b0100);
    CHECK((d | e) == 0b1101);
    CHECK((d & e) == 0b1001);
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
    CHECK((biogeosim::dist_t(0b1).log2()) == 1);
  }

  SECTION("addition") {
    auto f = d + 1;
    CHECK(f == 0b1010);
  }

  SECTION("negate bit") {
    for (size_t i = 0; i < regions; ++i) {
      auto tmp = d.negate_bit(i);
      CHECK(tmp != d);
      CHECK((tmp ^ d).popcount() == 1);
    }
  }
}

TEST_CASE("sample", "[sample]") {
  constexpr size_t regions = 4;
  constexpr double dis = 1.0;
  constexpr double ext = 1.0;

  constexpr double brlen = 1.0;

  biogeosim::substitution_model_t model(dis, ext, regions);
  biogeosim::dist_t init_dist = 0b0101;

  std::minstd_rand gen(Catch::getSeed());

  auto t = biogeosim::sample(init_dist, model, gen);
  CHECK(t.initial_state == init_dist);
  CHECK(t.initial_state != t.final_state);
}

TEST_CASE("stats", "[sample][stats]") {
  const biogeosim::dist_t init_dist = 0b1010;
  std::minstd_rand gen(Catch::getSeed());
  SECTION("e = d = t = 1.0") {
    constexpr size_t regions = 4;
    constexpr double dis = 1.0;
    constexpr double ext = 1.0;

    constexpr double brlen = 1.0;

    constexpr size_t iters = 1e6;
    constexpr double mu = 1.0;
    constexpr double sigma = mu;

    biogeosim::substitution_model_t model(dis, ext, regions);

    double mean = 0;
    double std = 0;

    for (size_t i = 0; i < iters; ++i) {
      auto trans = biogeosim::generate_samples(brlen, model, gen);
      double num = trans.size();
      mean += num;
      std += (num - mu) * (num - mu);
    }

    mean /= iters;
    std /= iters;

    CHECK_THAT(mean - mu, Catch::Matchers::WithinAbs(0.0, 0.01));
    CHECK_THAT(std, Catch::Matchers::WithinAbs(sigma, 0.01));
  }
  SECTION("e = d = 4.0, t = 1.0") {
    constexpr size_t regions = 4;
    constexpr double dis = 4.0;
    constexpr double ext = 4.0;

    constexpr double brlen = 1.0;

    constexpr size_t iters = 1e2;
    constexpr double mu = 4.0;
    constexpr double sigma = 1.0 / (mu * mu) / sqrt(iters);

    biogeosim::substitution_model_t model(dis, ext, regions);

    double sum = 0;
    double sumsq = 0;

    for (size_t i = 0; i < iters; ++i) {
      double count = biogeosim::generate_samples(brlen, model, gen).size();
      sum += count;
      sumsq += count * count;
    }

    INFO("Sum is:" << sum);
    INFO("sumsq is:" << sumsq);

    double mean = sum / iters;
    double std = (sumsq - sum * sum / iters) / (iters - 1);
    double t = (mean - mu) / (std / sqrt(iters));

    CHECK_THAT(mean - mu, Catch::Matchers::WithinAbs(0.0, 0.01));
    CHECK_THAT(std, Catch::Matchers::WithinAbs(sigma, 0.01));
    CHECK_THAT(t, Catch::Matchers::WithinAbs(0, 0.01));
  }
}
