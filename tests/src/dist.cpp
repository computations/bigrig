#include "model.hpp"
#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
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

TEST_CASE("stats for sample", "[sample][stats]") {
  constexpr size_t regions = 4;
  constexpr size_t iters = 1e4;
  std::minstd_rand gen(Catch::getSeed());

  biogeosim::dist_t init_dist =
      GENERATE(0b0001, 0b0010, 0b0100, 0b1000, 0b1010, 0b1010, 0b1110, 0b1111);

  double dis = GENERATE(0.25, 0.66, 1.0, 2.0);
  double ext = GENERATE(0.25, 0.66, 1.0, 2.0);

  double average_rate = dis * (regions - init_dist.popcount());

  if (init_dist.popcount() > 1) {
    average_rate += ext * init_dist.popcount();
  }

  double mu = 1 / (average_rate);
  double sigma = mu * mu;

  biogeosim::substitution_model_t model(dis, ext, regions);

  double brlen = GENERATE(1.0);
  INFO("dis: " << dis << " ext: " << ext << " brlen: " << brlen
               << " dist: " << init_dist);

  double sum = 0;
  double sumsq = 0;

  for (size_t i = 0; i < iters; ++i) {
    auto t = biogeosim::sample(init_dist, model, gen);
    sum += t.waiting_time;
    sumsq += t.waiting_time * t.waiting_time;
  }

  double mean = sum / iters;
  double std = (sumsq - sum * sum / iters) / (iters - 1);
  double t = (mean - mu) / (sqrt(sigma / iters));

  INFO("mean: " << mean << " mu: " << mu << " std: " << std
                << " sigma: " << sigma);

  // CHECK_THAT(mean - mu, Catch::Matchers::WithinAbs(0.0, 0.01));
  // CHECK_THAT(std, Catch::Matchers::WithinAbs(sigma, 0.01));
  CHECK_THAT(t, Catch::Matchers::WithinAbs(0, 4));
}

TEST_CASE("splitting", "[sample]") {
  constexpr size_t regions = 4;
  std::minstd_rand gen(Catch::getSeed());

  biogeosim::substitution_model_t model;
  model.set_params(1.0, 1.0).set_splitting_prob(0.5).set_region_count(4);

  SECTION("singleton") {
    biogeosim::dist_t init_dist = 0b1000;
    auto [d1, d2] = biogeosim::split_dist(init_dist, model, gen);
    CHECK(d1 == d2);
  }
}

/*
TEST_CASE("stats for generate_samples", "[sample][stats]") {
  const biogeosim::dist_t init_dist = 0b1010;
  std::minstd_rand gen(Catch::getSeed());
  constexpr size_t regions = 4;
  constexpr size_t iters = 1e6;

  double dis = GENERATE(1.0, 2.0);
  double ext = GENERATE(1.0, 2.0);

  double brlen = GENERATE(1.0, 2.0);

  INFO("dis: " << dis << " ext: " << ext << " brlen: " << brlen);

  double average_rate =
      ext * init_dist.popcount() + dis * (regions - init_dist.popcount());
  double mu = brlen * average_rate;
  double sigma = average_rate;

  biogeosim::substitution_model_t model(dis, ext, regions);

  double sum = 0;
  double sumsq = 0;

  for (size_t i = 0; i < iters; ++i) {
    double count = biogeosim::generate_samples(brlen, model, gen).size();
    sum += count;
    sumsq += count * count;
  }

  double mean = sum / iters;
  double std = (sumsq - sum * sum / iters) / (iters - 1);
  double t = (mean - mu) / (sqrt(std / iters));
  INFO("mean: " << mean << " mu: " << mu);

  CHECK_THAT(mean - mu, Catch::Matchers::WithinAbs(0.0, 0.01));
  CHECK_THAT(std, Catch::Matchers::WithinAbs(sigma, 0.01));
  CHECK_THAT(t, Catch::Matchers::WithinAbs(0, 2));
}
*/
