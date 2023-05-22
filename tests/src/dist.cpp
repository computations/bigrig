#include "model.hpp"
#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <dist.hpp>
#include <iostream>
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
  CHECK(t.waiting_time < brlen);
}
