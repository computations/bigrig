#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cmath>
#include <limits>
#include <period.hpp>
#include <rng.hpp>
#include <util.hpp>

using namespace bigrig;

TEST_CASE("period") {
  period_t p1{};
  SECTION("Null model") {
    REQUIRE(p1.model_ptr() == nullptr);

    period_t p2;
    p2 = p1;

    CHECK(p1.start() == p2.start());
    CHECK(p1.model_ptr() == p2.model_ptr());
    CHECK(p1.length() == p2.length());
    CHECK(p1.index() == p2.index());
  }

  SECTION("Model present") {
    p1.set_model(biogeo_model_t{});

    REQUIRE(p1.model_ptr() != nullptr);

    period_t p2;
    p2 = p1;

    CHECK(p1.start() == p2.start());
    CHECK(p1.model_ptr() != p2.model_ptr());
    CHECK(p1.length() == p2.length());
    CHECK(p1.index() == p2.index());
  }

  SECTION("setters") {
    REQUIRE(p1.model_ptr() == nullptr);
    REQUIRE(p1.length() == 0.0);
    REQUIRE(p1.start() == 0.0);
    REQUIRE(p1.index() == 0);

    p1.set_length(2.0);
    REQUIRE(p1.start() == 0.0);
    REQUIRE(p1.length() == 2.0);

    p1.set_start(1.0);
    REQUIRE(p1.start() == 1.0);
    REQUIRE(p1.length() == 2.0);

    p1.adjust_start(0.0);
    REQUIRE(p1.start() == 0.0);
    REQUIRE(p1.length() == 3.0);

    p1.set_end(1.0);
    REQUIRE(p1.start() == 0.0);
    REQUIRE(p1.length() == 1.0);

    SECTION("clamp") {
      double initial_start = 1.0;
      double initial_end   = 3.0;
      p1.set_start(1.0);
      p1.set_end(3.0);

      auto clamp_start = GENERATE(0.0, 0.5, 1.0, 2.0);
      auto clamp_end   = GENERATE(1.0, 2.0, 3.0, 4.0);

      INFO("clamp start: " << clamp_start << " clamp end: " << clamp_end);
      constexpr double initial_length = 3.0 - 1.0;

      if (clamp_end < clamp_start) {
        REQUIRE_THROWS(p1.clamp(clamp_start, clamp_end));
      } else {
        double expected_start = std::max(initial_start, clamp_start);
        double expected_end   = std::min(initial_end, clamp_end);

        double expected_length = expected_end - expected_start;

        p1.clamp(clamp_start, clamp_end);
        REQUIRE(p1.start() == expected_start);
        REQUIRE(p1.end() == expected_end);
        REQUIRE(p1.length() == expected_length);
      }
    }
  }
}

TEST_CASE("period list") {
  constexpr size_t region_count = 2;
  pcg64_fast       gen{Catch::getSeed()};
  SECTION("default constructor") {
    period_list_t pl1;
    REQUIRE(pl1.size() == 0);

    double get_time = GENERATE(0.0,
                               1.0,
                               2.0,
                               std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN());

    CHECK_THROWS(pl1.get(get_time));
  }

  SECTION("one period") {
    period_list_t pl1{{period_params_t{}},
                      bigrig::util::generate_area_names(region_count),
                      gen};
    REQUIRE(pl1.size() == 1);

    double get_time = GENERATE(0.0,
                               1.0,
                               2.0,
                               std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN());

    if (std::isfinite(get_time)) {
      auto p1 = pl1.get(get_time);
      CHECK(p1.start() == 0.0);
      CHECK(p1.end() == std::numeric_limits<double>::infinity());
    } else {
      CHECK_THROWS(pl1.get(get_time));
    }
  }

  SECTION("two periods") {
    period_params_t pp1{};
    period_params_t pp2{};

    pp1.start = 0.0;
    pp2.start = 1.0;

    period_list_t pl1{
        {pp1, pp2}, bigrig::util::generate_area_names(region_count), gen};

    REQUIRE(pl1.size() == 2);

    CHECK(pl1.get(0.0).model_ptr() != pl1.get(1.0).model_ptr());

    double get_time = GENERATE(0.0,
                               1.0,
                               2.0,
                               std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN());

    if (std::isfinite(get_time)) {
      auto p1 = pl1.get(get_time);

      CHECK(p1.start() <= get_time);
      CHECK(p1.end() > get_time);

    } else {
      CHECK_THROWS(pl1.get(get_time));
    }

    auto test_filter = [](const auto &p) { return p.start() > 0.5; };

    for (auto p : pl1 | std::ranges::views::filter(test_filter)) {
      CHECK(p.start() > 0.5);
    }
  }
}
