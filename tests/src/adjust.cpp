#include <adjustment.hpp>
#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <pcg_random.hpp>

TEST_CASE("adjustment matrix constructor", "[adjust]") {
  constexpr size_t regions = 4;

  pcg64_fast gen(Catch::getSeed());

  bigrig::adjustment_matrix_t am1;

  am1.simulate(regions, gen);

  for (size_t i = 0; i < regions; ++i) {
    for (size_t j = i; j < regions; ++j) {
      CHECK(am1.get_adjustment(i, j) == am1.get_adjustment(j, i));
      CHECK(am1.get_adjustment(i, j) >= 0);
    }
  }

  constexpr double expo = -1.0;

  auto am2{am1};
  am2.apply_exponent(expo);

  for (size_t i = 0; i < regions; ++i) {
    for (size_t j = i + 1; j < regions; ++j) {
      CHECK_THAT(
          am2.get_adjustment(i, j),
          Catch::Matchers::WithinRel(1 / am1.get_adjustment(i, j), 1e-4));
      CHECK(am2.get_adjustment(i, j) >= 0);
    }
  }
}
