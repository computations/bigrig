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

TEST_CASE("adjustment matrix from CSV-style arcs", "[adjust]") {
  const std::vector<std::string> region_names{"a", "b", "c"};
  pcg64_fast                     gen{42};

  SECTION("compact symmetric matrix") {
    bigrig::adjacency_graph_t graph{
        .adjacencies = {{.from = "b", .to = "a", .value = 2.0},
                        {.from = "c", .to = "a", .value = 3.0},
                        {.from = "c", .to = "b", .value = 4.0}},
        .type        = bigrig::adjustment_matrix_symmetry::symmetric};
    bigrig::adjustment_matrix_params_t params{.adjustments = graph};

    bigrig::adjustment_matrix_t matrix{params, region_names, gen};

    CHECK(matrix.is_symmetric());
    CHECK(matrix.get_adjustment(0, 0) == 0.0);
    CHECK(matrix.get_adjustment(0, 1) == 2.0);
    CHECK(matrix.get_adjustment(1, 0) == 2.0);
    CHECK(matrix.get_adjustment(0, 2) == 3.0);
    CHECK(matrix.get_adjustment(2, 0) == 3.0);
    CHECK(matrix.get_adjustment(1, 2) == 4.0);
    CHECK(matrix.get_adjustment(2, 1) == 4.0);
  }

  SECTION("directed matrix") {
    bigrig::adjacency_graph_t graph{
        .adjacencies = {{.from = "a", .to = "b", .value = 1.0},
                        {.from = "b", .to = "a", .value = 2.0},
                        {.from = "a", .to = "c", .value = 3.0},
                        {.from = "c", .to = "a", .value = 4.0},
                        {.from = "b", .to = "c", .value = 5.0},
                        {.from = "c", .to = "b", .value = 6.0}},
        .type        = bigrig::adjustment_matrix_symmetry::nonsymmetric};
    bigrig::adjustment_matrix_params_t params{.adjustments = graph};

    bigrig::adjustment_matrix_t matrix{params, region_names, gen};

    CHECK_FALSE(matrix.is_symmetric());
    CHECK(matrix.get_adjustment(0, 0) == 0.0);
    CHECK(matrix.get_adjustment(0, 1) == 1.0);
    CHECK(matrix.get_adjustment(1, 0) == 2.0);
    CHECK(matrix.get_adjustment(1, 2) == 5.0);
    CHECK(matrix.get_adjustment(2, 1) == 6.0);
  }

  SECTION("full symmetric matrix") {
    bigrig::adjacency_graph_t graph{
        .adjacencies = {{.from = "a", .to = "b", .value = 2.0},
                        {.from = "b", .to = "a", .value = 2.0},
                        {.from = "a", .to = "c", .value = 3.0},
                        {.from = "c", .to = "a", .value = 3.0},
                        {.from = "b", .to = "c", .value = 4.0},
                        {.from = "c", .to = "b", .value = 4.0}},
        .type        = bigrig::adjustment_matrix_symmetry::symmetric};
    bigrig::adjustment_matrix_params_t params{.adjustments = graph};

    bigrig::adjustment_matrix_t matrix{params, region_names, gen};

    CHECK(matrix.is_symmetric());
    CHECK(matrix.get_adjustment(0, 1) == matrix.get_adjustment(1, 0));
    CHECK(matrix.get_adjustment(1, 2) == matrix.get_adjustment(2, 1));
  }
}
