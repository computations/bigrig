#include <catch2/catch_test_macros.hpp>
#include <dist.hpp>
#include <model.hpp>

TEST_CASE("model init") {
  constexpr size_t REGIONS     = 4;
  bigrig::dist_t   sample_dist = {0b1100, REGIONS};
  SECTION("Basic constructor") { bigrig::biogeo_model_t model; }

  bigrig::biogeo_model_t model;

  SECTION("Setting parameters") {
    model.set_params(1.0, 1.0).set_cladogenesis_params(1.0, 1.0, 1.0, 1.0);

    CHECK(model.extinction_weight(sample_dist) == 2.0);
    CHECK(model.dispersion_weight(sample_dist) == 2.0);
    CHECK(model.total_rate_weight(sample_dist) == 4.0);

    /*
     * y = 2 * dist.full_region_count() == 2
     * j = 2 * dist.empty_region_count() == 6
     */
    CHECK(model.copy_weight(bigrig::make_singleton_dist(REGIONS)) == 2.0);
    CHECK(model.jump_weight(bigrig::make_singleton_dist(REGIONS)) == 6.0);
    CHECK(model.total_singleton_weight(bigrig::make_singleton_dist(REGIONS))
          == 8.0);

    /*
     * j = 2 * sample_dist.empty_region_count() == 4
     * v = 2 * sample_dist.full_region_count() - 2 == 2
     *   - 2 because this is a two region range, and we overcount in this case
     * s = 2 * sample_dist.full_region_count() == 4
     * total = 10
     */
    CHECK(model.jump_weight(sample_dist) == 4.0);
    CHECK(model.allopatry_weight(sample_dist) == 2.0);
    CHECK(model.sympatry_weight(sample_dist) == 4.0);
    CHECK(model.total_nonsingleton_weight(sample_dist) == 10.0);

    /*
     * j = 2 * sample_dist.empty_region_count() == 0
     * v = 2 * sample_dist.full_region_count() == 8
     * s = 2 * sample_dist.full_region_count() == 8
     * total = 16
     */
    CHECK(model.jump_weight(bigrig::make_full_dist(REGIONS)) == 0.0);
    CHECK(model.allopatry_weight(bigrig::make_full_dist(REGIONS)) == 8.0);
    CHECK(model.sympatry_weight(bigrig::make_full_dist(REGIONS)) == 8.0);
    CHECK(model.total_nonsingleton_weight(bigrig::make_full_dist(REGIONS))
          == 16.0);
  }

  SECTION("checks") {
    model.set_cladogenesis_params(0.0, 0.0, 0.0, 0.0);
    CHECK(!model.check_cladogenesis_params_ok(REGIONS));
    CHECK(!model.check_ok(REGIONS));

    model.set_cladogenesis_params(1.0, 1.0, 1.0, 0.0);
    CHECK(model.check_cladogenesis_params_ok(REGIONS));
    CHECK(model.check_ok(REGIONS));

    model.set_cladogenesis_params(1.0, 1.0, 1.0, 1.0);
    CHECK(model.check_cladogenesis_params_ok(REGIONS));
    CHECK(model.check_ok(REGIONS));

    model.set_cladogenesis_params(1.0, 0.0, 0.0, 1.0);
    CHECK(!model.check_cladogenesis_params_ok(REGIONS));
    CHECK(!model.check_ok(REGIONS));
  }
}
