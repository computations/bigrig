#include "adjustment.hpp"

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <dist.hpp>
#include <model.hpp>
#include <pcg_random.hpp>

TEST_CASE("dist operations", "[dist]") {
  constexpr size_t regions = 4;
  bigrig::dist_t   d       = {0b1001, regions};
  bigrig::dist_t   e       = {0b1101, regions};
  SECTION("bitwise ops") {
    CHECK((d ^ e) == bigrig::dist_t{0b0100, regions});
    CHECK((d.region_symmetric_difference(e))
          == bigrig::dist_t{0b0100, regions});
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

    BENCHMARK("access operator") { return e[0]; };
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
    bigrig::dist_t big = {
        0b100'0111'0110'0101'0110'1010'0110'0101'1100'1010'0111'1111'0101'0010'1111'0100,
        63};
    for (size_t i = 0; i < 63; ++i) {
      auto tmp = big.flip_region(i);
      CHECK(tmp != big);
      CHECK((tmp ^ big).full_region_count() == 1);
      CHECK(tmp.region_symmetric_difference(big).full_region_count() == 1);
      CHECK(tmp.region_symmetric_difference_size(big) == 1);
      CHECK(tmp.one_region_off(big));
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

  SECTION("popcount") {
    bigrig::dist_t dist{0b11'0011, 6};
    CHECK(dist.full_region_count() == 4);
    BENCHMARK("popcount speed") { return dist.full_region_count(); };
  }
}

TEST_CASE("spread", "[spread]") {
  constexpr double dis = 1.0;
  constexpr double ext = 1.0;

  uint16_t regions = GENERATE(4, 8, 16, 32);

  bigrig::biogeo_model_t model(dis, ext, true);
  bigrig::dist_t         init_dist = {0b0101, regions};

  pcg64_fast gen(Catch::getSeed());

  auto t = bigrig::spread(init_dist, model, gen);
  CHECK(t.initial_state == init_dist);
  CHECK(t.initial_state != t.final_state);

  BENCHMARK("spread " + std::to_string(regions) + " regions") {
    return bigrig::spread(init_dist, model, gen);
  };
}

TEST_CASE("stats for spread", "[spread][stats]") {
  constexpr size_t regions = 4;

#if D_RIGOROUS
  /* 99.999% confidence that error is less than 0.0001 */
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

  bigrig::biogeo_model_t model(dis, ext, true);

  INFO("dis: " << dis << " ext: " << ext << " dist: " << init_dist);

  double sum   = 0;
  double sumsq = 0;

  for (size_t i = 0; i < iters; ++i) {
    auto t  = bigrig::spread(init_dist, model, gen);
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

TEST_CASE("spread regression") {
  constexpr size_t regions = 4;

#if D_RIGOROUS
  /* 99.999% confidence that error is less than 0.0001 */
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

  bigrig::biogeo_model_t model(dis, ext, true);

  double rej_total = 0;
  double ana_total = 0;

  for (size_t i = 0; i < iters; ++i) {
    auto rej_res  = bigrig::spread_rejection(init_dist, model, gen);
    rej_total    += rej_res.waiting_time;
    auto ana_res  = bigrig::spread_analytic(init_dist, model, gen);
    ana_total    += ana_res.waiting_time;
  }

  double rej_mean = rej_total / iters;
  double ana_mean = ana_total / iters;

  CHECK_THAT(rej_mean - ana_mean, Catch::Matchers::WithinAbs(0, abs_tol));
}

TEST_CASE("spread index chi2 test", "[sample]") {
  constexpr size_t trials = 1e4;
  pcg64_fast       gen(Catch::getSeed());

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

  double dis = GENERATE(0.0, 1.0);
  double ext = GENERATE(0.0, 1.0);

  if (dis == 0.0 && ext == 0.0) { return; }

  INFO("dis:" << dis);
  INFO("ext:" << ext);

  bigrig::biogeo_model_t model(dis, ext, true);

  std::vector<size_t> index_counts;
  index_counts.resize(init_dist.regions());

  for (size_t i = 0; i < trials; ++i) {
    auto res      = bigrig::spread_analytic(init_dist, model, gen);
    auto tmp_dist = (res.final_state ^ res.initial_state);
    index_counts[tmp_dist.last_full_region() - 1] += 1;
  }

  size_t open_regions    = dis == 0.0 ? 0 : init_dist.empty_region_count();
  open_regions          += ext == 0.0 ? 0 : init_dist.full_region_count();
  double expected_count  = static_cast<double>(trials) / open_regions;
  size_t df              = open_regions - 1;

  double chi2 = 0.0;
  for (auto c : index_counts) {
    if (c == 0) { continue; }
    double num  = c - expected_count;
    num        *= num;
    chi2       += num / expected_count;
  }
  REQUIRE(chi2 >= 0);
  CHECK(chi2 < chi2_lut[df]);
}

TEST_CASE("spread regression with adjustment matrix", "[adjust]") {
  constexpr size_t regions = 4;

  /* for these two tests, I don't know what the error rates are, because they
   * seem to be higher due to the adjustment matrix simulation */
  constexpr size_t iters   = 688'609;
  constexpr double abs_tol = 1.0e-2;

  pcg64_fast gen(Catch::getSeed());

  bigrig::adjustment_matrix_t adjust;
  adjust.simulate(regions, gen);

  bigrig::dist_t init_dist = GENERATE(bigrig::dist_t{0b0001, regions},
                                      bigrig::dist_t{0b0100, regions},
                                      bigrig::dist_t{0b1010, regions},
                                      bigrig::dist_t{0b1110, regions},
                                      bigrig::dist_t{0b1111, regions});

  double dis = GENERATE(0.25, 0.66, 1.0, 2.0);
  double ext = GENERATE(0.25, 0.66, 1.0, 2.0);

  bigrig::biogeo_model_t model(dis, ext, true);
  model.set_adjustment_matrix(adjust);

  double rej_total = 0;
  double ana_total = 0;

  for (size_t i = 0; i < iters; ++i) {
    auto rej_res  = bigrig::spread_rejection(init_dist, model, gen);
    rej_total    += rej_res.waiting_time;
    auto ana_res  = bigrig::spread_analytic(init_dist, model, gen);
    ana_total    += ana_res.waiting_time;
  }

  double rej_mean = rej_total / iters;
  double ana_mean = ana_total / iters;

  INFO("rej_mean: " << rej_mean << "\nana_mean: " << ana_mean);

  CHECK_THAT(rej_mean - ana_mean, Catch::Matchers::WithinAbs(0, abs_tol));
}
