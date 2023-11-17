#include "tree.hpp"

#include "pcg_random.hpp"

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>

TEST_CASE("tree constructor", "[tree]") {
  const std::string tree_str = "((a,b),c);";

  biogeosim::tree_t tree(tree_str);

  CHECK(tree.to_newick() == "((a:0,b:0)1:0,c:0)0:0");
  CHECK(tree.node_count() == 5);
  CHECK(tree.leaf_count() == 3);
}

TEST_CASE("tree sample", "[tree]") {
  constexpr size_t regions = 4;
  constexpr double dis     = 1.0;
  constexpr double ext     = 1.0;

  const std::string tree_str = "((a:1,b:1),c:1);";

  biogeosim::substitution_model_t model(dis, ext, regions);
  biogeosim::dist_t               init_dist = {0b0101, regions};

  biogeosim::tree_t tree(tree_str);

  pcg64_fast gen(Catch::getSeed());

  BENCHMARK("sample on 3 taxa") { tree.sample(init_dist, model, gen); };
}
