#pragma once
#include "node.hpp"
#include <corax/corax.hpp>
#include <filesystem>
#include <memory>

namespace biogeosim {
class tree_t {
public:
  tree_t(const std::filesystem::path &tree_filename) {
    auto corax_tree = corax_utree_parse_newick_rooted(tree_filename.c_str());

    _tree = std::make_unique<node_t>(corax_tree->vroot);

    corax_utree_destroy(corax_tree, nullptr);
  }

  tree_t(const tree_t &) = delete;
  tree_t &operator=(const tree_t &) = delete;

  tree_t(tree_t &&) = default;
  constexpr tree_t &operator=(tree_t &&) = default;

private:
  std::unique_ptr<node_t> _tree;
};
} // namespace biogeosim
