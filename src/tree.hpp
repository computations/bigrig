#pragma once
#include "dist.hpp"
#include "exceptions.hpp"
#include "iterator.hpp"
#include "model.hpp"
#include "node.hpp"

#include <corax/corax.hpp>
#include <filesystem>
#include <memory>
#include <optional>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>

namespace bigrig {
/**
 * Top level class around `node_t`.
 *
 * Uses the coraxlib newick parser to parse the tree, but then immediatly
 * converts the tree into `node_t`s.
 */
class tree_t {
public:
  explicit tree_t(const std::filesystem::path &tree_filename);

  explicit tree_t(const std::string &tree_str);

  tree_t(const tree_t &)            = delete;
  tree_t &operator=(const tree_t &) = delete;

  tree_t(tree_t &&)            = default;
  tree_t &operator=(tree_t &&) = default;

  /**
   * Simulate the whole tree from an initial dist.
   */
  void simulate(dist_t                                  initial_distribution,
                const biogeo_model_t                   &model,
                std::uniform_random_bit_generator auto &gen) {
    if (!initial_distribution.valid_dist(model)) {
      throw invalid_dist{"Invalid dist provided as a start dist"};
    }
    LOG_DEBUG("Starting sample with init dist = %lb",
              static_cast<uint64_t>(initial_distribution));
    _tree->simulate(initial_distribution, model, gen);
  }

  std::optional<dist_t> get_dist_by_string_id(const std::string &key) const;

  std::string to_newick() const;

  std::string
  to_newick(std::function<void(std::ostream &, const node_t &)> cb) const;

  std::string to_phylip_body() const;

  std::string to_phylip_body_extended() const;

  std::ostream &to_phylip_body(std::ostream &os, bool all = false) const;

  size_t node_count() const;
  size_t leaf_count() const;

  bool is_binary() const;
  bool is_valid() const;

  preorder_iterator begin() const;
  preorder_iterator end() const;

private:
  void convert_tree(corax_utree_t *corax_tree);

  std::shared_ptr<node_t> _tree;
};
} // namespace bigrig
