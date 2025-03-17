#pragma once
#include "dist.hpp"
#include "iterator.hpp"
#include "model.hpp"
#include "node.hpp"
#include "util.hpp"

#include <corax/corax.hpp>
#include <functional>
#include <optional>
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
  explicit tree_t() : _tree{new node_t} {};

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
                std::uniform_random_bit_generator auto &gen) {
    LOG_DEBUG("Starting sample with init dist = %lb",
              static_cast<uint64_t>(initial_distribution));
    _tree->simulate(initial_distribution, gen, _mode);
  }

  void simulate_tree(dist_t               initial_distribution,
                     const period_list_t &periods,
                     double               tree_height,
                     std::uniform_random_bit_generator auto &gen) {
    _tree->simulate_tree(
        initial_distribution, tree_height, periods, gen, _mode);

    _tree->assign_id_root();
    set_tree_labels();
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
  bool is_ready(bool) const;

  size_t region_count() const;

  preorder_iterator begin() const;
  preorder_iterator end() const;

  void set_mode(operation_mode_e mode);

  void set_periods(const period_list_t &periods);
  void set_periods(const period_t &periods);

  dist_t get_root_range() const;

  void set_tree_labels() {
    size_t leaf_itr = 0;
    for (auto n : *this) {
      if (n->is_leaf()) {
        n->assign_label(util::compute_base26(leaf_itr));
        leaf_itr++;
      }
    }
  }

private:
  void convert_tree(corax_utree_t *corax_tree);

  std::shared_ptr<node_t> _tree;
  operation_mode_e        _mode = operation_mode_e::FAST;
};
} // namespace bigrig
