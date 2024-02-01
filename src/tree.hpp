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

namespace biogeosim {
class tree_t {
public:
  explicit tree_t(const std::filesystem::path &tree_filename) {
    auto corax_tree = corax_utree_parse_newick_rooted(tree_filename.c_str());

    convert_tree(corax_tree);
  }

  explicit tree_t(const std::string &tree_str) {
    auto corax_tree = corax_utree_parse_newick_string_rooted(tree_str.c_str());

    convert_tree(corax_tree);
  }

  tree_t(const tree_t &)            = delete;
  tree_t &operator=(const tree_t &) = delete;

  tree_t(tree_t &&)            = default;
  tree_t &operator=(tree_t &&) = default;

  void sample(dist_t                                  initial_distribution,
              const substitution_model_t             &model,
              std::uniform_random_bit_generator auto &gen) {
    if (!valid_dist(initial_distribution, model)) {
      throw invalid_dist{"Invalid dist provided as a start dist"};
    }
    LOG_DEBUG("Starting sample with init dist = %lb",
              static_cast<uint64_t>(initial_distribution));
    _tree->sample(initial_distribution, model, gen);
  }

  std::optional<dist_t> get_dist_by_string_id(const std::string &key) const {
    for (const auto &n : *this) {
      if (n->string_id() == key) { return n->final_state(); }
    }
    return {};
  }

  std::string to_newick() const {
    std::stringstream oss;
    _tree->to_newick(oss);
    return oss.str();
  }

  std::string
  to_newick(std::function<void(std::ostream &, const node_t &)> cb) const {
    std::stringstream oss;
    _tree->to_newick(oss, cb);
    return oss.str();
  }

  std::string to_phylip_body() const {
    std::stringstream oss;
    to_phylip_body(oss);
    return oss.str();
  }

  std::string to_phylip_body_extended() const {
    std::stringstream oss;
    to_phylip_body(oss, true);
    return oss.str();
  }

  std::ostream &to_phylip_body(std::ostream &os, bool all = false) const {
    size_t padding = _tree->get_string_id_len_max(all) + 1;
    _tree->to_phylip_line(os, padding, all);
    return os;
  }

  size_t node_count() const { return _tree->node_count(); }

  size_t leaf_count() const { return _tree->leaf_count(); }

  preorder_iterator begin() const { return preorder_iterator(_tree); }
  preorder_iterator end() const { return preorder_iterator(); }

private:
  void convert_tree(corax_utree_t *corax_tree) {
    _tree = std::make_unique<node_t>();
    _tree->add_child(std::make_shared<node_t>(corax_tree->vroot->back));
    _tree->add_child(std::make_shared<node_t>(corax_tree->vroot->next->back));
    if (corax_tree->vroot->label) {
      _tree->set_label(corax_tree->vroot->label);
    }
    _tree->assign_id_root();
    _tree->assign_abs_time_root();

    corax_utree_destroy(corax_tree, nullptr);
  }

  std::shared_ptr<node_t> _tree;
};
} // namespace biogeosim
