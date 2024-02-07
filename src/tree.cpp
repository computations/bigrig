#include "tree.hpp"

#include "logger.hpp"

#include <stdexcept>
#include <corax/core/common.h>

namespace bigrig {
tree_t::tree_t(const std::filesystem::path &tree_filename) {
  auto corax_tree = corax_utree_parse_newick_rooted(tree_filename.c_str());

  convert_tree(corax_tree);
}

tree_t::tree_t(const std::string &tree_str) {
  auto corax_tree = corax_utree_parse_newick_string_rooted(tree_str.c_str());

  convert_tree(corax_tree);
}

/**
 * Get a dist by "string_id" key. The string_id is either the label, if the node
 * has one, or a string version of the assigned id.
 */
std::optional<dist_t>
tree_t::get_dist_by_string_id(const std::string &key) const {
  for (const auto &n : *this) {
    if (n->string_id() == key) { return n->final_state(); }
  }
  return {};
}

std::string tree_t::to_newick() const {
  std::stringstream oss;
  _tree->to_newick(oss);
  return oss.str();
}

std::string tree_t::to_newick(
    std::function<void(std::ostream &, const node_t &)> cb) const {
  std::stringstream oss;
  _tree->to_newick(oss, cb);
  return oss.str();
}

std::string tree_t::to_phylip_body() const {
  std::stringstream oss;
  to_phylip_body(oss);
  return oss.str();
}

std::string tree_t::to_phylip_body_extended() const {
  std::stringstream oss;
  to_phylip_body(oss, true);
  return oss.str();
}

std::ostream &tree_t::to_phylip_body(std::ostream &os, bool all) const {
  size_t padding = _tree->get_string_id_len_max(all) + 1;
  for (const auto &c : *this) {
    c->to_phylip_line(os, padding, all);
    os << "\n";
  }
  // remove the last newline
  os.seekp(-1, std::ios_base::end);
  return os;
}

size_t tree_t::node_count() const { return _tree->node_count(); }
size_t tree_t::leaf_count() const { return _tree->leaf_count(); }

bool tree_t::is_binary() const { return _tree->is_binary(); }
bool tree_t::is_valid() const {
  if (!_tree) { return false; }
  return true;
}

preorder_iterator tree_t::begin() const { return preorder_iterator(_tree); }
preorder_iterator tree_t::end() const { return preorder_iterator(); }

/**
 * Converts a coraxlib tree into a local tree.
 *
 * Importantly, we don't try to resolve polytomies in this function. In reality,
 * once a tree is converted, we still need to check that it is binary.
 */
void tree_t::convert_tree(corax_utree_t *corax_tree) {
  if (corax_tree == nullptr) {
    LOG_ERROR("We failed to parse the tree: %s", corax_errmsg);
    return;
  }
  _tree = std::make_unique<node_t>();
  _tree->add_child(std::make_shared<node_t>(corax_tree->vroot->back));
  _tree->add_child(std::make_shared<node_t>(corax_tree->vroot->next->back));
  if (corax_tree->vroot->label) { _tree->set_label(corax_tree->vroot->label); }
  _tree->assign_id_root();
  _tree->assign_abs_time_root();

  corax_utree_destroy(corax_tree, nullptr);
}
} // namespace bigrig
