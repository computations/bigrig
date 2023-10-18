#pragma once
#include "dist.hpp"
#include "exceptions.hpp"
#include "model.hpp"
#include "node.hpp"

#include <corax/corax.hpp>
#include <filesystem>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace biogeosim {
class tree_t {
public:
  tree_t(const std::filesystem::path &tree_filename) {
    auto corax_tree = corax_utree_parse_newick_rooted(tree_filename.c_str());

    _tree = std::make_unique<node_t>();
    _tree->add_child(std::make_shared<node_t>(corax_tree->vroot->back));
    _tree->add_child(std::make_shared<node_t>(corax_tree->vroot->next->back));
    _tree->assign_id_root();

    corax_utree_destroy(corax_tree, nullptr);
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
    LOG_DEBUG("Starting sample with init dist = %b", initial_distribution);
    _tree->sample(initial_distribution, model, gen);
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
  std::ostream &to_phylip_body(std::ostream &os) const {
    _tree->to_phylip_line(os);
    return os;
  }

  size_t node_count() const { return _tree->node_count(); }

  size_t leaf_count() const { return _tree->leaf_count(); }

private:
  std::unique_ptr<node_t> _tree;
};
} // namespace biogeosim
