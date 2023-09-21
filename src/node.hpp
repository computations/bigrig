#pragma once

#include "corax/tree/utree.h"
#include "dist.hpp"
#include "model.hpp"
#include <logger.hpp>
#include <memory>
#include <random>
#include <sstream>
#include <vector>

namespace biogeosim {

class node_t {
public:
  node_t() = default;

  node_t(corax_unode_t *n) {
    _brlen = n->length;
    if (n->label) {
      _label = n->label;
    }
    if (n->next == nullptr) {
      return;
    }

    auto start = n;
    n = n->next;
    while (n != start) {
      _children.push_back(std::make_shared<node_t>(n->back));
      n = n->next;
    }
  }

  void add_child(const std::shared_ptr<node_t> &n) { _children.push_back(n); }

  void sample(dist_t initial_distribution, const substitution_model_t &model,
              std::uniform_random_bit_generator auto &gen) {
    _transitions = generate_samples(_brlen, model, gen);
    auto final_state = _transitions.back().final_state;
    for (auto &child : _children) {
      child->sample(final_state, model, gen);
    }
  }

  std::ostream &to_newick(std::ostream &os) const {
    if (!_children.empty()) {
      os << "(";
    }

    for (size_t i = 0; i < _children.size(); ++i) {
      const auto &c = _children[i];
      c->to_newick(os);
      if (i != _children.size() - 1) {
        os << ", ";
      }
    }

    if (!_children.empty()) {
      os << ")";
    }

    os << _label << ":" << _brlen;

    return os;
  }

private:
  double _brlen;
  std::string _label;
  std::vector<std::shared_ptr<node_t>> _children;
  std::vector<transition_t> _transitions;
};
} // namespace biogeosim
