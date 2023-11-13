#pragma once

#include "corax/tree/utree.h"
#include "dist.hpp"
#include "model.hpp"

#include <functional>
#include <logger.hpp>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace biogeosim {

class node_t {
public:
  node_t() = default;

  node_t(corax_unode_t *n) {
    _brlen = n->length;
    if (n->label) { _label = n->label; }
    if (n->next == nullptr) { return; }

    auto start = n;
    n          = n->next;
    while (n != start) {
      _children.push_back(std::make_shared<node_t>(n->back));
      n = n->next;
    }
  }

  void add_child(const std::shared_ptr<node_t> &n) { _children.push_back(n); }

  void sample(dist_t                                  initial_distribution,
              const substitution_model_t             &model,
              std::uniform_random_bit_generator auto &gen) {
    LOG_DEBUG("Node sampling with initial_distribution = %b",
              initial_distribution);
    _transitions = generate_samples(initial_distribution, _brlen, model, gen);
    if (_transitions.size() == 0) {
      _final_state = initial_distribution;
    } else {
      _final_state = _transitions.back().final_state;
    }
    _split = split_dist(_final_state, model, gen);

    if (!is_leaf()) {
      _children[0]->sample(_split.left, model, gen);
      _children[1]->sample(_split.right, model, gen);
    }
  }

  std::ostream &
  to_newick(std::ostream                                       &os,
            std::function<void(std::ostream &, const node_t &)> cb) const {
    if (!_children.empty()) { os << "("; }

    for (size_t i = 0; i < _children.size(); ++i) {
      const auto &c = _children[i];
      c->to_newick(os, cb);
      if (i != _children.size() - 1) { os << ", "; }
    }
    if (!_children.empty()) { os << ")"; }

    cb(os, *this);

    return os;
  }

  std::ostream &to_newick(std::ostream &os) const {
    auto cb = [](std::ostream &os, const node_t &n) {
      os << n.string_id() << ":" << n._brlen;
    };
    return to_newick(os, cb);
  }

  std::ostream &to_phylip_line(std::ostream &os, bool all = false) const {
    for (const auto child : _children) { child->to_phylip_line(os, all); }
    if (_children.size() == 0 || all) {
      os << string_id() << " " << _final_state << "\n";
    }
    return os;
  }

  inline bool is_leaf() const { return _children.size() == 0; }

  size_t leaf_count() const {
    if (is_leaf()) { return 1; }
    size_t count = 0;
    for (const auto child : _children) { count += child->leaf_count(); }
    return count;
  }

  size_t node_count() const {
    size_t count = 0;
    for (const auto child : _children) { count += child->node_count(); }
    return count + 1;
  }

  void assign_id_root() { assign_id(0); }

  size_t assign_id(size_t next) {
    if (is_leaf()) { return next; }
    _node_id = next++;
    for (const auto &c : _children) { next = c->assign_id(next); }
    return next;
  }

  std::string label() const { return _label; }
  double      brlen() const { return _brlen; }
  double      abs_time() const { return _abs_time; }
  size_t      node_id() const { return _node_id; }
  dist_t      final_state() const { return _final_state; }
  std::string string_id() const {
    return is_leaf() ? _label : std::to_string(_node_id);
  }

  split_t node_split() const { return _split; }

  std::vector<std::shared_ptr<node_t>> &children() { return _children; }

  std::vector<std::shared_ptr<node_t>> children() const { return _children; }

  std::vector<transition_t>  transitions() const { return _transitions; }
  std::vector<transition_t> &transitions() { return _transitions; }

  void assign_abs_time(double t) {
    _abs_time = t + _brlen;
    for (auto &c : children()) { c->assign_abs_time(_abs_time); }
  }

  void assign_abs_time_root() { assign_abs_time(0); }

private:
  double                               _brlen;
  double                               _abs_time;
  dist_t                               _final_state;
  split_t                              _split;
  std::string                          _label;
  std::vector<std::shared_ptr<node_t>> _children;
  std::vector<transition_t>            _transitions;
  size_t                               _node_id;
};
} // namespace biogeosim
