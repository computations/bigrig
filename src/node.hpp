#pragma once

#include "corax/tree/utree.h"
#include "dist.hpp"
#include "model.hpp"
#include <logger.hpp>
#include <memory>
#include <random>
#include <vector>

namespace biogeosim {

class node_t {
public:
  node_t() = default;

  node_t(corax_unode_t *n) {
    _brlen = n->length;
    if (n->next == nullptr) {
      return;
    }
    auto start = n;
    n = n->next;
    while (n != start) {
      _children.push_back(std::make_shared<node_t>(n));
    }
  }

  void add_child(const std::shared_ptr<node_t> &n) { _children.push_back(n); }

  void sample(dist_t initial_distribution, const substitution_model_t &model,
              std::uniform_random_bit_generator auto &gen) {}

private:
  std::vector<std::shared_ptr<node_t>> _children;
  double _brlen;
  std::vector<transition_t> _transitions;
};
} // namespace biogeosim
