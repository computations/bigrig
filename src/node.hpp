#pragma once

#include "corax/tree/utree.h"
#include "dist.hpp"
#include "model.hpp"
#include "split.hpp"

#include <functional>
#include <logger.hpp>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace bigrig {

class node_t {
public:
  node_t() = default;

  node_t(corax_unode_t *n);

  void add_child(const std::shared_ptr<node_t> &n);

  void sample(dist_t                                  initial_distribution,
              const biogeo_model_t                   &model,
              std::uniform_random_bit_generator auto &gen) {
    LOG_DEBUG("Node sampling with initial_distribution = %s",
              initial_distribution.to_str().c_str());
    _transitions = generate_samples(initial_distribution, _brlen, model, gen);
    LOG_DEBUG("Finished sampling with %lu transitions", _transitions.size());
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
            std::function<void(std::ostream &, const node_t &)> cb) const;

  std::ostream &to_newick(std::ostream &os) const;

  std::ostream &
  to_phylip_line(std::ostream &os, size_t pad_to = 0, bool all = false) const;

  inline bool is_leaf() const { return _children.size() == 0; }

  size_t leaf_count() const;

  size_t node_count() const;

  void assign_id_root();

  size_t assign_id(size_t next);

  size_t get_string_id_len_max(bool all);

  size_t get_string_id_len_max(size_t max, bool all);

  std::string label() const;
  double      brlen() const;
  double      abs_time() const;
  size_t      node_id() const;
  dist_t      final_state() const;
  std::string string_id() const;

  split_t node_split() const;

  std::vector<std::shared_ptr<node_t>> &children();

  std::vector<std::shared_ptr<node_t>> children() const;

  std::vector<transition_t>  transitions() const;
  std::vector<transition_t> &transitions();

  void assign_abs_time(double t);

  void assign_abs_time_root();

  void set_label(const std::string &str);

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
} // namespace bigrig
