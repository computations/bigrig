#pragma once

#include "dist.hpp"
#include "model.hpp"
#include "period.hpp"
#include "split.hpp"

#include <corax/tree/utree.h>
#include <functional>
#include <iterator>
#include <logger.hpp>
#include <string>
#include <vector>

namespace bigrig {

class node_t {
public:
  node_t() = default;

  node_t(corax_unode_t *n);

  void add_child(const std::shared_ptr<node_t> &n);

  /**
   * Run the simulation, given the initial distribution, model and RNG.
   *
   * Modifies the state of the node_t, so that once the function is finished,
   * the _transitions vector is full, _final_state is occupied, and _split is
   * computed.
   *
   * Also, this is a recursive function that is top down. After the results for
   * this node are computed, the results for the children are computed.
   */
  void simulate(dist_t                                  initial_distribution,
                std::uniform_random_bit_generator auto &gen,
                operation_mode_e mode = operation_mode_e::FAST) {
    LOG_DEBUG("Node sampling with initial_distribution = %s",
              initial_distribution.to_str().c_str());
    _transitions
        = simulate_transitions(initial_distribution, _periods, gen, mode);
    LOG_DEBUG("Finished sampling with %lu transitions", _transitions.size());
    if (_transitions.size() == 0) {
      _final_state = initial_distribution;
    } else {
      _final_state = _transitions.back().final_state;
    }
    _split = split_dist(_final_state, _periods.back().model(), gen, mode);

    if (!is_leaf()) {
      _children[0]->simulate(_split.left, gen, mode);
      _children[1]->simulate(_split.right, gen, mode);
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

  bool is_binary() const;
  bool is_valid() const;
  bool validate_periods() const;

  void assign_periods(const std::vector<period_t> &periods);
  void assign_id_root();

  size_t assign_id(size_t next);

  size_t get_string_id_len_max(bool all);

  size_t get_string_id_len_max(size_t max, bool all);

  std::string label() const;
  double      brlen() const;
  double      abs_time() const;
  double      abs_time_at_start() const;
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
  void     parse_periods(const std::vector<period_t> &periods);
  period_t clamp_period(const period_t &p) const;

  double                               _brlen;
  double                               _abs_time;
  dist_t                               _final_state;
  split_t                              _split;
  std::string                          _label;
  std::vector<std::shared_ptr<node_t>> _children;
  std::vector<transition_t>            _transitions;
  std::vector<period_t>                _periods;
  size_t                               _node_id;
};
} // namespace bigrig
