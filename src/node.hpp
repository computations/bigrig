#pragma once

#include "dist.hpp"
#include "model.hpp"
#include "period.hpp"
#include "split.hpp"
#include "spread.hpp"

#include <cmath>
#include <corax/tree/utree.h>
#include <functional>
#include <logger.hpp>
#include <random>
#include <string>
#include <vector>

namespace bigrig {

class node_t {
public:
  node_t()               = default;
  node_t(const node_t &) = default;
  node_t(node_t &&)      = default;

  node_t(corax_unode_t *n);

  void add_child(const std::shared_ptr<node_t> &n);

  std::shared_ptr<node_t> clone() const;

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
    _split.period_index = _periods.back().index();

    if (!is_leaf()) {
      _children[0]->simulate(_split.left, gen, mode);
      _children[1]->simulate(_split.right, gen, mode);
    }
  }

  void simulate_tree(dist_t               initial_distribution,
                     double               time_left,
                     const period_list_t &periods,
                     std::uniform_random_bit_generator auto &gen,
                     operation_mode_e mode = operation_mode_e::FAST) {
    auto   dist     = initial_distribution;
    double leftover = 0.0;

    while (true) {
      /*
       * two cases:
       *   - Dispersion/extinction event
       *   - Specitation event
       *
       * The best way to do this is to roll a time, and then which of the two
       * events it is, and then handle that case individually.
       */
      auto  period = periods.get(_abs_time + _brlen);
      auto &model  = period.model();

      auto total_rate = model.total_event_weight(dist);

      if (total_rate <= 0.0) {
        MESSAGE_ERROR("Rate while simulating the tree is invalid");
        throw std::runtime_error{"total_rate is not positive"};
      }

      if (!std::isfinite(_brlen)) {
        MESSAGE_ERROR("Simulation ran to infinity");
        throw std::runtime_error{"brlen is infinite"};
      }

      std::exponential_distribution<double> waiting_time_die(total_rate);

      double waiting_time = waiting_time_die(gen) + leftover;
      leftover            = 0.0;

      if (waiting_time + _brlen > time_left) {
        _brlen       = time_left;
        _final_state = dist;
        return;
      }

      /* check if the period is over */
      double period_time_left = period.length() - (_abs_time + _brlen);
      if (period_time_left < waiting_time) {
        leftover  = period_time_left;
        _brlen   += period_time_left;
        continue;
      }

      _brlen += waiting_time;

      double speciation_rate = model.total_speciation_weight(dist);

      std::bernoulli_distribution type_coin(speciation_rate / total_rate);
      if (type_coin(gen)) {
        LOG_DEBUG("Rolled a cladogenesis event. Time left %f", time_left);
        auto res   = split_dist(dist, model, gen);
        _split     = res;
        _abs_time += _brlen;

        _children.emplace_back(new node_t);
        _children.back()->_abs_time = _abs_time;
        _children.back()->simulate_tree(
            res.left, time_left - _brlen, periods, gen, mode);

        _children.emplace_back(new node_t);
        _children.back()->_abs_time = _abs_time;
        _children.back()->simulate_tree(
            res.right, time_left - _brlen, periods, gen, mode);

        _final_state = dist;
        return;
      } else {
        LOG_DEBUG("Rolled a transition event. Time left %f", time_left);
        auto res         = spread_flip_region(dist, model, gen);
        res.period_index = period.index();
        res.waiting_time = waiting_time;
        _transitions.push_back(res);

        if (res.final_state.empty()) {
          _final_state = _transitions.back().final_state;
          return;
        }
      }
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

  size_t reconstructed_leaf_count() const;
  size_t reconstructed_leaf_count(double height) const;

  double reconstructed_brlen_sum() const;
  double reconstructed_brlen_sum(double height) const;

  bool is_binary() const;
  bool is_valid() const;
  bool validate_periods() const;

  void assign_periods(const period_list_t &periods);
  void assign_id_root();

  size_t assign_id(size_t next);

  size_t get_string_id_len_max(bool all);

  size_t get_string_id_len_max(size_t max, bool all);

  std::string label() const;
  void        assign_label(const std::string &);
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

  void   set_label(const std::string &str);
  dist_t start_range() const;

  double brlen_sum() const;
  double max_tree_height() const;

  bool contractible() const;

  bool prunable() const;
  void prune();

private:
  void     parse_periods(const std::vector<period_t> &periods);
  period_t clamp_period(const period_t &p) const;

  double                               _brlen    = 0;
  double                               _abs_time = 0;
  dist_t                               _final_state;
  split_t                              _split;
  std::string                          _label;
  std::vector<std::shared_ptr<node_t>> _children;
  std::vector<transition_t>            _transitions;
  period_list_t                        _periods;
  size_t                               _node_id;
};
} // namespace bigrig
