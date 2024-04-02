#include "node.hpp"

#include <cstdlib>

namespace bigrig {

/**
 * Constructs a tree from a corax tree node.
 *
 * This function will copy the label and branch length, but nothing else.
 */
node_t::node_t(corax_unode_t *n) {
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

void node_t::add_child(const std::shared_ptr<node_t> &n) {
  _children.push_back(n);
}

/**
 * Convert the node_t into a newick string, with a formatting callback
 *
 * The callback is to format the label and branch length parameters, and any
 * other information that needs to be included. For example, the callback can
 * also construct NHX extension information.
 */
std::ostream &node_t::to_newick(
    std::ostream                                       &os,
    std::function<void(std::ostream &, const node_t &)> cb) const {
  if (!_children.empty()) { os << "("; }

  for (size_t i = 0; i < _children.size(); ++i) {
    const auto &c = _children[i];
    c->to_newick(os, cb);
    if (i != _children.size() - 1) { os << ","; }
  }
  if (!_children.empty()) { os << ")"; }

  cb(os, *this);

  return os;
}

/**
 * Convert the node_t to newick, standard method
 *
 * This will convert the node into a newick string and store it in `os`. This
 * function is just a wrapper around the "normal" callback, which is canned and
 * ready.
 */
std::ostream &node_t::to_newick(std::ostream &os) const {
  constexpr auto cb = [](std::ostream &os, const node_t &n) {
    os << n.string_id() << ":" << n._brlen;
  };
  return to_newick(os, cb);
}

/**
 * Converts the current node into a phylip row.
 *
 * When outputting results, the dist and taxa are stored as a row in a phylip
 * file. Here, we convert the current node into a row. We also output the
 * simulated dists for internal nodes. If this behavior is required, then all
 * should be set to true.
 */
std::ostream &
node_t::to_phylip_line(std::ostream &os, size_t pad_to, bool all) const {
  if (_children.size() == 0 || all) {
    auto tmp_name = string_id();
    os << tmp_name;
    if (pad_to != 0) {
      if (pad_to < tmp_name.size()) {
        throw std::runtime_error{"invalid padding, will overflow"};
      }
      pad_to -= tmp_name.size();
    }

    for (size_t i = 0; i < pad_to; ++i) { os << " "; }

    os << _final_state;
    os << "\n";
  }
  return os;
}

/**
 * Count the number of leaves in the tree.
 *
 * Linear complexity, use sparingly.
 */
size_t node_t::leaf_count() const {
  if (is_leaf()) { return 1; }
  size_t count = 0;
  for (const auto &child : _children) { count += child->leaf_count(); }
  return count;
}

/**
 * Count the number of nodes in the tree.
 *
 * Linear complexity, use sparingly.
 */
size_t node_t::node_count() const {
  size_t count = 0;
  for (const auto &child : _children) { count += child->node_count(); }
  return count + 1;
}

/**
 * Start assigning ids. Starts from index 0.
 */
void node_t::assign_id_root() { assign_id(0); }

/**
 * Recursively assign ids in a preorder fashion.
 */
size_t node_t::assign_id(size_t next) {
  if (is_leaf()) { return next; }
  _node_id = next++;
  for (const auto &c : _children) { next = c->assign_id(next); }
  return next;
}

size_t node_t::get_string_id_len_max(bool all) {
  return get_string_id_len_max(0, all);
}

/**
 * Get the length of the label as a string.
 *
 * Used to compute the padding for the phylip file
 */
size_t node_t::get_string_id_len_max(size_t max, bool all) {
  if (is_leaf() || all) {
    size_t label_size = string_id().size();
    max               = std::max(label_size, max);
  }
  for (const auto &c : _children) { max = c->get_string_id_len_max(max, all); }
  return max;
}

std::string node_t::label() const { return _label; }
double      node_t::brlen() const { return _brlen; }
double      node_t::abs_time() const { return _abs_time; }
double      node_t::abs_time_at_start() const { return _abs_time - brlen(); }
size_t      node_t::node_id() const { return _node_id; }
dist_t      node_t::final_state() const { return _final_state; }
std::string node_t::string_id() const {
  return is_leaf() ? _label : std::to_string(_node_id);
}

split_t node_t::node_split() const { return _split; }

std::vector<std::shared_ptr<node_t>> &node_t::children() { return _children; }

std::vector<std::shared_ptr<node_t>> node_t::children() const {
  return _children;
}

std::vector<transition_t>  node_t::transitions() const { return _transitions; }
std::vector<transition_t> &node_t::transitions() { return _transitions; }

/**
 * Compute the absolute time for the current node, as measured from the root.
 */
void node_t::assign_abs_time(double t) {
  _abs_time = t + _brlen;
  for (auto &c : children()) { c->assign_abs_time(_abs_time); }
}

void node_t::assign_abs_time_root() { assign_abs_time(0); }

void node_t::set_label(const std::string &str) { _label = str; }

bool node_t::is_binary() const {
  if (_children.size() != 0 && _children.size() != 2) { return false; }
  for (const auto &c : _children) {
    if (!c->is_binary()) { return false; }
  }
  return true;
}

bool node_t::is_valid() const {
  if (!validate_periods()) {
    LOG_ERROR("Failed to validate periods for node '%s'", string_id().c_str());
    return false;
  }

  for (const auto &c : _children) {
    if (!c->is_valid()) { return false; }
  }

  return true;
}

bool node_t::validate_periods() const {
  if (_periods.empty()) {
    LOG_ERROR("Period vector for node '%s' is empty", string_id().c_str());
    return false;
  }

  constexpr double abstol       = 1e-9;
  double           total_length = 0.0;
  for (auto p : _periods) { total_length += p.length(); }
  if (std::abs(total_length - _brlen) > abstol) {
    LOG_ERROR("Total period length for node '%s' is incorrect",
              string_id().c_str());
    return false;
  }

  return true;
}

void node_t::assign_periods(const std::vector<period_t> &periods) {
  parse_periods(periods);
  for (const auto &c : _children) { c->assign_periods(periods); }
}

void node_t::parse_periods(const std::vector<period_t> &periods) {
  _periods.clear();

  /* find the starting period */
  auto start_period_itr = periods.begin();
  while (start_period_itr != periods.end()
         && start_period_itr->end() < abs_time_at_start()) {
    start_period_itr++;
  }

  /* truncte the beginning of the period */
  auto  &start_period      = *start_period_itr;
  double start_time_offset = abs_time_at_start() - start_period.start();

  double new_length
      = std::min(start_period.length() - start_time_offset, brlen());

  _periods.emplace_back(
      abs_time_at_start(), new_length, start_period.model_ptr());

  /* find the last period */
  auto end_period_itr = start_period_itr;
  while (end_period_itr != periods.end()
         && end_period_itr->end() < abs_time()) {
    end_period_itr++;
  }

  /* insert up to the last period */
  for (start_period_itr++; start_period_itr != end_period_itr;
       start_period_itr++) {
    _periods.emplace_back(*start_period_itr);
  }

  /* truncate and insert the last period, if there is another period */
  auto  &final_period = *end_period_itr;
  double final_length
      = final_period.length() - (abs_time() - final_period.end());
  _periods.emplace_back(
      final_period.start(), final_length, final_period.model_ptr());
}

} // namespace bigrig
