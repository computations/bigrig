#include "node.hpp"

namespace bigrig {

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

std::ostream &node_t::to_newick(std::ostream &os) const {
  auto cb = [](std::ostream &os, const node_t &n) {
    os << n.string_id() << ":" << n._brlen;
  };
  return to_newick(os, cb);
}

std::ostream &
node_t::to_phylip_line(std::ostream &os, size_t pad_to, bool all) const {
  for (const auto &child : _children) {
    child->to_phylip_line(os, pad_to, all);
  }
  if (_children.size() == 0 || all) {
    // os << string_id() << " " << _final_state << "\n";
    auto tmp_name = string_id();
    os << tmp_name;
    if (pad_to != 0) {
      if (pad_to < tmp_name.size()) {
        throw std::runtime_error{"invalid padding, will overflow"};
      }
      pad_to -= tmp_name.size();
    }

    for (size_t i = 0; i < pad_to; ++i) { os << " "; }

    os << _final_state << "\n";
  }
  return os;
}

size_t node_t::leaf_count() const {
  if (is_leaf()) { return 1; }
  size_t count = 0;
  for (const auto &child : _children) { count += child->leaf_count(); }
  return count;
}

size_t node_t::node_count() const {
  size_t count = 0;
  for (const auto &child : _children) { count += child->node_count(); }
  return count + 1;
}

void node_t::assign_id_root() { assign_id(0); }

size_t node_t::assign_id(size_t next) {
  if (is_leaf()) { return next; }
  _node_id = next++;
  for (const auto &c : _children) { next = c->assign_id(next); }
  return next;
}

size_t node_t::get_string_id_len_max(bool all) {
  return get_string_id_len_max(0, all);
}

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

void node_t::assign_abs_time(double t) {
  _abs_time = t + _brlen;
  for (auto &c : children()) { c->assign_abs_time(_abs_time); }
}

void node_t::assign_abs_time_root() { assign_abs_time(0); }

void node_t::set_label(const std::string &str) { _label = str; }
} // namespace bigrig
