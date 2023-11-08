#pragma once
#include "node.hpp"

#include <stack>

namespace biogeosim {

class preorder_iterator {
public:
  using difference_type = std::ptrdiff_t;
  using value_type      = std::shared_ptr<node_t>;

  preorder_iterator() = default;
  explicit preorder_iterator(const std::shared_ptr<node_t> &n) {
    _stack.push(n);
  }
  value_type  operator*() const { return _stack.top(); }
  preorder_iterator &operator++() {
    auto tmp_node = _stack.top();
    _stack.pop();
    for (auto &c : tmp_node->children()) { _stack.push(c); }
    return *this;
  }
  preorder_iterator operator++(int) {
    auto tmp = *this;
    ++*this;
    return tmp;
  }

  bool operator==(const preorder_iterator &it) const {
    return (it._stack.empty() && _stack.empty())
        || (!_stack.empty() && it._stack.top() == _stack.top());
  }
  bool operator!=(const preorder_iterator &it) const { return !(it == *this); }

  std::shared_ptr<node_t> &node() { return _stack.top(); }

private:
  std::stack<std::shared_ptr<node_t>> _stack;
};

} // namespace biogeosim
