#include "split.hpp"

namespace bigrig {
std::string split_t::to_nhx_string() const {
  std::ostringstream oss;
  oss << "init-dist=" << top << ":"
      << "left-split=" << left << ":"
      << "right-split=" << right << ":";
  oss << "split-type=";

  return oss.str();
}

std::string type_string(const split_type_e &st) {
  switch (st) {
  case split_type_e::singleton:
    return "singleton";
  case split_type_e::allopatric:
    return "allopatric";
  case split_type_e::sympatric:
    return "sympatric";
  case split_type_e::jump:
    return "jump";
  case split_type_e::invalid:
    return "invalid";
  }
  throw std::runtime_error{"Did not cover all cases"};
}

std::string split_t::to_type_string() const { return type_string(type); }

/**
 * Given a triplet of `dist_t`s determine what type of split it is between:
 *  - singleton (i.e. copy)
 *  - sympatric
 *  - allopatric
 *  - invalid
 * Generally speaking, this function is pretty slow, and so should be avoided.
 * If possible
 */
split_type_e
determine_split_type(dist_t init_dist, dist_t left_dist, dist_t right_dist) {
  if (left_dist.full_region_count() > right_dist.full_region_count()) {
    std::swap(left_dist, right_dist);
  }

  if ((left_dist | right_dist) == init_dist) {
    if (left_dist == right_dist && left_dist.full_region_count() == 1) {
      return split_type_e::singleton;
    }
    if (!left_dist.singleton() && !right_dist.singleton()) {
      return split_type_e::invalid;
    }
    if ((left_dist & right_dist).full_region_count() == 1) {
      return split_type_e::sympatric;
    }
    if ((left_dist & right_dist).full_region_count() == 0) {
      return split_type_e::allopatric;
    }
  } else if (left_dist.singleton()
             && ((left_dist & init_dist).full_region_count() == 0)
             && (right_dist == init_dist)) {
    return split_type_e::jump;
  }
  return split_type_e::invalid;
}
} // namespace bigrig
