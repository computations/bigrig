#include "split.hpp"

namespace biogeosim {
std::string split_t::to_nhx_string() const {
  std::ostringstream oss;
  oss << "left-split=" << left << ":"
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
} // namespace biogeosim
