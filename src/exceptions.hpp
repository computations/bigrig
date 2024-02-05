#pragma once

#include <stdexcept>
#include <string>

namespace bigrig {
class invalid_dist : public std::invalid_argument {
public:
  invalid_dist(const std::string &msg) : std::invalid_argument{msg} {}
};
} // namespace bigrig
