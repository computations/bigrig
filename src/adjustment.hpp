#pragma once

#include "errors.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <filesystem>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

namespace bigrig {

enum class adjustment_matrix_symmetry {
  symmetric,
  nonsymmetric,
  unknown
};

struct adjacency_arc_t {
  std::string from;
  std::string to;
  double      value;

  [[nodiscard]] bool reverse(const adjacency_arc_t &a) const {
    return from == a.to && to == a.from;
  }
};


struct adjustment_matrix_params_t {
  std::expected<std::filesystem::path, io_err> matrix_filename;
  std::optional<std::vector<adjacency_arc_t>>          adjustments;
  std::optional<double>                                exponent;
  std::optional<bool>                                  simulate;
};

typedef std::vector<double> region_adjustment_map_t;

class adjustment_matrix_t {
public:
  adjustment_matrix_t() = default;

  adjustment_matrix_t(const adjustment_matrix_params_t &params,
                      size_t                            region_count)
      : _region_count{region_count} {
    if (params.adjustments.has_value()) {
      auto &matrix = params.adjustments.value();
      if (matrix_size() != matrix.size()) {
        throw std::runtime_error{"The matrix is not the correct size"};
      }
      //_map = matrix;
    }

    if (params.exponent.has_value()) {
      apply_exponent(params.exponent.value());
    }
  }

  adjustment_matrix_t(const adjustment_matrix_t &rhs)         = default;
  //: _map{rhs._map}, _region_count{rhs._region_count} {}
  adjustment_matrix_t &operator=(const adjustment_matrix_t &) = default;

  double get_adjustment(size_t from, size_t to) const {
    return _map[from * _region_count + to];
  }

  void apply_exponent(double exponent) {
    std::for_each(_map.begin(), _map.end(), [exponent](double &adj) {
      if (adj != 0.0) { adj = std::pow(adj, exponent); }
    });
  }

  void simulate(size_t                                  region_count,
                std::uniform_random_bit_generator auto &gen) {
    simulate(region_count, 2.0, 2.0, gen);
  }

  void simulate(size_t                                  region_count,
                double                                  alpha,
                double                                  beta,
                std::uniform_random_bit_generator auto &gen) {
    _region_count = region_count;

    _map.resize(_region_count * _region_count);

    std::gamma_distribution<> dis(alpha, beta);

    /* Make a symmetric matrix */
    for (size_t i = 0; i < _region_count; ++i) {
      for (size_t j = i + 1; j < _region_count; ++j) {
        _map[i * _region_count + j] = _map[j * _region_count + i] = dis(gen);
      }
    }
  }

  size_t matrix_size() const { return _region_count * _region_count; }

private:
  region_adjustment_map_t _map;
  size_t                  _region_count;
};

} // namespace bigrig
