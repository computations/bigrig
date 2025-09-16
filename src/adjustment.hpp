#pragma once

#include "errors.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <filesystem>
#include <logger.hpp>
#include <optional>
#include <random>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace bigrig {

class adjustment_matrix_t;

enum class adjustment_matrix_symmetry { symmetric, nonsymmetric, unknown };

struct adjacency_arc_t {
  std::string from;
  std::string to;
  double      value;

  [[nodiscard]] bool reverse(const adjacency_arc_t &a) const {
    return from == a.to && to == a.from;
  }
};

struct adjacency_graph_t {
  std::vector<adjacency_arc_t> adjacencies;
  adjustment_matrix_symmetry   type;

  size_t size() const { return adjacencies.size(); }
};

struct adjustment_matrix_params_t {
  std::expected<std::filesystem::path, io_err> matrix_filename;
  std::optional<adjacency_graph_t>             adjustments;
  std::optional<double>                        exponent;
  std::optional<bool>                          simulate;
};

typedef std::vector<double> region_adjustment_map_t;

class adjustment_matrix_t {
public:
  adjustment_matrix_t() = default;

  adjustment_matrix_t(const adjustment_matrix_params_t       &params,
                      const std::vector<std::string>         &area_names,
                      std::uniform_random_bit_generator auto &gen)
      : _region_count{area_names.size()} {
    if (params.adjustments && params.simulate.value_or(false)) {
      MESSAGE_ERROR(
          "Both an adjustment graph file {} and the simulate flag were given");
      throw std::runtime_error{"failed to setup adjustment matrix"};
    }

    if (auto matrix = params.adjustments; matrix) {
      _type = matrix->type;
      if (matrix_size() != matrix->size()) {
        throw std::runtime_error{"The matrix is not the correct size"};
      }
      _map.resize(matrix_size());
      auto area_name_inverse_map
          = std::views::enumerate(area_names)
          | std::views::transform(
                [](const auto &i) -> std::pair<std::string, size_t> {
                  auto [a, b] = i;
                  return {b, a};
                })
          | std::ranges::to<std::unordered_map>();

      for (auto &arc : matrix->adjacencies) {
        auto from_idx = area_name_inverse_map.at(arc.from);
        auto to_idx   = area_name_inverse_map.at(arc.to);

        _map[get_index(from_idx, to_idx)] = arc.value;
        if (is_symmetric()) { _map[get_index(to_idx, from_idx)] = arc.value; }
      }
    }

    if (params.simulate.value_or(false)) { simulate(gen); }

    if (auto expo = params.exponent; expo) { apply_exponent(*expo); }
  }

  adjustment_matrix_t(const adjustment_matrix_t &rhs)         = default;
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
    _region_count = region_count;
    simulate(gen);
  }

  void simulate(std::uniform_random_bit_generator auto &gen) {
    constexpr double alpha_default = 2.0;
    constexpr double beta_default  = 2.0;
    simulate(alpha_default, beta_default, gen);
  }


  void simulate(double                                  alpha,
                double                                  beta,
                std::uniform_random_bit_generator auto &gen) {
    _map.resize(_region_count * _region_count);

    std::gamma_distribution<> dis(alpha, beta);

    /* Make a symmetric matrix */
    for (size_t i = 0; i < _region_count; ++i) {
      for (size_t j = i + 1; j < _region_count; ++j) {
        _map[i * _region_count + j] = _map[j * _region_count + i]
            = dis(gen) + 1;
      }
    }

    apply_exponent(-1.0);
  }

  size_t matrix_size() const {
    if (_type == adjustment_matrix_symmetry::symmetric) {
      return (_region_count * (_region_count - 1)) / 2;
    }
    return _region_count * _region_count;
  }

  bool is_symmetric() const {
    return _type == adjustment_matrix_symmetry::symmetric;
  }

  size_t get_row_size() const { return _region_count; }

private:
  size_t get_index(size_t from, size_t to) const {
    if (from == to) { return 0.0; }
    if (is_symmetric() && from > to) { std::swap(from, to); }
    size_t triangle_adjustment = (_type == adjustment_matrix_symmetry::symmetric
                                      ? (from + 1) * (from + 2) / 2
                                      : 0);
    return from * _region_count + to - triangle_adjustment;
  }

  /**
   * We are going to use this function to setup a n(n-1)/2 vector representing
   * the upper triangular form.
   */
  void set_symmetric_map(const adjacency_graph_t &) {
    throw std::runtime_error{"Not implemented"};
  }

  void set_unsymmetric_map(const adjacency_graph_t &) {
    throw std::runtime_error{"Not implemented"};
  }

  region_adjustment_map_t    _map;
  adjustment_matrix_symmetry _type;
  size_t                     _region_count;
};

} // namespace bigrig
