#pragma once

#include "adjustment.hpp"
#include "dist.hpp"

#include <cstddef>
#include <optional>
#include <sstream>

namespace bigrig {

struct rate_params_t {
  double dis;
  double ext;
};

/**
 * This data structure is a bit wacky. It really is just an array with named
 * elements (which the compiler is aware of). Sometimes I want to perform
 * operations on this as if it was an array. So, there are some WACKY functions
 * which basically cast the structure to an array, and iterate over the array to
 * compute a sum. Similarly, we use the same technique to normalize the
 * parameters.
 */
struct cladogenesis_params_t {
  using data_type = double;

  data_type allopatry;
  data_type sympatry;
  data_type copy;
  data_type jump;

  /**
   * Part of the wacky stuff that we do :). Computes the constant size of the
   * struct. Hypothetically, this should be baked in at compile time.
   */
  static constexpr size_t size() {
    static_assert(sizeof(cladogenesis_params_t) % sizeof(data_type) == 0,
                  "something is fucked");
    return sizeof(cladogenesis_params_t) / sizeof(data_type);
  }

  /**
   * Compute the sum of the rates. Again, uses the wacky trick of converting the
   * struct to a pointer and treating it as an array. Should work most of the
   * time :)
   */
  data_type sum() const {
    auto   start = as_ptr();
    double acc   = 0;
    for (size_t i = 0; i < size(); ++i) { acc += start[i]; }
    return acc;
  }

  /**
   * Produced a normalized version of the current parameters. Should be copy
   * elided, since everything is dumb data.
   */
  cladogenesis_params_t normalize() const {
    cladogenesis_params_t tmp{*this};

    tmp.normalize(sum());
    return tmp;
  }

  /**
   * Produce a helper debug string.
   */
  std::string to_debug_string() const {
    std::ostringstream oss;
    oss << "copy: " << copy << ", sympatry: " << sympatry
        << ", allopatry: " << allopatry << ", jump: " << jump;
    return oss.str();
  }

private:
  using array_type       = data_type *;
  using const_array_type = data_type const *;

  /**
   * Converts the current struct to an pointer, which can then be treated as an
   * array.
   */
  const_array_type as_ptr() const {
    return reinterpret_cast<const_array_type>(this);
  }

  array_type as_ptr() { return reinterpret_cast<array_type>(this); }

  /**
   * Actually normalizes the data, as a modification of the current instance.
   */
  void normalize(data_type d) {
    auto start = as_ptr();
    for (size_t i = 0; i < size(); ++i) { start[i] /= d; }
  }
};

/* Assert that the fields are in the right order */
static_assert(offsetof(cladogenesis_params_t, allopatry)
                  / sizeof(cladogenesis_params_t::data_type)
              == 0);
static_assert(offsetof(cladogenesis_params_t, sympatry)
                  / sizeof(cladogenesis_params_t::data_type)
              == 1);
static_assert(offsetof(cladogenesis_params_t, copy)
                  / sizeof(cladogenesis_params_t::data_type)
              == 2);
static_assert(offsetof(cladogenesis_params_t, jump)
                  / sizeof(cladogenesis_params_t::data_type)
              == 3);
static_assert(alignof(cladogenesis_params_t) == 8);
static_assert(sizeof(cladogenesis_params_t) == 32);

/* Check that jump is the last field */
static_assert(offsetof(cladogenesis_params_t, jump)
                  + sizeof(cladogenesis_params_t::data_type)
              == sizeof(cladogenesis_params_t));

struct tree_params_t {
  double cladogenesis;
};

struct per_region_params_t {
  std::variant<std::monostate, dist_t, std::string, size_t> region_id;
  std::optional<rate_params_t>                              rates;
  std::optional<cladogenesis_params_t>                      cladogenesis;
};

/**
 * Class containing the model parameters, which includes:
 * - Rate parameters,
 * - Cladogenesis parameters, and
 * - Region count.
 * Together, these make up the total amount of parameters for the model. In all,
 * this class mostly stores data, and doesn't do much other than that. In
 * particular, this class is intended to hold _static_ data. The assumption here
 * is that the data does not change during the runtime.
 */
class biogeo_model_t {
public:
  biogeo_model_t() = default;

  biogeo_model_t(double d, double e, bool duplicity)
      : _rate_params{.dis = d, .ext = e},
        _clad_params{
            .allopatry = 1.0, .sympatry = 1.0, .copy = 1.0, .jump = 0.0},
        _duplicity{duplicity} {}

  biogeo_model_t(const rate_params_t         &rp,
                 const cladogenesis_params_t &cp,
                 bool                         duplicity)
      : _rate_params{rp}, _clad_params{cp}, _duplicity{duplicity} {}

  /**
   * Returns the pair (e,d) for a simple 2 parameter dec model.
   */
  inline rate_params_t rates() const { return _rate_params; }

  biogeo_model_t &set_rate_params(rate_params_t p);

  biogeo_model_t &set_rate_params(double d, double e);

  biogeo_model_t &set_per_region_rate_params(size_t        region_index,
                                             rate_params_t p);

  biogeo_model_t &
  set_cladogenesis_params(double v, double s, double y, double j);

  biogeo_model_t &set_cladogenesis_params(const cladogenesis_params_t &p);

  biogeo_model_t &
  set_per_region_cladogenesis_params(size_t                       region_index,
                                     const cladogenesis_params_t &p);

  biogeo_model_t &set_two_region_duplicity(bool d);

  biogeo_model_t &set_extinction(bool e);

  biogeo_model_t &set_tree_params(const tree_params_t &);

  biogeo_model_t &set_adjustment_matrix(const adjustment_matrix_t &);

  bool has_adjustment_matrix() const { return _adjustment_matrix.has_value(); }

  adjustment_matrix_t const &get_adjustment_matrix() const {
    return *_adjustment_matrix;
  }

  inline cladogenesis_params_t cladogenesis_params() const {
    return _clad_params;
  }

  inline cladogenesis_params_t cladogenesis_params(const dist_t &parent) const {
    return {.allopatry = allopatry_weight(parent),
            .sympatry  = sympatry_weight(parent),
            .copy      = copy_weight(parent),
            .jump      = jump_weight(parent)};
  }

  cladogenesis_params_t normalized_cladogenesis_params() const;

  cladogenesis_params_t
  normalized_cladogenesis_params(const dist_t &dist) const;

  size_t dispersion_count(const dist_t &dist) const;
  size_t extinction_count(const dist_t &dist) const;

  double dispersion_weight(const dist_t &dist) const;
  double extinction_weight(const dist_t &dist) const;

  double dispersion_rate(size_t from, size_t to) const;
  double dispersion_weight_for_index(const dist_t &dist, size_t index) const;

  constexpr double dispersion_weight_with_adj(const dist_t &dist) const {
    auto   inv_dist = ~dist;
    double sum      = 0.0;
    auto  &adj_mat  = _adjustment_matrix.value();
    for (size_t i = 0; i < dist.regions(); ++i) {
      if (!inv_dist[i]) { continue; }
      double dis = dispersion_rate_for_region(i);
      for (size_t j = 0; j < dist.regions(); ++j) {
        sum += dis * adj_mat.get_adjustment(j, i) * dist[j];
      }
    }
    return sum;
  }

  constexpr double dispersion_rate_for_region(size_t region_index) const {
    return _has_per_region_params && _per_region_params[region_index].rates
             ? _per_region_params[region_index].rates->dis
             : _rate_params.dis;
  }

  inline bool extinction_allowed() const { return _extinction; }

  double total_rate_weight(const dist_t &dist) const;

  size_t jump_count(const dist_t &dist) const;
  size_t allopatry_count(const dist_t &dist) const;
  size_t sympatry_count(const dist_t &dist) const;
  size_t copy_count(const dist_t &dist) const;

  constexpr double jump_weight(const dist_t &dist) const {
    if (!_has_per_region_params && _clad_params.jump == 0.0) { return 0.0; }

    if (!_has_per_region_params && !_adjustment_matrix.has_value()) {
      return jump_count(dist) * _clad_params.jump;
    }

    double sum = 0.0;

    for (size_t i = 0; i < dist.regions(); ++i) {
      if (dist[i]) { continue; }
      double rate         = _clad_params.jump;
      size_t region_index = i;
      if (!_per_region_params.empty()) {
        auto &prp = _per_region_params[i];
        rate = prp.cladogenesis ? prp.cladogenesis->jump : _clad_params.jump;
      }
      double adj = 0.0;
      if (_adjustment_matrix) {
        for (size_t j = 0; j < dist.regions(); ++j) {
          adj += _adjustment_matrix->get_adjustment(region_index, j);
        }
        sum += rate * adj;
      } else {
        sum += rate;
      }
    }

    return sum;
  }

  constexpr double copy_weight(const dist_t &dist) const {
    if (!dist.singleton()) { return 0.0; }
    if (_has_per_region_params) {
      for (size_t i = 0; i < _per_region_params.size(); ++i) {
        auto &prp = _per_region_params[i];
        if (dist[i]) { return prp.cladogenesis->copy; }
      }
    }

    return copy_count(dist) * _clad_params.copy;
  }

  constexpr double allopatry_weight(const dist_t &dist) const {
    if (dist.singleton()) { return 0.0; }
    if (!_has_per_region_params) {
      return allopatry_count(dist) * _clad_params.allopatry;
    }

    double sum = 0.0;
    for (size_t i = 0; i < _per_region_params.size(); ++i) {
      const auto &prp = _per_region_params[i];

      sum += (prp.cladogenesis ? prp.cladogenesis->allopatry
                               : _clad_params.allopatry)
           * dist[i];
    }
    return sum;
  }

  constexpr double sympatry_weight(const dist_t &dist) const {
    if (dist.singleton()) { return 0.0; }
    if (!_has_per_region_params) {
      return sympatry_count(dist) * _clad_params.sympatry;
    }

    double sum = 0.0;
    for (size_t i = 0; i < _per_region_params.size(); ++i) {
      const auto &prp = _per_region_params[i];

      sum += (prp.cladogenesis ? prp.cladogenesis->sympatry
                               : _clad_params.sympatry)
           * dist[i];
    }
    return sum;
  }

  double total_singleton_weight(const dist_t &dist) const;
  double total_nonsingleton_weight(const dist_t &dist) const;

  double total_event_weight(const dist_t &dist) const;

  double total_speciation_weight(const dist_t &dist) const;

  double adjustment_prob(size_t from, size_t to) const;

  /**
   * Returns if jumps are enabled in the current model. This is done by
   * comparing the rate to zero. If the jump rate is not optimized, then
   * everything should be fine.
   */
  constexpr bool jumps_ok() const { return _clad_params.jump != 0.0; }

  bool check_cladogenesis_params_ok(size_t region_count) const;
  bool check_per_region_params_ok(size_t region_count) const;
  bool check_ok(size_t region_count) const;

private:
  rate_params_t                      _rate_params;
  cladogenesis_params_t              _clad_params;
  std::vector<per_region_params_t>   _per_region_params;
  std::optional<tree_params_t>       _tree_params;
  std::optional<adjustment_matrix_t> _adjustment_matrix;

  bool _duplicity             = false;
  bool _extinction            = false;
  bool _has_per_region_params = false;
  bool _has_adj_matrix        = false;
};
} // namespace bigrig
