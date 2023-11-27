#pragma once

#include <cstdint>
#include <utility>
#include <vector>

namespace biogeosim {

class dist_t;

struct rate_params_t {
  double dis;
  double ext;
};

struct cladogenesis_params_t {
  double copy;
  double sympatry;
  double allopatry;
  double jump;
};

class substitution_model_t {
public:
  substitution_model_t() = default;

  substitution_model_t(double d, double e, size_t r)
      : _rate_params{.dis = d, .ext = e},
        _clad_params{
            .copy = 1.0, .sympatry = 1.0, .allopatry = 1.0, .jump = 0.0},
        _region_count{r} {}

  /**
   * Returns the pair (e,d) for a simple 2 parameter dec model.
   */
  inline rate_params_t rates() const { return _rate_params; }

  inline size_t region_count() const { return _region_count; }

  double compute_denominator(size_t active_regions) const;

  substitution_model_t &set_params(rate_params_t p) {
    _rate_params = p;
    return *this;
  }

  substitution_model_t &set_params(double d, double e) {
    return set_params({.dis = d, .ext = e});
  }

  substitution_model_t &set_region_count(size_t regions) {
    _region_count = regions;
    return *this;
  }

  substitution_model_t &
  set_cladogenesis_params(double y, double s, double v, double j) {
    _clad_params.copy      = y;
    _clad_params.sympatry  = s;
    _clad_params.allopatry = v;
    _clad_params.jump      = j;

    return *this;
  }

  cladogenesis_params_t cladogenesis_params() const { return _clad_params; }

  cladogenesis_params_t normalized_cladogenesis_params() const;

  size_t jump_count(const dist_t &dist) const;
  size_t allopatry_count(const dist_t &dist) const;
  size_t sympatry_count(const dist_t &dist) const;
  size_t copy_count(const dist_t &dist) const;

  double jump_weight(const dist_t &dist) const;
  double allopatry_weight(const dist_t &dist) const;
  double sympatry_weight(const dist_t &dist) const;
  double copy_weight(const dist_t &dist) const;

  uint64_t valid_mask() const {
    uint64_t mask = 0;
    for (size_t i = 0; i < _region_count; ++i) { mask |= 1ul << i; }
    return mask;
  }

  bool jumps_ok() const { return _clad_params.jump != 0.0; }

private:
  rate_params_t         _rate_params;
  cladogenesis_params_t _clad_params;

  bool _duplicity = true;

  size_t _region_count;
};
} // namespace biogeosim
