#pragma once
#include <cstdint>
#include <utility>
#include <vector>

namespace biogeosim {

struct rate_params_t {
  double dis;
  double ext;
};

class substitution_model_t {
public:
  substitution_model_t() = default;

  substitution_model_t(double d, double e, size_t r)
      : _rate_params{.dis = d, .ext = e},
        _splitting_prob{0.5},
        _region_count{r} {}

  /**
   * Returns the pair (e,d) for a simple 2 parameter dec model.
   */
  rate_params_t rates() const { return _rate_params; }

  double splitting_prob() const { return _splitting_prob; }

  size_t region_count() const { return _region_count; }

  double compute_denominator(size_t active_regions) const {
    if (active_regions == 1) [[unlikely]] {
      return (region_count() - active_regions) * _rate_params.dis;
    }
    return active_regions * _rate_params.ext
         + (region_count() - active_regions) * _rate_params.dis;
  }

  substitution_model_t &set_params(rate_params_t p) {
    _rate_params = p;
    return *this;
  }

  substitution_model_t &set_params(double d, double e) {
    return set_params({.dis = d, .ext = e});
  }

  substitution_model_t &set_splitting_prob(double prob) {
    _splitting_prob = prob;
    return *this;
  }

  substitution_model_t &set_region_count(size_t regions) {
    _region_count = regions;
    return *this;
  }

  uint64_t valid_mask() const {
    uint64_t mask = 0;
    for (size_t i = 0; i < _region_count; ++i) { mask |= 1ul << i; }
    return mask;
  }

private:
  rate_params_t _rate_params;
  double        _splitting_prob; // Probability of Allopatry
  size_t        _region_count;
};
} // namespace biogeosim
