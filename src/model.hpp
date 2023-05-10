#pragma once
#include <utility>
#include <vector>

class biogeosim_substitution_model_t {
public:
  /**
   * Returns the pair (e,d) for a simple 2 parameter dec model.
   */
  std::pair<double, double> rates() const { return {_ext, _dis}; }

  size_t region_count() const { return _region_count; }

  double compute_denominator(size_t active_regions) const {
    if (active_regions == 1) {
      return (region_count() - active_regions) * _dis;
    }
    return active_regions * _ext + (region_count() - active_regions) * _dis;
  }

private:
  double _ext;
  double _dis;
  size_t _region_count;
};
