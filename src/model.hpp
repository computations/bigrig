#pragma once
#include <utility>
#include <vector>

namespace biogeosim {
class substitution_model_t {
public:
  substitution_model_t(double d, double e, size_t r)
      : _dis{d}, _ext{e}, _region_count{r} {}

  /**
   * Returns the pair (e,d) for a simple 2 parameter dec model.
   */
  std::pair<double, double> rates() const { return {_ext, _dis}; }

  size_t region_count() const { return _region_count; }

  double compute_denominator(size_t active_regions) const {
    [[assume(active_regions != 0)]];
    if (active_regions == 1) [[unlikely]] {
      return (region_count() - active_regions) * _dis;
    }
    return active_regions * _ext + (region_count() - active_regions) * _dis;
  }

private:
  double _dis;
  double _ext;
  size_t _region_count;
};
} // namespace biogeosim
