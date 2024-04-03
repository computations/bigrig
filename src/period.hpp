#pragma once

#include "model.hpp"

#include <memory>

namespace bigrig {

class period_t {
public:
  period_t()                            = default;
  period_t(const period_t &)            = default;
  period_t &operator=(const period_t &) = default;
  period_t(double                       start,
           double                       length,
           const rate_params_t         &rp,
           const cladogenesis_params_t &cp,
           bool                         two_region_duplicity,
           size_t                       index)
      : _start{start},
        _length{length},
        _model{std::make_shared<biogeo_model_t>(rp, cp, two_region_duplicity)},
        _index{index} {}

  period_t(double                                 start,
           double                                 length,
           const std::shared_ptr<biogeo_model_t> &model)
      : _start{start}, _length{length}, _model{model} {}

  double start() const { return _start; }
  double length() const { return _length; }
  double end() const { return start() + length(); }
  void   set_length(double l) { _length = l; }
  void   set_start(double s) { _start = s; }
  void   adjust_start(double s) {
    double new_length = _length - (s - _start);
    _start            = s;
    _length           = new_length;
  }
  void adjust_end(double e) { _length = e - _start; }

  const biogeo_model_t &model() const { return *_model; }

  std::shared_ptr<biogeo_model_t> model_ptr() const { return _model; }

private:
  double                          _start;
  double                          _length;
  std::shared_ptr<biogeo_model_t> _model;
  size_t                          _index;
};
} // namespace bigrig
