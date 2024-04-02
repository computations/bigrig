#pragma once

#include "model.hpp"

#include <memory>

namespace bigrig {

class period_t {
public:
  period_t()                 = default;
  period_t(const period_t &) = default;
  period_t(double                                 start,
           double                                 length,
           const std::shared_ptr<biogeo_model_t> &model)
      : _start{start}, _length{length}, _model{model} {}

  double start() const { return _start; }
  double length() const { return _length; }
  double end() const { return start() + length(); }

  const biogeo_model_t &model() const { return *_model; }

  std::shared_ptr<biogeo_model_t> model_ptr() const { return _model; }

private:
  double                          _start;
  double                          _length;
  std::shared_ptr<biogeo_model_t> _model;
};
} // namespace bigrig
