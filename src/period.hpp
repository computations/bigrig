#pragma once

#include "model.hpp"

#include <memory>
#include <ranges>
#include <vector>

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
  void   adjust_end(double e) { _length = e - _start; }
  size_t index() const { return _index; }

  void clamp(double s, double e) {
    if (start() < s) { adjust_start(s); }
    if (end() > e) { adjust_end(e); }
  }

  const biogeo_model_t &model() const { return *_model; }

  std::shared_ptr<biogeo_model_t> model_ptr() const { return _model; }

private:
  double                          _start;
  double                          _length;
  std::shared_ptr<biogeo_model_t> _model;
  size_t                          _index;
};

class period_list_t {
public:
  period_list_t() = default;

  period_list_t(const std::vector<period_t> &periods) : _periods{periods} {}

  period_list_t(const period_list_t &other, double start, double end) {
    auto period_filter = [&](const period_t &p) {
      return !((p.end() < start) && (p.start() > end));
    };

    for (auto &p : other._periods | std::ranges::views::filter(period_filter)) {
      _periods.push_back(p);
    }
    clamp_periods(start, end);
  }

  using period_iter       = std::vector<period_t>::iterator;
  using const_period_iter = std::vector<period_t>::const_iterator;

  period_iter       begin() { return _periods.begin(); }
  const_period_iter begin() const { return _periods.begin(); }

  period_iter       end() { return _periods.end(); }
  const_period_iter end() const { return _periods.end(); }

  period_t back() const { return _periods.back(); }

  bool empty() const { return _periods.empty(); }

  period_t get(double d) const {
    for (auto &p : _periods) {
      if (p.start() <= d && d <= p.end()) { return p; }
    }
    throw std::runtime_error{"Period not found"};
  }

private:
  void clamp_periods(double start, double end) {
    for (auto &p : _periods) { p.clamp(start, end); }
  }
  std::vector<period_t> _periods;
};

static_assert(std::ranges::range<period_list_t>);
} // namespace bigrig
