#pragma once

#include "model.hpp"

#include <logger.hpp>
#include <memory>
#include <ranges>
#include <vector>

namespace bigrig {

struct period_params_t {
  bigrig::rate_params_t                rates;
  bigrig::cladogenesis_params_t        clado;
  std::optional<bigrig::tree_params_t> tree;
  double                               start;
  std::optional<bool>                  extinction;
};

class period_t {
public:
  period_t() : _start{0.0}, _length{0.0}, _model{nullptr}, _index{0} {}
  period_t(const period_t &other)
      : _start{other._start},
        _length{other._length},
        _model{other._model ? std::make_unique<biogeo_model_t>(*other._model)
                            : nullptr},
        _index{other._index} {}

  period_t &operator=(const period_t &other) {
    _start  = other._start;
    _length = other._length;
    _model  = other._model ? std::make_unique<biogeo_model_t>(*other._model)
                           : nullptr;
    _index  = other._index;

    return *this;
  }

  period_t(double                start,
           double                length,
           const biogeo_model_t &model,
           size_t                index)
      : _start{start},
        _length{length},
        _model{new biogeo_model_t{model}},
        _index{index} {}

  period_t(double                                 start,
           double                                 length,
           const std::unique_ptr<biogeo_model_t> &model,
           size_t                                 index)
      : period_t{start, length, *model, index} {}

  period_t(double                                 start,
           double                                 length,
           const std::unique_ptr<biogeo_model_t> &model)
      : period_t{start, length, *model, 0} {}

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
  void   set_end(double e) { _length = e - _start; }
  size_t index() const { return _index; }

  void clamp(double s, double e) {
    if (e < s) {
      throw std::runtime_error{
          "Cannot clamp with an end time earlier than the start time"};
    }
    if (start() < s) { adjust_start(s); }
    if (end() > e) { set_end(e); }
  }

  const biogeo_model_t &model() const { return *_model; }

  std::unique_ptr<biogeo_model_t> const &model_ptr() const { return _model; }

  void set_model(const biogeo_model_t &model) {
    _model = std::make_unique<biogeo_model_t>(model);
  }

  void set_model(const std::unique_ptr<biogeo_model_t> &model) {
    set_model(*model);
  }

private:
  double                          _start;
  double                          _length;
  std::unique_ptr<biogeo_model_t> _model;
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

  period_list_t(const std::vector<period_params_t> &params) {
    size_t index = 0;
    for (auto &param : params) {
      biogeo_model_t model{};
      model.set_rate_params(param.rates);
      model.set_cladogenesis_params(param.clado);
      if (param.tree) { model.set_tree_params(param.tree.value()); }
      model.set_two_region_duplicity(false);
      if (param.extinction) { model.set_extinction(param.extinction.value()); }

      _periods.emplace_back(param.start, 0, model, index);
      index++;
    }

    for (size_t i = 0; i < _periods.size() - 1; ++i) {
      auto &p1 = _periods[i];
      auto &p2 = _periods[i + 1];

      p1.set_length(p2.start() - p1.start());
    }
    _periods.back().set_length(std::numeric_limits<double>::infinity());
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
      if (p.start() <= d && d < p.end()) { return p; }
    }
    throw std::runtime_error{"Period not found"};
  }

  bool validate(size_t region_count) const {
    bool ok = true;

    for (auto &p : _periods) {
      if (!p.model().check_ok(region_count)) {
        LOG_ERROR("There is an issue with the model for period '%lu', we can't "
                  "continue",
                  p.index());
      }
    }

    return ok;
  }

  void set_extinction(bool e_ok) {
    for (auto &p : _periods) { p.model_ptr()->set_extinction(e_ok); }
  }

  size_t size() const { return _periods.size(); }

private:
  void clamp_periods(double start, double end) {
    for (auto &p : _periods) { p.clamp(start, end); }
  }
  std::vector<period_t> _periods;
};

static_assert(std::ranges::range<period_list_t>);
} // namespace bigrig
