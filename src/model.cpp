#include "model.hpp"

#include "dist.hpp"

namespace bigrig {

double substitution_model_t::compute_denominator(size_t active_regions) const {
  if (active_regions == 1) [[unlikely]] {
    return (region_count() - active_regions) * _rate_params.dis;
  }
  return active_regions * _rate_params.ext
       + (region_count() - active_regions) * _rate_params.dis;
}

cladogenesis_params_t
substitution_model_t::normalized_cladogenesis_params() const {
  auto   tmp         = _clad_params;
  double denominator = _clad_params.sum();

  tmp.copy      /= denominator;
  tmp.sympatry  /= denominator;
  tmp.allopatry /= denominator;
  tmp.jump      /= denominator;

  return _clad_params;
}

size_t substitution_model_t::jump_count(const dist_t &dist) const {
  return (dist.regions() - dist.popcount()) * 2;
}

size_t substitution_model_t::allopatry_count(const dist_t &dist) const {
  return dist.popcount() * 2 - _duplicity * (dist.popcount() == 2) * 2;
}

size_t substitution_model_t::sympatry_count(const dist_t &dist) const {
  return dist.popcount() * 2;
}

size_t substitution_model_t::copy_count(const dist_t &dist) const {
  return (dist.popcount() == 1) * (2 * _duplicity);
}

double substitution_model_t::jump_weight(const dist_t &dist) const {
  return jump_count(dist) * _clad_params.jump;
}

double substitution_model_t::sympatry_weight(const dist_t &dist) const {
  return sympatry_count(dist) * _clad_params.sympatry;
}

double substitution_model_t::allopatry_weight(const dist_t &dist) const {
  return allopatry_count(dist) * _clad_params.allopatry;
}

double substitution_model_t::copy_weight(const dist_t &dist) const {
  return copy_count(dist) * _clad_params.copy;
}

substitution_model_t &substitution_model_t::set_params(rate_params_t p) {
  _rate_params = p;
  return *this;
}
substitution_model_t &substitution_model_t::set_params(double d, double e) {
  return set_params({.dis = d, .ext = e});
}

substitution_model_t &substitution_model_t::set_region_count(size_t regions) {
  _region_count = regions;
  return *this;
}

substitution_model_t &substitution_model_t::set_cladogenesis_params(double y,
                                                                    double s,
                                                                    double v,
                                                                    double j) {
  _clad_params.copy      = y;
  _clad_params.sympatry  = s;
  _clad_params.allopatry = v;
  _clad_params.jump      = j;

  return *this;
}

substitution_model_t &
substitution_model_t::set_cladogenesis_params(const cladogenesis_params_t &p) {
  _clad_params = p;

  return *this;
}

substitution_model_t &substitution_model_t::set_two_region_duplicity(bool d) {
  _duplicity = d;
  return *this;
}

} // namespace bigrig
