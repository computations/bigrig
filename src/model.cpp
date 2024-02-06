#include "model.hpp"

#include "dist.hpp"

namespace bigrig {

cladogenesis_params_t biogeo_model_t::normalized_cladogenesis_params() const {
  return _clad_params.normalize();
}

size_t biogeo_model_t::dispersion_count(const dist_t &dist) const {
  return dist.empty_region_count();
}

size_t biogeo_model_t::extinction_count(const dist_t &dist) const {
  return dist.full_region_count();
}

double biogeo_model_t::dispersion_weight(const dist_t &dist) const {
  return _rate_params.dis * dispersion_count(dist);
}

double biogeo_model_t::extinction_weight(const dist_t &dist) const {
  return dist.singleton() ? 0.0 : _rate_params.ext * extinction_count(dist);
}

double biogeo_model_t::total_rate_weight(const dist_t &dist) const {
  return extinction_weight(dist) + dispersion_weight(dist);
}

size_t biogeo_model_t::jump_count(const dist_t &dist) const {
  return dist.empty_region_count() * 2;
}

size_t biogeo_model_t::allopatry_count(const dist_t &dist) const {
  return dist.full_region_count() * 2
       - _duplicity * (dist.full_region_count() == 2) * 2;
}

size_t biogeo_model_t::sympatry_count(const dist_t &dist) const {
  return dist.full_region_count() * 2;
}

size_t biogeo_model_t::copy_count(const dist_t &dist) const {
  return (dist.full_region_count() == 1) * (2 * _duplicity);
}

double biogeo_model_t::jump_weight(const dist_t &dist) const {
  return jump_count(dist) * _clad_params.jump;
}

double biogeo_model_t::sympatry_weight(const dist_t &dist) const {
  return sympatry_count(dist) * _clad_params.sympatry;
}

double biogeo_model_t::allopatry_weight(const dist_t &dist) const {
  return allopatry_count(dist) * _clad_params.allopatry;
}

double biogeo_model_t::copy_weight(const dist_t &dist) const {
  return copy_count(dist) * _clad_params.copy;
}

biogeo_model_t &biogeo_model_t::set_params(rate_params_t p) {
  _rate_params = p;
  return *this;
}
biogeo_model_t &biogeo_model_t::set_params(double d, double e) {
  return set_params({.dis = d, .ext = e});
}

biogeo_model_t &biogeo_model_t::set_region_count(size_t regions) {
  _region_count = regions;
  return *this;
}

biogeo_model_t &biogeo_model_t::set_cladogenesis_params(double y,
                                                        double s,
                                                        double v,
                                                        double j) {
  _clad_params.copy      = y;
  _clad_params.sympatry  = s;
  _clad_params.allopatry = v;
  _clad_params.jump      = j;

  return *this;
}

biogeo_model_t &
biogeo_model_t::set_cladogenesis_params(const cladogenesis_params_t &p) {
  _clad_params = p;

  return *this;
}

biogeo_model_t &biogeo_model_t::set_two_region_duplicity(bool d) {
  _duplicity = d;
  return *this;
}

} // namespace bigrig
