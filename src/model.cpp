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

/**
 * Count the number of jump splits possible.
 *
 * A jump split is
 *
 *     10
 *    /  \
 *   01  10
 *
 * In this case, it is the number of empty regions times 2. Because the region
 * that is "jump started" can be either on the left or right, and the left and
 * right branches are distinct.
 */
size_t biogeo_model_t::jump_count(const dist_t &dist) const {
  return dist.empty_region_count() * 2;
}

/**
 * Count the number of allopatric splits for a given dist.
 *
 * An allopatric split is
 *
 *     111
 *    /   \
 *  011   100
 *
 * The number of allopatric splits is 2 times the number of full regions. The
 * idea is that we pick a full region to split, and then pick which branch the
 * split gets assigned to.
 *
 * Because the split in this case is disjoint, there is a weird duplicity case
 * when there are 2 regions. If you count by _outcomes_ you get two less
 * possible allopatric splits the if you count by _process_. By default, we
 * don't allow this as this is what was done by Matzke for DEC+J. However, we do
 * support the other method.
 */
size_t biogeo_model_t::allopatry_count(const dist_t &dist) const {
  return dist.full_region_count() * 2
       - ((!_duplicity && dist.full_region_count() == 2) ? 2 : 0);
}

/**
 * Compute the number of sympatric splits
 *
 * A sympatric split is
 *
 *     11
 *    /  \
 *   11  10
 *
 * In this case, the number of splits is 2 times the full splits. The idea is
 * that we pick a full region to split, and then pick the left or right branch.
 */
size_t biogeo_model_t::sympatry_count(const dist_t &dist) const {
  return dist.full_region_count() * 2;
}

/**
 * Compute the number of copy splits
 *
 * A copy split is
 *
 *     10
 *    /  \
 *   10  10
 *
 * Copies can only occur on singleton regions. Again, there is a duplicity issue
 * here, though I am less confident here. In one sense, there is no process, so
 * the left and right branches being distinguished doesn't matter. In another
 * sense, we distinguish the branches in all other cases, so why should we stop
 * for this one.
 */
size_t biogeo_model_t::copy_count(const dist_t &dist) const {
  return (dist.full_region_count() == 1) * (_duplicity ? 1 : 2);
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

double biogeo_model_t::total_singleton_weight(const dist_t &dist) const {
  return copy_weight(dist) + jump_weight(dist);
}
double biogeo_model_t::total_nonsingleton_weight(const dist_t &dist) const {
  return sympatry_weight(dist) + allopatry_weight(dist) + jump_weight(dist);
}

biogeo_model_t &biogeo_model_t::set_params(rate_params_t p) {
  _rate_params = p;
  return *this;
}
biogeo_model_t &biogeo_model_t::set_params(double d, double e) {
  return set_params({.dis = d, .ext = e});
}

biogeo_model_t &biogeo_model_t::set_cladogenesis_params(double v,
                                                        double s,
                                                        double y,
                                                        double j) {
  _clad_params = {.allopatry = v, .sympatry = s, .copy = y, .jump = j};

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

bool biogeo_model_t::check_cladogenesis_params_ok(size_t region_count) const {
  bool ok = true;
  if (total_nonsingleton_weight(make_full_dist(region_count)) == 0.0) {
    MESSAGE_ERROR("The sympatry, allopatry, or jump weights are invalid");
    ok = false;
  }
  if (total_singleton_weight(make_singleton_dist(region_count)) == 0.0) {
    MESSAGE_ERROR("The copy or jump weights are invalid");
    ok = false;
  }

  return ok;
}

bool biogeo_model_t::check_ok(size_t region_count) const {
  return check_cladogenesis_params_ok(region_count);
}
} // namespace bigrig
