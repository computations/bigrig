#pragma once

#include <cstddef>
#include <cstdint>
#include <sstream>
#include <utility>
#include <vector>

namespace bigrig {

class dist_t;

struct rate_params_t {
  double dis;
  double ext;
};

/**
 * This data structure is a bit wacky. It really is just an array with named
 * elements (which the compiler is aware of. Sometimes I want to perform
 * operations on this as if it was an array. So, there are some WACKY functions
 * which basically cast the structure to an array, and iterate over the array to
 * compute a sum. Similarly, we use the same technique to normalize the
 * parameters.
 */
struct cladogenesis_params_t {
  using data_type = double;

  data_type copy;
  data_type sympatry;
  data_type allopatry;
  data_type jump;

  static constexpr size_t size() {
    static_assert(sizeof(cladogenesis_params_t) % sizeof(data_type) == 0,
                  "something is fucked");
    return sizeof(cladogenesis_params_t) / sizeof(data_type);
  }

  data_type sum() const {
    auto   start = as_ptr();
    double acc   = 0;
    for (size_t i = 0; i < size(); ++i) { acc += start[i]; }
    return acc;
  }

  cladogenesis_params_t normalize() const {
    cladogenesis_params_t tmp{*this};

    tmp.normalize(sum());
    return tmp;
  }

  std::string to_debug_string() const {
    std::ostringstream oss;
    oss << "copy: " << copy << ", sympatry: " << sympatry
        << ", allopatry: " << allopatry << ", jump: " << jump;
    return oss.str();
  }

private:
  using array_type       = data_type *;
  using const_array_type = data_type const *;

  const_array_type as_ptr() const {
    return reinterpret_cast<const_array_type>(this);
  }

  array_type as_ptr() { return reinterpret_cast<array_type>(this); }

  void normalize(data_type d) {
    auto start = as_ptr();
    for (size_t i = 0; i < size(); ++i) { start[i] /= d; }
  }
};

/* Assert that the fields are in the right order */
static_assert(offsetof(cladogenesis_params_t, copy)
                  / sizeof(cladogenesis_params_t::data_type)
              == 0);
static_assert(offsetof(cladogenesis_params_t, sympatry)
                  / sizeof(cladogenesis_params_t::data_type)
              == 1);
static_assert(offsetof(cladogenesis_params_t, allopatry)
                  / sizeof(cladogenesis_params_t::data_type)
              == 2);
static_assert(offsetof(cladogenesis_params_t, jump)
                  / sizeof(cladogenesis_params_t::data_type)
              == 3);

/* Check that jump is the last field */
static_assert(offsetof(cladogenesis_params_t, jump)
                  + sizeof(cladogenesis_params_t::data_type)
              == sizeof(cladogenesis_params_t));

class substitution_model_t {
public:
  substitution_model_t() = default;

  substitution_model_t(double d, double e, size_t r, bool duplicity)
      : _rate_params{.dis = d, .ext = e},
        _clad_params{
            .copy = 1.0, .sympatry = 1.0, .allopatry = 1.0, .jump = 0.0},
        _duplicity{duplicity},
        _region_count{r} {}

  /**
   * Returns the pair (e,d) for a simple 2 parameter dec model.
   */
  inline rate_params_t rates() const { return _rate_params; }

  inline size_t region_count() const { return _region_count; }

  double compute_denominator(size_t active_regions) const;

  substitution_model_t &set_params(rate_params_t p);

  substitution_model_t &set_params(double d, double e);

  substitution_model_t &set_region_count(size_t regions);

  substitution_model_t &
  set_cladogenesis_params(double y, double s, double v, double j);

  substitution_model_t &set_cladogenesis_params(const cladogenesis_params_t &p);

  substitution_model_t &set_two_region_duplicity(bool d);

  cladogenesis_params_t cladogenesis_params() const { return _clad_params; }

  cladogenesis_params_t normalized_cladogenesis_params() const;

  size_t jump_count(const dist_t &dist) const;
  size_t allopatry_count(const dist_t &dist) const;
  size_t sympatry_count(const dist_t &dist) const;
  size_t copy_count(const dist_t &dist) const;

  double jump_weight(const dist_t &dist) const;
  double allopatry_weight(const dist_t &dist) const;
  double sympatry_weight(const dist_t &dist) const;
  double copy_weight(const dist_t &dist) const;

  constexpr uint64_t valid_mask() const {
    uint64_t mask = 0;
    for (size_t i = 0; i < _region_count; ++i) { mask |= 1ul << i; }
    return mask;
  }

  constexpr bool jumps_ok() const { return _clad_params.jump != 0.0; }

private:
  rate_params_t         _rate_params;
  cladogenesis_params_t _clad_params;

  bool _duplicity = true;

  size_t _region_count;
};
} // namespace bigrig
