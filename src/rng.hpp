#pragma once

#include "pcg_extras.hpp"
#include "pcg_random.hpp"

#include <memory>
#include <random>

namespace bigrig {
class rng_wrapper_t {
public:
  static rng_wrapper_t &get_instance() {
    static rng_wrapper_t instance;

    return instance;
  }

  static pcg64_fast &rng() { return *get_instance()._rng; }

  static void seed() {
    rng().seed(pcg_extras::seed_seq_from<std::random_device>{});
  }

  static void seed(uint64_t seed) { rng().seed(seed); }

private:
  rng_wrapper_t() { _rng = std::make_unique<pcg64_fast>(); }

  std::unique_ptr<pcg64_fast> _rng;
};
} // namespace bigrig
