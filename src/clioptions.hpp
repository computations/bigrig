#pragma once
#include "dist.hpp"

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <filesystem>
#include <logger.hpp>
#include <nlohmann/json.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <yaml-cpp/yaml.h>

/**
 * Type-safe enum for the output file format.
 */
enum class output_format_type_e { JSON, YAML };

class cli_option_missing_required_yaml_option : std::invalid_argument {
public:
  cli_option_missing_required_yaml_option(const std::string &msg)
      : std::invalid_argument{msg} {}
};

class cli_option_invalid_parameter : std::invalid_argument {
public:
  cli_option_invalid_parameter(const std::string &msg)
      : std::invalid_argument{msg} {}
};

/**
 * Semi-smart struct containing the parsed program options. Normally these are
 * passed on the command line, but they can also be provided via a config file.
 * Includes some logic relating to merging CLI and file options.
 *
 * Many of the member variables are wrapped in a `std::optional`. This is
 * because CLI11, the cli parsing library, infers that a switch is optional
 * based on the type. In order to support both a CLI interface as well as a
 * config file, we need to allow for all of the CLI options to be optional.
 *
 * Its possible that we could do this via some other method with CLI11, but
 * honestly it's easier to handle the merging logic myself, rather than fight
 * with CLI11.
 */
struct cli_options_t {
  /**
   * Path to file containing config.
   */
  std::optional<std::filesystem::path> config_filename;

  /**
   * Path to a file containing a newick tree which will be used for simulations.
   */
  std::optional<std::filesystem::path> tree_filename;

  /**
   * Prefix used for output files. If its not specified on the CLI, it will be
   * instead be the tree filename.
   */
  std::optional<std::filesystem::path> prefix;

  /**
   * Flag to enable a debug log, which will be very verbose.
   */
  std::optional<bool> debug_log;

  /**
   * Output format enum. Right now valid options are
   * - JSON
   * - YAML
   */
  std::optional<output_format_type_e> output_format_type;

  /**
   * Starting distribution for the simulation.
   */
  std::optional<bigrig::dist_t> root_distribution;

  /**
   * Dispersion rate which was provided by the user.
   */
  std::optional<double> dispersion_rate;

  /**
   * Extinction rate which was provided by the user.
   */
  std::optional<double> extinction_rate;

  /**
   * Allopatry rate for cladogenesis which was provide by the user.
   */
  std::optional<double> allopatry_rate;

  /**
   * Sympatry rate for cladogenesis which was provide by the user.
   */
  std::optional<double> sympatry_rate;

  /**
   * Copy rate for cladogenesis which was provide by the user.
   */
  std::optional<double> copy_rate;

  /**
   * Jump rate for cladogenesis which was provide by the user.
   */
  std::optional<double> jump_rate;

  /**
   * Enable overwriting the exsisting result files.
   */
  std::optional<bool> redo;

  /**
   * There is a case that I don't know how to solve, so I decided to offer it as
   * an option. In the case of 2 regions, there are 2 ways to count the number
   * of splits. One of them "double counts", because the method of production
   * produces equivalent results. I didn't quite know what to do, so I decided
   * to support both ways of counting.
   */
  std::optional<bool> two_region_duplicity;

  std::optional<bigrig::operation_mode_e> mode;

  std::optional<uint64_t> rng_seed;

  std::filesystem::path phylip_filename() const;

  std::filesystem::path yaml_filename() const;

  std::filesystem::path json_filename() const;

  bool cli_arg_specified() const;

  bool yaml_file_set() const;

  bool json_file_set() const;

  void merge(const cli_options_t &other);

  cli_options_t() = default;

  /**
   * Construct a cli_option_t from a YAML file. Please see the YAML file schema
   * for the details.
   */
  cli_options_t(const YAML::Node &yaml)
      : tree_filename{get_tree_filename(yaml)},
        prefix{get_prefix(yaml)},
        debug_log{get_debug_log(yaml)},
        output_format_type{get_output_format(yaml)},
        root_distribution{get_root_range(yaml)},
        redo{get_redo(yaml)},
        two_region_duplicity{get_two_region_duplicity(yaml)},
        mode{get_mode(yaml)},
        rng_seed{get_seed(yaml)} {
    std::tie(dispersion_rate, extinction_rate) = get_rates(yaml);
    std::tie(allopatry_rate, sympatry_rate, copy_rate, jump_rate)
        = get_cladogenesis(yaml);
  }

private:
  static std::filesystem::path get_tree_filename(const YAML::Node &);
  static std::optional<std::filesystem::path> get_prefix(const YAML::Node &);
  static std::optional<bool>                  get_debug_log(const YAML::Node &);
  static std::optional<output_format_type_e>
                             get_output_format(const YAML::Node &);
  static std::optional<bool> get_two_region_duplicity(const YAML::Node &);
  static std::optional<bigrig::dist_t> get_root_range(const YAML::Node &);
  static std::tuple<std::optional<double>, std::optional<double>>
  get_rates(const YAML::Node &yaml);

  static std::tuple<std::optional<double>,
                    std::optional<double>,
                    std::optional<double>,
                    std::optional<double>>
  get_cladogenesis(const YAML::Node &yaml);

  static std::optional<bool>                     get_redo(const YAML::Node &);
  static std::optional<bigrig::operation_mode_e> get_mode(const YAML::Node &);
  static std::optional<uint64_t>                 get_seed(const YAML::Node &);
};
