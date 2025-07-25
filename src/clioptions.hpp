#pragma once

#include "dist.hpp"
#include "errors.hpp"
#include "model.hpp"
#include "period.hpp"
#include "rng.hpp"

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <chrono>
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
enum class output_format_type_e { JSON, YAML, CSV };

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

struct program_stats_t {
  using time_point_t
      = std::chrono::time_point<std::chrono::high_resolution_clock>;
  double execution_time_in_seconds() const {
    return std::chrono::duration<double>(end_time - config_time).count();
  }

  double config_time_in_seconds() const {
    return std::chrono::duration<double>(config_time - start_time).count();
  }

  double total_time_in_seconds() const {
    return std::chrono::duration<double>(end_time - start_time).count();
  }

  time_point_t start_time;
  time_point_t config_time;
  time_point_t end_time;
};

[[nodiscard]] bool verify_path_is_readable(const std::filesystem::path &path);

[[nodiscard]] bool verify_path_is_writable(const std::filesystem::path &path);

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
  std::optional<bigrig::dist_t> root_range;

  /**
   * Number of regions to simulate.
   *
   * Optional, should be specified if no root_distribution is present.
   */
  std::optional<size_t> region_count;

  /**
   * Name of the regions. 
   *
   * If not given, region names will be generated as needed.
   */
  std::optional<std::vector<std::string>> region_names;

  /**
   * Rates which have been provided by the user.
   */
  std::vector<bigrig::period_params_t> periods;

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

  /**
   * Option to simulate the tree _alongside_ the ranges.
   */
  std::optional<bool> simulate_tree;

  /**
   * Height of the simmulated tree.
   */
  std::optional<double> tree_height;

  std::filesystem::path phylip_filename() const;

  std::filesystem::path yaml_filename() const;

  std::filesystem::path json_filename() const;

  std::filesystem::path csv_splits_filename() const;
  std::filesystem::path csv_events_filename() const;
  std::filesystem::path csv_periods_filename() const;
  std::filesystem::path csv_matrix_filename() const;
  std::filesystem::path csv_region_names_filename() const;
  std::filesystem::path csv_program_stats_filename() const;

  std::vector<std::filesystem::path> csv_file_vector() const;

  std::vector<std::filesystem::path> result_filename_vector() const;

  pcg64_fast            &get_rng();
  bigrig::rng_wrapper_t &get_rng_wrapper();

  size_t compute_region_count() const {
    if (root_range.value()) { return root_range->regions(); }
    if (region_count && region_names) {
      if (region_count != region_names->size()) {
        MESSAGE_ERROR(
            "The number of regions provided in names differs from the "
            "region count");
        throw std::runtime_error{"Failed to compute the region count"};
      }
      return region_count.value();
    } else if (region_count || region_names) {
      return region_count ? region_count.value() : region_names->size();
    }

    MESSAGE_ERROR("There was an issue with the root region");
    throw std::runtime_error{"Failed to compute the region count"};
  }

  bool cli_arg_specified() const;

  bool yaml_file_set() const;

  bool json_file_set() const;

  bool csv_file_set() const;

  void merge(const cli_options_t &other);

  [[nodiscard]] bool convert_cli_parameters(std::optional<double> dis,
                                            std::optional<double> ext,
                                            std::optional<double> allo,
                                            std::optional<double> symp,
                                            std::optional<double> copy,
                                            std::optional<double> jump);

  bigrig::period_list_t
  make_periods(std::uniform_random_bit_generator auto &gen) const {
    return {periods, *region_names, gen};
  }

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
        root_range{get_root_range(yaml)},
        region_count{get_region_count(yaml)},
        region_names{get_region_names(yaml)},
        periods{get_periods(yaml)},
        redo{get_redo(yaml)},
        two_region_duplicity{get_two_region_duplicity(yaml)},
        mode{get_mode(yaml)},
        rng_seed{get_seed(yaml)},
        simulate_tree{get_simulate_tree(yaml)},
        tree_height{get_tree_height(yaml)} {}

private:
  static std::optional<std::filesystem::path>
  get_tree_filename(const YAML::Node &);

  static std::optional<std::filesystem::path> get_prefix(const YAML::Node &);
  static std::optional<bool>                  get_debug_log(const YAML::Node &);
  static std::optional<output_format_type_e>
                             get_output_format(const YAML::Node &);
  static std::optional<bool> get_two_region_duplicity(const YAML::Node &);
  static std::optional<bigrig::dist_t> get_root_range(const YAML::Node &);
  static std::optional<size_t>         get_region_count(const YAML::Node &yaml);

  static std::optional<std::vector<std::string>>
  get_region_names(const YAML::Node &);

  static std::optional<bigrig::rate_params_t> get_rates(const YAML::Node &yaml);

  static std::optional<bigrig::cladogenesis_params_t>
  get_cladogenesis(const YAML::Node &yaml);

  static std::optional<bigrig::period_params_t>
  get_period(const YAML::Node &yaml);

  static std::vector<bigrig::period_params_t>
  get_periods(const YAML::Node &yaml);

  static std::optional<bool>                     get_redo(const YAML::Node &);
  static std::optional<bigrig::operation_mode_e> get_mode(const YAML::Node &);
  static std::optional<uint64_t>                 get_seed(const YAML::Node &);

  static std::optional<bool>   get_simulate_tree(const YAML::Node &);
  static std::optional<double> get_tree_height(const YAML::Node &);

  static std::optional<bigrig::adjustment_matrix_params_t>
  get_adjustment_matrix_parameters(const YAML::Node &yaml);

  static std::optional<bigrig::adjustment_matrix_t>
  get_adjustment_matrix(const YAML::Node &yaml);

  static std::expected<std::filesystem::path, bigrig::io_err>
  get_adjustment_matrix_filename(const YAML::Node &);

  static std::optional<double> get_distance_exponent(const YAML::Node &);

  static std::optional<bool> get_simulate_adjustment_matrix(const YAML::Node &);

  std::vector<std::string> generate_region_names();
};
