#include "clioptions.hpp"

#include "adjustment.hpp"
#include "dist.hpp"
#include "logger.hpp"
#include "period.hpp"
#include "rng.hpp"
#include "util.hpp"

#include <algorithm>
#include <expected>
#include <filesystem>
#include <string>

[[nodiscard]] bool verify_path_is_readable(const std::filesystem::path &path) {
  bool ok     = true;
  auto status = std::filesystem::status(path);
  if ((std::filesystem::perms::owner_read & status.permissions())
      == std::filesystem::perms::none) {
    LOG_ERROR("The path '{}' is not readable", path.c_str());
    ok = false;
  }
  return ok;
}

[[nodiscard]] bool verify_path_is_writable(const std::filesystem::path &path) {
  bool ok             = true;
  auto status         = std::filesystem::status(path);
  auto required_perms = std::filesystem::perms::owner_write
                      | std::filesystem::perms::owner_exec;
  if ((required_perms & status.permissions()) == std::filesystem::perms::none) {
    LOG_ERROR("The path '{}' is not writable", path.c_str());
    ok = false;
  }
  return ok;
}

template <typename T>
std::optional<T> get_yaml_val_or_nothing(const YAML::Node &yaml,
                                         const char       *KEY) {
  if (yaml[KEY]) { return yaml[KEY].as<T>(); }
  return {};
}

pcg64_fast &cli_options_t::get_rng() { return bigrig::rng_wrapper_t::rng(); }
bigrig::rng_wrapper_t &cli_options_t::get_rng_wrapper() {
  return bigrig::rng_wrapper_t::get_instance();
}

std::filesystem::path cli_options_t::phylip_filename() const {
  auto tmp  = prefix.value();
  tmp      += bigrig::util::PHYILP_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::yaml_filename() const {
  auto tmp  = prefix.value();
  tmp      += bigrig::util::YAML_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::json_filename() const {
  auto tmp  = prefix.value();
  tmp      += bigrig::util::JSON_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::csv_splits_filename() const {
  constexpr auto splits_subprefix  = ".splits";
  auto           tmp               = prefix.value();
  tmp                             += splits_subprefix;
  tmp                             += bigrig::util::CSV_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::csv_events_filename() const {
  constexpr auto events_subprefix  = ".events";
  auto           tmp               = prefix.value();
  tmp                             += events_subprefix;
  tmp                             += bigrig::util::CSV_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::csv_periods_filename() const {
  constexpr auto periods_subprefix  = ".periods";
  auto           tmp                = prefix.value();
  tmp                              += periods_subprefix;
  tmp                              += bigrig::util::CSV_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::csv_matrix_filename() const {
  constexpr auto matrix_subprefix  = ".periods.matrix";
  auto           tmp               = prefix.value();
  tmp                             += matrix_subprefix;
  tmp                             += bigrig::util::CSV_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::csv_region_names_filename() const {
  constexpr auto regions_subprefix  = ".regions";
  auto           tmp                = prefix.value();
  tmp                              += regions_subprefix;
  tmp                              += bigrig::util::CSV_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::csv_program_stats_filename() const {
  constexpr auto stats_subprefix  = ".program-stats";
  auto           tmp              = prefix.value();
  tmp                            += stats_subprefix;
  tmp                            += bigrig::util::CSV_EXT;
  return tmp;
}

std::vector<std::filesystem::path> cli_options_t::csv_file_vector() const {
  return {csv_splits_filename(),
          csv_events_filename(),
          csv_periods_filename(),
          csv_matrix_filename(),
          csv_region_names_filename(),
          csv_program_stats_filename()};
}

std::vector<std::filesystem::path>
cli_options_t::result_filename_vector() const {
  if (yaml_file_set()) { return {yaml_filename()}; }
  if (json_file_set()) { return {json_filename()}; }
  if (csv_file_set()) { return csv_file_vector(); }
  LOG_ERROR("Results files are ill-configured");
  return {};
}

/**
 * Checks if all the required args for the CLI have been specified.
 */
bool cli_options_t::cli_arg_specified() const {
  return tree_filename.has_value() || prefix.has_value()
      || debug_log.has_value() || tree_filename.has_value();
}

/**
 * Checks if the output format is a YAML file.
 */
bool cli_options_t::yaml_file_set() const {
  return output_format_type.has_value()
      && output_format_type.value() == output_format_type_e::YAML;
}

/**
 * Checks if the output format is a JSON file.
 */
bool cli_options_t::json_file_set() const {
  return output_format_type.has_value()
      && output_format_type.value() == output_format_type_e::JSON;
}

bool cli_options_t::csv_file_set() const {
  return output_format_type.has_value()
      && output_format_type.value() == output_format_type_e::CSV;
}

/**
 * Merges a `cli_options_t` with the current value. Specifically, it overwrites
 * the current values with the values from the passed `cli_options_t`. Values
 * affected are:
 *  - `config_filename`
 *  - `tree_filename`
 *  - `debug_log`
 *  - `output_format_type`
 *  - `root_distribution`
 *  - `dispersion_rate`
 *  - `extinction_rate`
 *  - `redo`
 *  - `two_region_duplicity`
 */

void print_config_cli_warning(const char *option_name) {
  LOG_WARNING("The '{}' option is specified in both the config file and "
              "the command line. Using the value from the command line",
              option_name);
}

template <typename T>
void merge_variable(std::optional<T>       &ours,
                    const std::optional<T> &theirs,
                    const char             *name) {
  if (theirs.has_value()) {
    if (ours.has_value()) {
      print_config_cli_warning(name);
    } else {
      ours = theirs;
    }
  }
}

void cli_options_t::merge(const cli_options_t &other) {
  merge_variable(tree_filename, other.tree_filename, "tree");
  merge_variable(prefix, other.prefix, "prefix");
  merge_variable(debug_log, other.debug_log, "debug-log");
  merge_variable(output_format_type, other.output_format_type, "output-format");
  merge_variable(root_range, other.root_range, "root-range");
  merge_variable(region_count, other.region_count, "region-count");

  /* periods are a bit different, so we handle them differently */
  if (!other.periods.empty()) {
    if (!periods.empty()) {
      print_config_cli_warning("periods");
    } else {
      periods = other.periods;
    }
  }

  merge_variable(redo, other.redo, "redo");
  merge_variable(
      two_region_duplicity, other.two_region_duplicity, "two-region-duplicity");
  merge_variable(mode, other.mode, "mode");
  merge_variable(rng_seed, other.rng_seed, "seed");
  merge_variable(simulate_tree, other.simulate_tree, "simulate-tree");
  merge_variable(tree_height, other.tree_height, "tree-height");
}

std::optional<std::filesystem::path>
cli_options_t::get_tree_filename(const YAML::Node &yaml) {
  constexpr auto TREE_KEY = "tree";
  return get_yaml_val_or_nothing<std::string>(yaml, TREE_KEY);
}

std::optional<std::filesystem::path>
cli_options_t::get_prefix(const YAML::Node &yaml) {
  constexpr auto PREFIX_KEY = "prefix";
  return get_yaml_val_or_nothing<std::string>(yaml, PREFIX_KEY);
}

std::optional<bool> cli_options_t::get_debug_log(const YAML::Node &yaml) {
  constexpr auto DEBUG_LOG_KEY = "debug-log";
  return get_yaml_val_or_nothing<bool>(yaml, DEBUG_LOG_KEY);
}

std::optional<output_format_type_e>
cli_options_t::get_output_format(const YAML::Node &yaml) {
  constexpr auto OUPUT_FORMAT_KEY = "output-format";
  if (yaml[OUPUT_FORMAT_KEY]) {
    auto value = yaml[OUPUT_FORMAT_KEY].as<std::string>();
    std::transform(value.begin(),
                   value.end(),
                   value.begin(),
                   [](char c) -> char { return std::tolower(c); });

    if (value == "json") { return output_format_type_e::JSON; }
    if (value == "yaml") { return output_format_type_e::YAML; }
    if (value == "csv") { return output_format_type_e::CSV; }
  }
  return {};
}

std::optional<bool>
cli_options_t::get_two_region_duplicity(const YAML::Node &yaml) {
  constexpr auto DUPLICITY_KEY = "two-region-duplicity";
  return get_yaml_val_or_nothing<bool>(yaml, DUPLICITY_KEY);
}

std::optional<bigrig::dist_t>
cli_options_t::get_root_range(const YAML::Node &yaml) {
  constexpr auto ROOT_DIST_KEY = "root-range";
  return get_yaml_val_or_nothing<std::string>(yaml, ROOT_DIST_KEY);
}

std::optional<size_t> cli_options_t::get_region_count(const YAML::Node &yaml) {
  constexpr auto REGION_COUNT_KEY = "region-count";
  return get_yaml_val_or_nothing<size_t>(yaml, REGION_COUNT_KEY);
}

std::optional<std::vector<std::string>>
cli_options_t::get_region_names(const YAML::Node &yaml) {
  constexpr auto REGION_NAMES_KEY = "region-names";
  return get_yaml_val_or_nothing<std::vector<std::string>>(yaml,
                                                           REGION_NAMES_KEY);
}

std::optional<bigrig::rate_params_t>
cli_options_t::get_rates(const YAML::Node &yaml) {
  double         dis = 0.0, ext = 0.0;
  constexpr auto RATE_DICT_KEY = "rates";
  if (yaml[RATE_DICT_KEY]) {
    auto rates = yaml[RATE_DICT_KEY];

    constexpr auto DISPERSION_KEY = "dispersion";
    if (rates[DISPERSION_KEY]) {
      dis = rates[DISPERSION_KEY].as<double>();
    } else {
      LOG_ERROR("No dispersion parameter provided for a period");
    }

    constexpr auto EXTINCTION_KEY = "extinction";
    if (rates[EXTINCTION_KEY]) {
      ext = rates[EXTINCTION_KEY].as<double>();
    } else {
      LOG_ERROR("No extinction parameter provided for a period");
    }
    return {{dis, ext}};
  }
  return {};
}

std::optional<bigrig::cladogenesis_params_t>
cli_options_t::get_cladogenesis(const YAML::Node &yaml) {
  double allopatry_rate = 0.0, sympatry_rate = 0.0, copy_rate = 0.0,
         jump_rate              = 0.0;
  bool           ok             = true;
  constexpr auto CLADO_DICT_KEY = "cladogenesis";
  if (yaml[CLADO_DICT_KEY]) {
    auto clado = yaml[CLADO_DICT_KEY];

    constexpr auto ALLOPATRY_KEY = "allopatry";
    if (clado[ALLOPATRY_KEY]) {
      allopatry_rate = clado[ALLOPATRY_KEY].as<double>();
    } else {
      LOG_ERROR("No allopatry parameter provided for a period");
      ok = false;
    }

    constexpr auto SYMPATRY_KEY = "sympatry";
    if (clado[SYMPATRY_KEY]) {
      sympatry_rate = clado[SYMPATRY_KEY].as<double>();
    } else {
      LOG_ERROR("No sympatry parameter provided for a period");
      ok = false;
    }

    constexpr auto COPY_KEY = "copy";
    if (clado[COPY_KEY]) {
      copy_rate = clado[COPY_KEY].as<double>();
    } else {
      LOG_ERROR("No copy parameter provided for a period");
      ok = false;
    }

    constexpr auto JUMP_KEY = "jump";
    if (clado[JUMP_KEY]) {
      jump_rate = clado[JUMP_KEY].as<double>();
    } else {
      LOG_ERROR("No jump parameter provided for a period");
      ok = false;
    }

    if (ok) { return {{allopatry_rate, sympatry_rate, copy_rate, jump_rate}}; }
  }
  return {};
}

std::optional<bigrig::tree_params_t> get_tree_params(const YAML::Node &yaml) {
  constexpr auto LAMBDA_KEY = "lambda";
  if (yaml[LAMBDA_KEY]) {
    return {
        bigrig::tree_params_t{.cladogenesis = yaml[LAMBDA_KEY].as<double>()}};
  }
  return {};
}

std::optional<bool> get_period_extinction(const YAML::Node &yaml) {
  constexpr auto EXTINCTION_KEY = "allow-extinction";
  return get_yaml_val_or_nothing<bool>(yaml, EXTINCTION_KEY);
}

std::optional<bigrig::per_region_params_t>
cli_options_t::get_per_region_params(const YAML::Node &yaml) {
  bigrig::per_region_params_t ret;
  constexpr auto              NAME_KEY  = "name";
  constexpr auto              DIST_KEY  = "dist";
  constexpr auto              INDEX_KEY = "index";

  if (yaml[NAME_KEY]) { ret.region_id = yaml[NAME_KEY].as<std::string>(); }
  if (yaml[DIST_KEY]) {
    ret.region_id = bigrig::dist_t{yaml[DIST_KEY].as<std::string>()};
  }
  if (yaml[INDEX_KEY]) { ret.region_id = yaml[INDEX_KEY].as<size_t>(); }

  ret.rates        = get_rates(yaml);
  ret.cladogenesis = get_cladogenesis(yaml);

  return ret;
}

std::optional<std::vector<bigrig::per_region_params_t>>
cli_options_t::get_per_region_params_list(const YAML::Node &yaml) {
  constexpr auto REGIONS_KEY = "regions";
  if (!yaml[REGIONS_KEY]) { return {}; }

  std::vector<bigrig::per_region_params_t> params;

  for (auto node : yaml[REGIONS_KEY]) {
    auto ret = get_per_region_params(node);
    if (ret) { params.push_back(*ret); }
  }

  return params;
}

std::optional<bigrig::period_params_t>
cli_options_t::get_period(const YAML::Node &yaml) {
  bigrig::period_params_t period_params;
  bool                    ok        = true;
  constexpr auto          START_KEY = "start";

  if (yaml[START_KEY]) {
    period_params.start = yaml[START_KEY].as<double>();
  } else {
    LOG_ERROR("No start time provided for a period");
    ok = false;
  }

  auto rates = get_rates(yaml);
  if (rates.has_value()) {
    period_params.rates = rates.value();
  } else {
    LOG_ERROR("Rates for a period are malformed");
    ok = false;
  }

  auto clado_params = get_cladogenesis(yaml);
  if (clado_params.has_value()) {
    period_params.clado = clado_params.value();
  } else {
    LOG_ERROR("Cladogenesis parameters for a period are malformed");
    ok = false;
  }

  period_params.tree       = get_tree_params(yaml);
  period_params.extinction = get_period_extinction(yaml);

  period_params.adjustment_matrix = get_adjustment_matrix_parameters(yaml);

  period_params.per_region_params
      = get_per_region_params_list(yaml).value_or({});

  if (ok) { return period_params; }
  return {};
}

std::vector<bigrig::period_params_t>
cli_options_t::get_periods(const YAML::Node &yaml) {
  std::vector<bigrig::period_params_t> ret;

  constexpr auto PERIODS_KEY = "periods";
  if (yaml[PERIODS_KEY]) {
    auto   list  = yaml[PERIODS_KEY];
    size_t index = 0;
    for (auto y : list) {
      auto period = get_period(y);
      if (!period.has_value()) {
        LOG_ERROR("Period {} is malformed", index);
        return {};
      }
      index++;
      ret.push_back(period.value());
    }
  } else {
    auto rates             = get_rates(yaml);
    auto clado_params      = get_cladogenesis(yaml);
    auto tree_params       = get_tree_params(yaml);
    auto extinction        = get_period_extinction(yaml);
    auto matrix            = get_adjustment_matrix_parameters(yaml);
    auto per_region_params = get_per_region_params_list(yaml);
    bool ok                = true;

    if (!per_region_params.has_value()) {
      if (!rates.has_value()) {
        ok = false;
        LOG_ERROR("Failed to find rates in the config file");
      }

      if (!clado_params.has_value()) {
        ok = false;
        LOG_ERROR("Failed to find cladogensis parameters in the config file");
      }
    }
    if (ok) {
      return {bigrig::period_params_t{
          .rates             = rates.value_or({}),
          .clado             = clado_params.value_or({}),
          .start             = 0.0,
          .tree              = tree_params,
          .per_region_params = per_region_params.value_or({}),
          .extinction        = extinction,
          .adjustment_matrix = {},
      }};
    }
  }
  return ret;
}

bigrig::period_params_t cli_options_t::default_period_params() {
  return {.rates = {.dis = 1.0, .ext = 1.0},
          .clado
          = {.allopatry = 1.0, .sympatry = 1.0, .copy = 1.0, .jump = 1.0},
          .start             = 0,
          .tree              = {},
          .per_region_params = {},
          .extinction        = {},
          .adjustment_matrix = {}};
}

std::optional<bool> cli_options_t::get_redo(const YAML::Node &yaml) {
  constexpr auto REDO_KEY = "redo";
  return get_yaml_val_or_nothing<bool>(yaml, REDO_KEY);
}

std::optional<bigrig::operation_mode_e>
cli_options_t::get_mode(const YAML::Node &yaml) {
  constexpr auto MODE_KEY = "mode";
  if (yaml[MODE_KEY]) {
    std::string value = yaml[MODE_KEY].as<std::string>();
    std::transform(value.begin(),
                   value.end(),
                   value.begin(),
                   [](char c) -> char { return std::tolower(c); });
    if (value == "fast") {
      return bigrig::operation_mode_e::FAST;
    } else if (value == "sim") {
      return bigrig::operation_mode_e::SIM;
    } else {
      throw cli_option_invalid_parameter{
          "Failed to recognize the run mode in the config file"};
    }
  }
  return {};
}

std::optional<uint64_t> cli_options_t::get_seed(const YAML::Node &yaml) {
  constexpr auto SEED_KEY = "seed";
  return get_yaml_val_or_nothing<uint64_t>(yaml, SEED_KEY);
}

std::optional<bool> cli_options_t::get_simulate_tree(const YAML::Node &yaml) {
  constexpr auto SIMULATE_TREE_KEY = "simulate-tree";
  return get_yaml_val_or_nothing<bool>(yaml, SIMULATE_TREE_KEY);
}

std::optional<double> cli_options_t::get_tree_height(const YAML::Node &yaml) {
  constexpr auto TREE_HEIGHT_KEY = "tree-height";
  return get_yaml_val_or_nothing<double>(yaml, TREE_HEIGHT_KEY);
}

std::optional<bigrig::adjustment_matrix_params_t>
cli_options_t::get_adjustment_matrix_parameters(const YAML::Node &yaml) {
  constexpr auto ADJUSTMENT_PARAMS_KEY = "adjust";
  if (!yaml[ADJUSTMENT_PARAMS_KEY]) { return {}; }

  auto sub_yaml = yaml[ADJUSTMENT_PARAMS_KEY];
  return {bigrig::adjustment_matrix_params_t{
      .matrix_filename = get_adjustment_matrix_filename(sub_yaml),
      .adjustments     = std::nullopt,
      .exponent        = get_distance_exponent(sub_yaml),
      .simulate        = get_simulate_adjustment_matrix(sub_yaml),
  }};
}

std::expected<std::filesystem::path, bigrig::io_err>
cli_options_t::get_adjustment_matrix_filename(const YAML::Node &yaml) {
  constexpr auto DISTANCE_MATRIX_FILENAME_KEY = "file";
  if (!yaml[DISTANCE_MATRIX_FILENAME_KEY]) {
    return std::unexpected{bigrig::io_err::KeyNotFound};
  }
  return {yaml[DISTANCE_MATRIX_FILENAME_KEY].as<std::string>()};
}

std::optional<double>
cli_options_t::get_distance_exponent(const YAML::Node &yaml) {
  constexpr auto DISTANCE_MATRIX_EXPONENT_KEY = "exponent";
  return get_yaml_val_or_nothing<double>(yaml, DISTANCE_MATRIX_EXPONENT_KEY);
}

std::optional<bool>
cli_options_t::get_simulate_adjustment_matrix(const YAML::Node &yaml) {
  constexpr auto SIMULATE_DISTANCE_MATRIX = "simulate";
  return get_yaml_val_or_nothing<bool>(yaml, SIMULATE_DISTANCE_MATRIX);
}

template <typename T>
[[nodiscard]] bool check_passed_cli_parameter(const std::optional<T> &o,
                                              const char             *name) {
  if (!o.has_value()) {
    LOG_ERROR("No value provided for parameter '{}' when other values were "
              "provided on the command line",
              name);
    return false;
  }
  return true;
}

[[nodiscard]] bool
cli_options_t::convert_cli_parameters(std::optional<double> dis,
                                      std::optional<double> ext,
                                      std::optional<double> allo,
                                      std::optional<double> symp,
                                      std::optional<double> copy,
                                      std::optional<double> jump) {
  if (!(dis.has_value() || ext.has_value() || allo.has_value()
        || symp.has_value() || copy.has_value() || jump.has_value())) {
    return true;
  }
  bool ok  = true;
  ok      &= check_passed_cli_parameter(dis, "dispersion");
  ok      &= check_passed_cli_parameter(ext, "extinction");
  ok      &= check_passed_cli_parameter(allo, "allopatry");
  ok      &= check_passed_cli_parameter(symp, "sympatry");
  ok      &= check_passed_cli_parameter(copy, "copy");
  ok      &= check_passed_cli_parameter(jump, "jump");

  if (ok) {
    bigrig::period_params_t period{};
    period.start = 0.0;
    period.rates = {.dis = dis.value(), .ext = ext.value()};
    period.clado = {.allopatry = allo.value(),
                    .sympatry  = ext.value(),
                    .copy      = copy.value(),
                    .jump      = jump.value()};
    periods.push_back(period);
    return true;
  }
  return false;
}
