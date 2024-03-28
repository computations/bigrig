#include "clioptions.hpp"

#include "dist.hpp"
#include "logger.hpp"
#include "util.hpp"

#include <algorithm>
#include <tuple>

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
  LOG_WARNING("The '%s' option is specified in both the config file and "
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
  merge_variable(debug_log, other.debug_log, "debug-log");
  merge_variable(output_format_type, other.output_format_type, "output-format");
  merge_variable(root_distribution, other.root_distribution, "root-range");
  merge_variable(dispersion_rate, other.dispersion_rate, "rates:dispersion");
  merge_variable(extinction_rate, other.extinction_rate, "rates:extinction");
  merge_variable(
      allopatry_rate, other.allopatry_rate, "cladogenesis:allopatry");
  merge_variable(sympatry_rate, other.sympatry_rate, "cladogenesis:sympatry");
  merge_variable(copy_rate, other.copy_rate, "cladogenesis:copy");
  merge_variable(jump_rate, other.jump_rate, "cladogenesis:jump");
  merge_variable(redo, other.redo, "redo");
  merge_variable(
      two_region_duplicity, other.two_region_duplicity, "two-region-duplicity");
  merge_variable(mode, other.mode, "mode");
  merge_variable(rng_seed, other.rng_seed, "seed");
}

std::filesystem::path cli_options_t::get_tree_filename(const YAML::Node &yaml) {
  constexpr auto TREE_KEY = "tree";
  return yaml[TREE_KEY].as<std::string>();
}

std::optional<std::filesystem::path>
cli_options_t::get_prefix(const YAML::Node &yaml) {
  constexpr auto PREFIX_KEY = "prefix";
  if (yaml[PREFIX_KEY]) { return yaml[PREFIX_KEY].as<std::string>(); }
  return {};
}

std::optional<bool> cli_options_t::get_debug_log(const YAML::Node &yaml) {
  constexpr auto DEBUG_LOG_KEY = "debug-log";
  if (yaml[DEBUG_LOG_KEY]) { return yaml[DEBUG_LOG_KEY].as<bool>(); }
  return {};
}

std::optional<output_format_type_e>
cli_options_t::get_output_format(const YAML::Node &yaml) {
  constexpr auto OUPUT_FORMAT_KEY = "output-format";
  if (yaml[OUPUT_FORMAT_KEY]) {
    if (yaml[OUPUT_FORMAT_KEY].as<std::string>() == "json") {
      return output_format_type_e::JSON;
    }
    if (yaml[OUPUT_FORMAT_KEY].as<std::string>() == "yaml") {
      return output_format_type_e::YAML;
    }
  }
  return {};
}

std::optional<bool>
cli_options_t::get_two_region_duplicity(const YAML::Node &yaml) {
  constexpr auto DUPLICITY_KEY = "two-region-duplicity";
  if (yaml[DUPLICITY_KEY]) { return yaml[DUPLICITY_KEY].as<bool>(); }
  return {};
}

std::optional<bigrig::dist_t>
cli_options_t::get_root_range(const YAML::Node &yaml) {
  constexpr auto ROOT_DIST_KEY = "root-range";
  if (yaml[ROOT_DIST_KEY]) { return yaml[ROOT_DIST_KEY].as<std::string>(); }
  return {};
}

std::tuple<std::optional<double>, std::optional<double>>
cli_options_t::get_rates(const YAML::Node &yaml) {
  std::optional<double> dis, ext;
  constexpr auto        RATE_DICT_KEY = "rates";
  if (yaml[RATE_DICT_KEY]) {
    auto rates = yaml[RATE_DICT_KEY];

    constexpr auto DISPERSION_KEY = "dispersion";
    if (rates[DISPERSION_KEY]) { dis = rates[DISPERSION_KEY].as<double>(); }

    constexpr auto EXTINCTION_KEY = "extinction";
    if (rates[EXTINCTION_KEY]) { ext = rates[EXTINCTION_KEY].as<double>(); }
  }
  return {dis, ext};
}

std::tuple<std::optional<double>,
           std::optional<double>,
           std::optional<double>,
           std::optional<double>>
cli_options_t::get_cladogenesis(const YAML::Node &yaml) {
  std::optional<double> allopatry_rate, sympatry_rate, copy_rate, jump_rate;
  constexpr auto        CLADO_DICT_KEY = "cladogenesis";
  if (yaml[CLADO_DICT_KEY]) {
    auto clado = yaml[CLADO_DICT_KEY];

    constexpr auto ALLOPATRY_KEY = "allopatry";
    if (clado[ALLOPATRY_KEY]) {
      allopatry_rate = clado[ALLOPATRY_KEY].as<double>();
    }

    constexpr auto SYMPATRY_KEY = "sympatry";
    if (clado[SYMPATRY_KEY]) {
      sympatry_rate = clado[SYMPATRY_KEY].as<double>();
    }

    constexpr auto COPY_KEY = "copy";
    if (clado[COPY_KEY]) { copy_rate = clado[COPY_KEY].as<double>(); }

    constexpr auto JUMP_KEY = "jump";
    if (clado[JUMP_KEY]) { jump_rate = clado[JUMP_KEY].as<double>(); }
  }
  return {allopatry_rate, sympatry_rate, copy_rate, jump_rate};
}

std::optional<bool> cli_options_t::get_redo(const YAML::Node &yaml) {
  constexpr auto REDO_KEY = "redo";
  if (yaml[REDO_KEY]) { return yaml[REDO_KEY].as<bool>(); }
  return {};
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
  if (yaml[SEED_KEY]) { return yaml[SEED_KEY].as<uint64_t>(); }
  return {};
}
