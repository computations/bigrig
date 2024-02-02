#include "clioptions.hpp"

std::filesystem::path cli_options_t::phylip_filename() const {
  auto tmp  = prefix.value();
  tmp      += biogeosim::util::PHYILP_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::yaml_filename() const {
  auto tmp  = prefix.value();
  tmp      += biogeosim::util::YAML_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::json_filename() const {
  auto tmp  = prefix.value();
  tmp      += biogeosim::util::JSON_EXT;
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
void cli_options_t::merge(const cli_options_t &other) {
  if (other.config_filename.has_value()) {
    config_filename = other.config_filename;
  }
  if (other.tree_filename.has_value()) { tree_filename = other.tree_filename; }
  if (other.debug_log.has_value()) { debug_log = other.debug_log; }
  if (other.output_format_type.has_value()) {
    output_format_type = other.output_format_type;
  }
  if (other.root_distribution.has_value()) {
    root_distribution = other.root_distribution;
  }
  if (other.dispersion_rate.has_value()) {
    dispersion_rate = other.dispersion_rate;
  }
  if (other.extinction_rate.has_value()) {
    extinction_rate = other.extinction_rate;
  }
  if (other.redo.has_value()) { redo = other.redo; }
  if (other.two_region_duplicity.has_value()) {
    two_region_duplicity = other.two_region_duplicity;
  }
}

/**
 * Construct a cli_option_t from a YAML file. Please see the YAML file schema
 * for the details.
 */
cli_options_t::cli_options_t(const YAML::Node &yaml) {
  tree_filename = yaml["tree"].as<std::string>();
  if (yaml["prefix"]) { prefix = yaml["prefix"].as<std::string>(); }
  if (yaml["debug-log"]) { debug_log = yaml["debug-log"].as<bool>(); }
  if (yaml["output-format"]) {
    if (yaml["output-format"].as<std::string>() == "json") {
      output_format_type = output_format_type_e::JSON;
    }
    if (yaml["output-format"].as<std::string>() == "yaml") {
      output_format_type = output_format_type_e::YAML;
    }
  }
  if (output_format_type.has_value()) {
    MESSAGE_WARNING(
        "Output format specified in both the config file and on the command "
        "line. Using the specification from the config file.");
  }
  if (yaml["two-region-duplicity"]) {
    two_region_duplicity = yaml["two-region-duplicity"].as<bool>();
  }
  if (yaml["root-dist"]) {
    root_distribution = yaml["root-dist"].as<std::string>();
  }
  if (yaml["dispersion"]) { dispersion_rate = yaml["dispersion"].as<double>(); }
  if (yaml["extinction"]) { extinction_rate = yaml["extinction"].as<double>(); }
}
