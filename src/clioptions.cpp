#include "clioptions.hpp"

#include "logger.hpp"

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

  if (other.dispersion_rate.has_value() && !dispersion_rate.has_value()) {
    dispersion_rate = other.dispersion_rate;
  }
  if (other.extinction_rate.has_value() && !extinction_rate.has_value()) {
    extinction_rate = other.extinction_rate;
  }

  if (other.allopatry_rate.has_value() && !allopatry_rate.has_value()) {
    allopatry_rate = other.allopatry_rate;
  }
  if (other.sympatry_rate.has_value() && !sympatry_rate.has_value()) {
    sympatry_rate = other.sympatry_rate;
  }
  if (other.copy_rate.has_value() && !copy_rate.has_value()) {
    copy_rate = other.copy_rate;
  }
  if (other.jump_rate.has_value() && !jump_rate.has_value()) {
    jump_rate = other.jump_rate;
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
  constexpr auto TREE_KEY = "tree";
  tree_filename           = yaml[TREE_KEY].as<std::string>();

  constexpr auto PREFIX_KEY = "prefix";
  if (yaml[PREFIX_KEY]) { prefix = yaml[PREFIX_KEY].as<std::string>(); }

  constexpr auto DEBUG_LOG_KEY = "debug-log";
  if (yaml[DEBUG_LOG_KEY]) { debug_log = yaml[DEBUG_LOG_KEY].as<bool>(); }

  constexpr auto OUPUT_FORMAT_KEY = "output-format";
  if (yaml[OUPUT_FORMAT_KEY]) {
    if (yaml[OUPUT_FORMAT_KEY].as<std::string>() == "json") {
      output_format_type = output_format_type_e::JSON;
    }
    if (yaml[OUPUT_FORMAT_KEY].as<std::string>() == "yaml") {
      output_format_type = output_format_type_e::YAML;
    }
  }
  if (output_format_type.has_value()) {
    MESSAGE_WARNING(
        "Output format specified in both the config file and on the command "
        "line. Using the specification from the config file.");
  }

  constexpr auto DUPLICITY_KEY = "two-region-duplicity";
  if (yaml[DUPLICITY_KEY]) {
    two_region_duplicity = yaml[DUPLICITY_KEY].as<bool>();
  }

  constexpr auto ROOT_DIST_KEY = "root-dist";
  if (yaml[ROOT_DIST_KEY]) {
    root_distribution = yaml[ROOT_DIST_KEY].as<std::string>();
  }

  constexpr auto RATE_DICT_KEY = "rates";
  if (yaml[RATE_DICT_KEY]) {
    auto rates = yaml[RATE_DICT_KEY];

    constexpr auto DISPERSION_KEY = "dispersion";
    if (rates[DISPERSION_KEY]) {
      dispersion_rate = rates[DISPERSION_KEY].as<double>();
    }

    constexpr auto EXTINCTION_KEY = "extinction";
    if (rates[EXTINCTION_KEY]) {
      extinction_rate = rates[EXTINCTION_KEY].as<double>();
    }
  }

  constexpr auto CLADO_DICT_KEY = "cladogenesis";
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
}
