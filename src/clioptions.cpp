#include "clioptions.hpp"

#include "dist.hpp"
#include "logger.hpp"
#include "util.hpp"

#include <algorithm>

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
  if (other.tree_filename.has_value()) {
    if (tree_filename.has_value()) {
      MESSAGE_WARNING(
          "A tree file is specified in both the config file and "
          "the command line. Using the value from the command line");
    } else {
      tree_filename = other.tree_filename;
    }
  }

  if (other.debug_log.has_value()) {
    if (debug_log.has_value()) {
      MESSAGE_WARNING(
          "The debug log option is specified in both the config file and the "
          "command line. Using the value from the command line");
    } else {
      debug_log = other.debug_log;
    }
  }

  if (other.output_format_type.has_value()) {
    if (output_format_type.has_value()) {
      MESSAGE_WARNING(
          "The output format is specified in both the config file and "
          "the command line. Using the value from the command line");
    } else {
      output_format_type = other.output_format_type;
    }
  }

  if (other.root_distribution.has_value()) {
    if (root_distribution.has_value()) {
      MESSAGE_WARNING(
          "The root range is specified in both the config file and "
          "the command line. Using the value from the command line");
    } else {
      root_distribution = other.root_distribution;
    }
  }

  if (other.dispersion_rate.has_value()) {
    if (dispersion_rate.has_value()) {
      MESSAGE_WARNING(
          "The dispersion rate is specified in both the config file and "
          "the command line. Using the value from the command line");
    } else {
      dispersion_rate = other.dispersion_rate;
    }
  }

  if (other.extinction_rate.has_value()) {
    if (extinction_rate.has_value()) {
      MESSAGE_WARNING(
          "The extinction rate is specified in both the config file and "
          "the command line. Using the value from the command line");
    } else {
      extinction_rate = other.extinction_rate;
    }
  }

  if (other.allopatry_rate.has_value()) {
    if (allopatry_rate.has_value()) {
      MESSAGE_WARNING(
          "The allopatry rate is specified in both the config file and "
          "the command line. Using the value from the command line");
    } else {
      allopatry_rate = other.allopatry_rate;
    }
  }

  if (other.sympatry_rate.has_value()) {
    if (sympatry_rate.has_value()) {
      MESSAGE_WARNING(
          "The sympatry rate is specified in both the config file and "
          "the command line. Using the value from the command line");
    } else {
      sympatry_rate = other.sympatry_rate;
    }
  }

  if (other.copy_rate.has_value()) {
    if (copy_rate.has_value()) {
      MESSAGE_WARNING(
          "The copy rate is specified in both the config file and "
          "the command line. Using the value from the command line");
    } else {
      copy_rate = other.copy_rate;
    }
  }

  if (other.jump_rate.has_value()) {
    if (jump_rate.has_value()) {
      MESSAGE_WARNING(
          "The jump rate is specified in both the config file and "
          "the command line. Using the value from the command line");
    } else {
      jump_rate = other.jump_rate;
    }
  }

  if (other.redo.has_value()) {
    if (redo.has_value()) {
      MESSAGE_WARNING(
          "Redo is specified in both the config file and the command "
          "line. Using the value from the command line");
    } else {
      redo = other.redo;
    }
  }

  if (other.two_region_duplicity.has_value()) {
    if (two_region_duplicity.has_value()) {
      MESSAGE_WARNING(
          "The option two value duplicity is specified in both the config file "
          "and the command line. Using the value from the command line");
    } else {
      two_region_duplicity = other.two_region_duplicity;
    }
  }

  if (other.mode.has_value()) {
    if (mode.has_value()) {
      MESSAGE_WARNING(
          "The option mode is specified in both the config file "
          "and the command line. Using the value from the command line");
    } else {
      mode = other.mode;
    }
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

  constexpr auto ROOT_DIST_KEY = "root-range";
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

  constexpr auto REDO_KEY = "redo";
  if (yaml[REDO_KEY]) { redo = yaml[REDO_KEY].as<bool>(); }

  constexpr auto MODE_KEY = "mode";
  if (yaml[MODE_KEY]) {
    std::string value = yaml[MODE_KEY].as<std::string>();
    std::transform(value.begin(),
                   value.end(),
                   value.begin(),
                   [](char c) -> char { return std::tolower(c); });
    if (value == "fast") {
      mode = bigrig::operation_mode_e::FAST;
    } else if (value == "sim") {
      mode = bigrig::operation_mode_e::SIM;
    } else {
      throw cli_option_invalid_parameter{
          "Failed to recognize the run mode in the config file"};
    }
  }
}
