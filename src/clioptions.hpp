#pragma once
#include "dist.hpp"

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <logger.hpp>
#include <nlohmann/json.hpp>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <yaml-cpp/exceptions.h>
#include <yaml-cpp/yaml.h>

constexpr auto PHYILP_EXT = ".phy";
constexpr auto NEWICK_EXT = ".nwk";
constexpr auto YAML_EXT   = ".yaml";
constexpr auto JSON_EXT   = ".json";

enum class output_format_type_e { JSON, YAML };

class cli_option_missing_required_yaml_option : std::invalid_argument {
public:
  cli_option_missing_required_yaml_option(const std::string &msg)
      : std::invalid_argument{msg} {}
};

struct cli_options_t {
  std::optional<std::filesystem::path> config_filename;
  std::optional<std::filesystem::path> tree_filename;
  std::optional<std::filesystem::path> prefix;
  std::optional<bool>                  debug_log;
  std::optional<output_format_type_e>  output_format_type;
  std::optional<biogeosim::dist_t>     root_distribution;
  std::optional<double>                dispersion_rate;
  std::optional<double>                extinction_rate;
  std::optional<bool>                  redo;
  std::optional<bool>                  two_region_duplicity;

  std::filesystem::path phylip_filename() const {
    auto tmp  = prefix.value();
    tmp      += PHYILP_EXT;
    return tmp;
  }

  std::filesystem::path yaml_filename() const {
    auto tmp  = prefix.value();
    tmp      += YAML_EXT;
    return tmp;
  }

  std::filesystem::path json_filename() const {
    auto tmp  = prefix.value();
    tmp      += JSON_EXT;
    return tmp;
  }

  bool cli_arg_specified() const {
    return tree_filename.has_value() || prefix.has_value()
        || debug_log.has_value() || tree_filename.has_value();
  }

  bool yaml_file_set() const {
    return output_format_type.has_value()
        && output_format_type.value() == output_format_type_e::YAML;
  }

  bool json_file_set() const {
    return output_format_type.has_value()
        && output_format_type.value() == output_format_type_e::JSON;
  }

  void merge(const cli_options_t &other) {
    if (other.config_filename.has_value()) {
      config_filename = other.config_filename;
    }
    if (other.tree_filename.has_value()) {
      tree_filename = other.tree_filename;
    }
    if (other.debug_log.has_value()) { debug_log = other.debug_log; }
    if (other.output_format_type.has_value()) {
      output_format_type = other.output_format_type;
    }
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

  cli_options_t() = default;
  cli_options_t(const YAML::Node &yaml) {
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
    if (yaml["dispersion"]) {
      dispersion_rate = yaml["dispersion"].as<double>();
    }
    if (yaml["extinction"]) {
      extinction_rate = yaml["extinction"].as<double>();
    }
  }
};
