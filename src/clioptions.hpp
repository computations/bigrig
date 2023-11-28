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
#include <string>
#include <yaml-cpp/exceptions.h>
#include <yaml-cpp/yaml.h>

constexpr auto PHYILP_EXT = ".phy";
constexpr auto NEWICK_EXT = ".nwk";
constexpr auto YAML_EXT   = ".yaml";
constexpr auto JSON_EXT   = ".json";

enum class output_format_type_e { JSON, YAML };

struct cli_options_t {
  std::optional<std::filesystem::path> config_filename;
  std::optional<std::filesystem::path> tree_filename;
  std::optional<std::filesystem::path> prefix;
  std::optional<bool>                  debug_log;
  std::optional<output_format_type_e>  output_format_type;
  biogeosim::dist_t                    root_distribution;
  double                               dispersion_rate;
  double                               extinction_rate;
  bool                                 redo;
  bool                                 two_region_duplicity;

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
    root_distribution = yaml["root-dist"].as<std::string>();
    dispersion_rate   = yaml["dispersion"].as<double>();
    extinction_rate   = yaml["extinction"].as<double>();
  }
};
