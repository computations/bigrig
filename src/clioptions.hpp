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

  std::filesystem::path phylip_filename() const {
    auto tmp  = prefix.value();
    tmp      += PHYILP_EXT;
    return tmp;
  }

  bool cli_arg_specified() const {
    return tree_filename.has_value() || prefix.has_value()
        || debug_log.has_value() || tree_filename.has_value();
  }

  cli_options_t() = default;
  cli_options_t(const YAML::Node &yaml) {
    tree_filename = yaml["tree"].as<std::string>();
    if (yaml["prefix"]) { prefix = yaml["prefix"].as<std::string>(); }
    if (yaml["debug-log"]) { debug_log = yaml["debug-log"].as<bool>(); }
    if (yaml["output-format"]
        && yaml["output-format"].as<std::string>() == "json") {
      output_format_type = output_format_type_e::JSON;
    }
    if (yaml["output-format"]
        && yaml["output-format"].as<std::string>() == "yaml") {
      output_format_type = output_format_type_e::YAML;
    }
    root_distribution = yaml["root-dist"].as<std::string>();
    dispersion_rate   = yaml["dispersion"].as<double>();
    extinction_rate   = yaml["extinction"].as<double>();
  }
};
