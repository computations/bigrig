#include "model.hpp"
#include "node.hpp"
#include "pcg_random.hpp"
#include "tree.hpp"
#include "yaml-cpp/exceptions.h"

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <corax/corax.hpp>
#include <filesystem>
#include <fstream>
#include <logger.hpp>
#include <optional>
#include <sstream>
#include <string>
#include <yaml-cpp/yaml.h>

constexpr auto PHYILP_EXT = ".phy";

struct cli_options_t {
  std::optional<std::filesystem::path> config_filename;
  std::optional<std::filesystem::path> tree_filename;
  std::optional<std::filesystem::path> prefix;
  std::optional<bool>                  debug_log;
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
    root_distribution = yaml["root-dist"].as<std::string>();
    dispersion_rate   = yaml["dispersion"].as<double>();
    extinction_rate   = yaml["extinction"].as<double>();
  }
};

std::string to_phylip(const biogeosim::tree_t              &tree,
                      const biogeosim::substitution_model_t model) {
  std::ostringstream oss;
  oss << std::to_string(tree.leaf_count()) << " " << model.region_count()
      << "\n";

  tree.to_phylip_body(oss);

  return oss.str();
}

std::string to_phylip_extended(const biogeosim::tree_t              &tree,
                               const biogeosim::substitution_model_t model) {
  std::ostringstream oss;
  oss << std::to_string(tree.node_count()) << " " << model.region_count()
      << "\n";

  tree.to_phylip_body(oss, true);

  return oss.str();
}

[[nodiscard]] bool verify_is_readable(const std::filesystem::path &path) {
  bool ok     = true;
  auto status = std::filesystem::status(path);
  if ((std::filesystem::perms::owner_read & status.permissions())
      == std::filesystem::perms::none) {
    LOG_ERROR("The path '%s' is not readable", path.c_str());
    ok = false;
  }
  return ok;
}

[[nodiscard]] bool verify_is_writable(const std::filesystem::path &path) {
  bool ok             = true;
  auto status         = std::filesystem::status(path);
  auto required_perms = std::filesystem::perms::owner_write
                      | std::filesystem::perms::owner_exec;
  if ((required_perms & status.permissions()) == std::filesystem::perms::none) {
    LOG_ERROR("The path '%s' is not writable", path.c_str());
    ok = false;
  }
  return ok;
}

[[nodiscard]] bool
validate_tree_filename(const std::filesystem::path &tree_filename) {
  bool ok = true;
  if (!std::filesystem::exists(tree_filename)) {
    LOG_ERROR("The tree file '%s' does not exist", tree_filename.c_str());
    ok = false;
  } else if (!std::filesystem::is_regular_file(tree_filename)) {
    LOG_ERROR("The tree file '%s' is not a file that we can read",
              tree_filename.c_str());
    ok = false;
  }
  if (!verify_is_readable(tree_filename)) {
    LOG_ERROR("The tree file '%s' can't be read by us as we don't have the "
              "permissions",
              tree_filename.c_str());
    ok = false;
  }
  return ok;
}

[[nodiscard]] bool validate_prefix(const std::filesystem::path &prefix) {
  bool ok = true;
  if (!std::filesystem::exists(prefix.parent_path())) {
    LOG_WARNING("The path '%s' does not exist", prefix.parent_path().c_str());

    try {
      std::filesystem::create_directories(prefix.parent_path());
    } catch (std::filesystem::filesystem_error err) {
      LOG_ERROR("%s", err.what());
      ok = false;
    }
  } else if (!verify_is_writable(prefix.parent_path())) {
    LOG_ERROR("The prefix '%s' is not writable", prefix.c_str());
    ok = false;
  }

  return ok;
}

[[nodiscard]] bool
verify_config_file(const std::filesystem::path &config_filename) {
  bool ok = true;
  if (!std::filesystem::exists(config_filename)) {
    LOG_ERROR("The config file %s does not exist", config_filename.c_str());
    ok = false;
  } else if (!verify_is_readable(config_filename)) {
    LOG_ERROR("We dont' have the permissions to read the config file %s",
              config_filename.c_str());
    ok = false;
  }

  return ok;
}

[[nodiscard]] bool validate_cli_options(const cli_options_t &cli_options) {
  bool ok = true;

  ok &= validate_tree_filename(cli_options.tree_filename.value());
  ok &= validate_prefix(cli_options.prefix.value());

  if (cli_options.root_distribution.regions() > 64) {
    LOG_ERROR("Simulating with %u regions is unsupported. Please choose a "
              "number less than or equal to 64",
              cli_options.root_distribution.regions());
    ok = false;
  }
  if (cli_options.dispersion_rate < 0) {
    LOG_ERROR("Simulating with dispersion rate = %f is not valid, please pick "
              "a positive number",
              cli_options.dispersion_rate);
    ok = false;
  }
  if (cli_options.extinction_rate < 0) {
    LOG_ERROR("Simulating with dispersion rate = %f is not valid, please pick "
              "a positive number",
              cli_options.extinction_rate);
    ok = false;
  }

  return ok;
}

[[nodiscard]] bool config_or_cli(const cli_options_t &cli_options) {
  bool ok = true;

  if (cli_options.config_filename.has_value()
      && cli_options.cli_arg_specified()) {
    MESSAGE_ERROR("Both config and CLI arguements have been specified. Please "
                  "only specify one");
    ok = false;
  }
  if (cli_options.config_filename.has_value()) {
    auto &config_filename = cli_options.config_filename.value();
    if (!verify_config_file(config_filename)) {
      LOG_ERROR("There was an issue with the config file %s",
                config_filename.c_str());
      return false;
    }
  }

  return ok;
}

void write_header(const cli_options_t &cli_options) {
  MESSAGE_INFO("Running simulation with the following parameters:");
  LOG_INFO("   Tree file: %s", cli_options.tree_filename.value().c_str());
  LOG_INFO("   Prefix: %s", cli_options.prefix.value().c_str());
  LOG_INFO("   Root Distribution: %s",
           cli_options.root_distribution.to_str().c_str());
}

bool normalize_paths(cli_options_t &cli_options) {
  bool ok = true;
  // std::filesystem::canonical will throw an error here, and we might want to
  // make the paths later on. So we instead call absolute and then
  // weakly canonicalize them
  try {
    cli_options.tree_filename = std::filesystem::weakly_canonical(
        std::filesystem::absolute(cli_options.tree_filename.value()));
  } catch (const std::filesystem::filesystem_error &err) {
    LOG_ERROR("Failed to canonicalize '%s' because '%s'",
              cli_options.tree_filename.value().c_str(),
              err.what());
    ok = false;
  }
  try {
    cli_options.prefix = std::filesystem::weakly_canonical(
        std::filesystem::absolute(cli_options.prefix.value()));
  } catch (const std::filesystem::filesystem_error &err) {
    LOG_ERROR("Failed to canonicalize '%s' because '%s'",
              cli_options.prefix.value().c_str(),
              err.what());
    ok = false;
  }
  return ok;
}

[[nodiscard]] bool check_existing_results(const cli_options_t &cli_options) {
  bool ok = true;

  if (std::filesystem::exists(cli_options.phylip_filename())) {
    LOG_WARNING("Results file %s exists already",
                cli_options.phylip_filename().c_str());
    ok = false;
  }
  return ok;
}

cli_options_t parse_yaml_options(const std::filesystem::path &config_filename) {
  auto          yaml = YAML::LoadFile(config_filename);
  cli_options_t cli_options(yaml);

  return cli_options;
}

int main(int argc, char **argv) {
  logger::get_log_states().add_stream(
      stdout,
      logger::log_level::info | logger::log_level::warning
          | logger::log_level::important | logger::log_level::error
          | logger::log_level::progress);

  CLI::App app{"Biogeosim"};

  cli_options_t cli_options;

  app.add_option(
      "--config", cli_options.config_filename, "Config file with options");
  app.add_option("--tree",
                 cli_options.tree_filename,
                 "A file containing a newick encoded tree which will be used "
                 "to perform the simulation");
  app.add_option("--prefix", cli_options.prefix, "prefix for the output files");
  app.add_option("--root-dist",
                 cli_options.root_distribution,
                 "Distribution at the root at the start of the simulation");
  app.add_option("--d", cli_options.dispersion_rate, "Dispersion rate");
  app.add_option("--e", cli_options.extinction_rate, "Extinction rate");
  app.add_flag("--redo", cli_options.redo, "Extinction rate")
      ->default_val(false);
  app.add_flag("--debug-log",
               cli_options.debug_log,
               "Create a file in the prefix that contains the debug log");

  CLI11_PARSE(app);

  if (cli_options.config_filename.has_value() && config_or_cli(cli_options)) {
    try {
      auto cli_options_tmp
          = parse_yaml_options(cli_options.config_filename.value());
      cli_options_tmp.redo            = cli_options.redo;
      cli_options_tmp.config_filename = cli_options.config_filename;
      std::swap(cli_options, cli_options_tmp);
    } catch (const YAML::Exception &e) {
      LOG_ERROR("Failed to parse the config file: %s", e.what());
      return 1;
    }

    LOG_INFO("tree: %s", cli_options.tree_filename.value().c_str());
  }

  if (!cli_options.prefix.has_value() || cli_options.prefix.value().empty()) {
    cli_options.prefix = cli_options.tree_filename;
  }

  normalize_paths(cli_options);
  write_header(cli_options);

  if (!validate_cli_options(cli_options)) {
    MESSAGE_ERROR(
        "We can't continue with the current options, exiting instead");
    return 1;
  }

  if (!check_existing_results(cli_options)) {
    if (!cli_options.redo) {
      MESSAGE_ERROR("Refusing to run with existing results. Please specify the "
                    "--redo option if you want to overwrite existig results");
      return 1;
    }
  }

  if (cli_options.debug_log) {
    std::filesystem::path debug_filename  = cli_options.prefix.value();
    debug_filename                       += ".debug.log";
    LOG_INFO("Logging debug information to %s", debug_filename.c_str());
    logger::get_log_states().add_file_stream(
        debug_filename.c_str(),
        logger::log_level::info | logger::log_level::warning
            | logger::log_level::important | logger::log_level::error
            | logger::log_level::debug);
  }

  auto tree = biogeosim::tree_t(cli_options.tree_filename.value().c_str());

  pcg_extras::seed_seq_from<std::random_device> seed_source;
  pcg64_fast                                    gen(seed_source);

  biogeosim::substitution_model_t model(
      cli_options.dispersion_rate,
      cli_options.extinction_rate,
      cli_options.root_distribution.regions());

  tree.sample(cli_options.root_distribution, model, gen);

  auto phylip_filename  = cli_options.prefix.value();
  phylip_filename      += ".phy";
  std::ofstream phylip_file(phylip_filename);
  phylip_file << to_phylip(tree, model);

  auto phylip_all_filename  = cli_options.prefix.value();
  phylip_all_filename      += ".all.phy";
  std::ofstream phylip_all_file(phylip_all_filename);
  phylip_all_file << to_phylip_extended(tree, model);

  auto cb = [](std::ostream &os, biogeosim::node_t n) {
    if (n.is_leaf()) {
      os << n.label();
    } else {
      os << n.node_id();
      os << "[&&NHX:";
      os << n.node_split().to_nhx_string();
      os << "]";
    }
  };

  auto annotated_tree_filename  = cli_options.prefix.value();
  annotated_tree_filename      += ".annotated.nwk";
  std::ofstream annotated_tree_file(annotated_tree_filename);
  annotated_tree_file << tree.to_newick(cb) << std::endl;

  return 0;
}
