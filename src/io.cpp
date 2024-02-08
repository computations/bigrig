#include "io.hpp"

#include "clioptions.hpp"
#include "logger.hpp"

#include <filesystem>
#include <format>
#include <optional>

/**
 * Print the message at the start of the run.
 */
void write_header(const cli_options_t &cli_options) {
  MESSAGE_INFO("Running simulation with the following parameters:");
  LOG_INFO("   Tree file: %s", cli_options.tree_filename.value().c_str());
  LOG_INFO("   Prefix: %s", cli_options.prefix.value().c_str());
  LOG_INFO("   Root Distribution: %s",
           cli_options.root_distribution.value().to_str().c_str());
}

/**
 * Produce a phylip file as a string.
 */
std::string to_phylip(const bigrig::tree_t        &tree,
                      const bigrig::biogeo_model_t model) {
  std::ostringstream oss;
  oss << std::to_string(tree.leaf_count()) << " " << model.region_count()
      << "\n";

  tree.to_phylip_body(oss);

  return oss.str();
}

/**
 * Produce a phylip file as a string, including inner nodes.
 */
std::string to_phylip_all_nodes(const bigrig::tree_t        &tree,
                                const bigrig::biogeo_model_t model) {
  std::ostringstream oss;
  oss << std::to_string(tree.node_count()) << " " << model.region_count()
      << "\n";

  tree.to_phylip_body(oss, true);

  return oss.str();
}

[[nodiscard]] bool verify_path_is_readable(const std::filesystem::path &path) {
  bool ok     = true;
  auto status = std::filesystem::status(path);
  if ((std::filesystem::perms::owner_read & status.permissions())
      == std::filesystem::perms::none) {
    LOG_ERROR("The path '%s' is not readable", path.c_str());
    ok = false;
  }
  return ok;
}

/**
 * Check that the tree file path is ok to use.
 */
[[nodiscard]] bool validate_tree_filename(
    const std::optional<std::filesystem::path> &tree_filename_option) {
  if (!tree_filename_option.has_value()) {
    MESSAGE_ERROR("No tree file was provided");
    return false;
  }
  const auto &tree_filename = tree_filename_option.value();
  bool        ok            = true;
  if (!std::filesystem::exists(tree_filename)) {
    LOG_ERROR("The tree file '%s' does not exist", tree_filename.c_str());
    ok = false;
  } else if (!std::filesystem::is_regular_file(tree_filename)) {
    LOG_ERROR("The tree file '%s' is not a file that we can read",
              tree_filename.c_str());
    ok = false;
  }
  if (!verify_path_is_readable(tree_filename)) {
    LOG_ERROR("The tree file '%s' can't be read by us as we don't have the "
              "permissions",
              tree_filename.c_str());
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
    LOG_ERROR("The path '%s' is not writable", path.c_str());
    ok = false;
  }
  return ok;
}

/**
 * Checks that the prefix is alright, and then make the directories.
 *
 * When the user passes a prefix, we need to check that its all kosher, and then
 * _actually_ make the directories. We make them here.
 */
[[nodiscard]] bool validate_and_make_prefix(
    const std::optional<std::filesystem::path> &prefix_option) {
  if (!prefix_option.has_value()) {
    MESSAGE_ERROR("No prefix was provided");
    return false;
  }
  auto prefix = prefix_option.value();
  bool ok     = true;
  if (!std::filesystem::exists(prefix.parent_path())) {
    LOG_WARNING("The path '%s' does not exist", prefix.parent_path().c_str());

    try {
      std::filesystem::create_directories(prefix.parent_path());
    } catch (const std::filesystem::filesystem_error &err) {
      LOG_ERROR("%s", err.what());
      ok = false;
    }
  } else if (!verify_path_is_writable(prefix.parent_path())) {
    LOG_ERROR("The prefix '%s' is not writable", prefix.c_str());
    ok = false;
  }

  return ok;
}

/**
 * Check that the program options are valid
 *
 * Check that the program options work. Note, it tries to be thorough when
 * checking, as the loop of "change a thing, find something else wrong" is
 * annoying. So, this function tries to check as much as it can, and not to bail
 * out at the first error.
 */
[[nodiscard]] bool validate_cli_options(const cli_options_t &cli_options) {
  bool ok = true;

  ok &= validate_tree_filename(cli_options.tree_filename);
  ok &= validate_and_make_prefix(cli_options.prefix);

  if (!cli_options.root_distribution.has_value()) {
    MESSAGE_ERROR("The root distribution was not provided. Please provide a "
                  "value for the root distribution");
    ok = false;
  } else if (cli_options.root_distribution.value().regions() > 64) {
    LOG_ERROR("Simulating with %u regions is unsupported. Please choose a "
              "number less than or equal to 64",
              cli_options.root_distribution.value().regions());
    ok = false;
  }
  if (!cli_options.dispersion_rate.has_value()) {
    MESSAGE_ERROR("Dispersion rate was not set. Please provide a value for "
                  "the dispersion rate");
    ok = false;
  } else if (cli_options.dispersion_rate.value() < 0) {
    LOG_ERROR("Simulating with dispersion rate = %f is not valid, please pick "
              "a positive number",
              cli_options.dispersion_rate.value());
    ok = false;
  }
  if (!cli_options.extinction_rate.has_value()) {
    MESSAGE_ERROR("Extinction rate was not set. Please provide a value for "
                  "the extinction rate");
  } else if (cli_options.extinction_rate.value() < 0) {
    LOG_ERROR("Simulating with dispersion rate = %f is not valid, please pick "
              "a positive number",
              cli_options.extinction_rate.value());
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
  } else if (!verify_path_is_readable(config_filename)) {
    LOG_ERROR("We don't have the permissions to read the config file %s",
              config_filename.c_str());
    ok = false;
  }

  return ok;
}

/**
 * Make the paths in a cli_options_t absolute, or at least simpler
 */
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
  } catch (const std::bad_optional_access &err) {
    return false;
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

/**
 * Check if results files exist already.
 *
 * Again, we use the philosophy that we should check as much as possible. So, we
 * check all the possible outputs, as long as they are set.
 */
[[nodiscard]] bool check_existing_results(const cli_options_t &cli_options) {
  bool ok = true;

  if (std::filesystem::exists(cli_options.phylip_filename())) {
    LOG_WARNING("Results file %s exists already",
                cli_options.phylip_filename().c_str());
    ok = false;
  }
  if (cli_options.yaml_file_set()
      && std::filesystem::exists(cli_options.yaml_filename())) {
    LOG_WARNING("Results file %s exists already",
                cli_options.yaml_filename().c_str());
    ok = false;
  }
  if (cli_options.json_file_set()
      && std::filesystem::exists(cli_options.json_filename())) {
    LOG_WARNING("Results file %s exists already",
                cli_options.json_filename().c_str());
    ok = false;
  }
  return ok;
}

/**
 * Convert a yaml file into cli_options_t
 */
cli_options_t parse_yaml_options(const std::filesystem::path &config_filename) {
  auto          yaml = YAML::LoadFile(config_filename);
  cli_options_t cli_options(yaml);

  return cli_options;
}

/**
 * Write the output as a YAML file.
 */
void write_yaml_file(std::ostream &os, const bigrig::tree_t &tree) {
  YAML::Emitter yaml;
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "tree" << YAML::Value << tree.to_newick();

  yaml << YAML::Key << "align";
  yaml << YAML::BeginMap;

  for (const auto &n : tree) {
    if (!n->is_leaf()) { continue; }
    yaml << YAML::Key << n->label() << YAML::Value << n->final_state().to_str();
  }
  yaml << YAML::EndMap;

  yaml << YAML::Key << "splits";
  yaml << YAML::BeginMap;
  for (const auto &n : tree) {
    if (n->is_leaf()) { continue; }
    yaml << YAML::Key << n->node_id();
    yaml << YAML::BeginMap;
    yaml << YAML::Key << "left" << YAML::Value << n->node_split().left.to_str();
    yaml << YAML::Key << "right" << YAML::Value
         << n->node_split().right.to_str();
    yaml << YAML::Key << "type" << YAML::Value
         << n->node_split().to_type_string();
    yaml << YAML::EndMap;
  }

  yaml << YAML::Key << "events";
  yaml << YAML::BeginMap;
  for (auto const &n : tree) {
    if (n->is_leaf()) { continue; }

    for (const auto &c : n->children()) {
      yaml << YAML::Key
           << std::format("{} -> {}", n->string_id(), c->string_id());
      yaml << YAML::BeginSeq;
      double total_time = 0;
      for (auto const &t : c->transitions()) {
        total_time += t.waiting_time;
        yaml << YAML::BeginMap;
        yaml << YAML::Key << "abs-time" << YAML::Value
             << n->abs_time() + total_time;
        yaml << YAML::Key << "waiting-time" << YAML::Value << t.waiting_time;
        yaml << YAML::Key << "initial-state" << YAML::Value
             << t.initial_state.to_str();
        yaml << YAML::Key << "final-state" << YAML::Value
             << t.final_state.to_str();
        yaml << YAML::EndMap;
      }
      yaml << YAML::EndSeq;
    }
  }
  yaml << YAML::EndMap;

  yaml << YAML::EndMap;
  os << yaml.c_str() << std::endl;
}

void write_json_file(std::ostream &os, const bigrig::tree_t &tree) {
  nlohmann::json j;

  j["tree"] = tree.to_newick();

  for (const auto &n : tree) {
    if (!n->is_leaf()) { continue; }

    j["align"][n->label()] = n->final_state().to_str();
  }

  for (const auto &n : tree) {
    if (n->is_leaf()) { continue; }
    auto split                                = n->node_split();
    j["splits"][std::to_string(n->node_id())] = {
        {"left", split.left.to_str()},
        {"right", split.right.to_str()},
        {"type", split.to_type_string()},
    };
  }

  for (const auto &n : tree) {
    if (n->is_leaf()) { continue; }
    for (const auto &c : n->children()) {
      auto   node_key = std::format("{} -> {}", n->string_id(), c->string_id());
      double total_time = 0;
      for (const auto &t : c->transitions()) {
        total_time += t.waiting_time;
        j["events"][node_key].push_back(
            {{"abs-time", n->abs_time() + total_time},
             {"waiting_time", t.waiting_time},
             {"initial-state", t.initial_state.to_str()},
             {"final-state", t.final_state.to_str()}});
      }
    }
  }

  os << j.dump() << std::endl;
}

/**
 * Write the output files given a sampled tree and model.
 *
 * Automatically selects which outputs need to be created based on
 * `cli_options_t`.
 */
void write_output_files(const cli_options_t          &cli_options,
                        const bigrig::tree_t         &tree,
                        const bigrig::biogeo_model_t &model) {
  auto phylip_filename  = cli_options.prefix.value();
  phylip_filename      += ".phy";
  std::ofstream phylip_file(phylip_filename);
  phylip_file << to_phylip(tree, model);

  auto phylip_all_filename  = cli_options.prefix.value();
  phylip_all_filename      += ".all.phy";
  std::ofstream phylip_all_file(phylip_all_filename);
  phylip_all_file << to_phylip_all_nodes(tree, model);

  auto cb = [](std::ostream &os, bigrig::node_t n) {
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
  annotated_tree_filename      += ".annotated";
  annotated_tree_filename      += bigrig::util::NEWICK_EXT;
  std::ofstream annotated_tree_file(annotated_tree_filename);
  annotated_tree_file << tree.to_newick(cb) << std::endl;

  if (cli_options.yaml_file_set()) {
    auto          output_yaml_filename = cli_options.yaml_filename();
    std::ofstream output_yaml_file(output_yaml_filename);
    write_yaml_file(output_yaml_file, tree);
  }
  if (cli_options.json_file_set()) {
    auto          output_json_filename = cli_options.json_filename();
    std::ofstream output_json_file(output_json_filename);
    write_json_file(output_json_file, tree);
  }
}

/**
 * Validate the CLI options, and merge the config file if passed
 *
 * This does much of the setup for the runtime of the program, including:
 * - Merging the config file with the current `cli_options_t`.
 * - Making directories
 * - Normalizing paths
 * - Checking existing results
 */
bool validate_options(cli_options_t &cli_options) {
  if (cli_options.config_filename.has_value()
      && config_compatible(cli_options)) {
    try {
      auto cli_options_tmp
          = parse_yaml_options(cli_options.config_filename.value());
      cli_options.merge(cli_options_tmp);
    } catch (const YAML::Exception &e) {
      LOG_ERROR("Failed to parse the config file: %s", e.what());
      return false;
    }
  }

  if (!cli_options.prefix.has_value() || cli_options.prefix.value().empty()) {
    cli_options.prefix = cli_options.tree_filename;
  }

  normalize_paths(cli_options);

  if (!validate_cli_options(cli_options)) {
    MESSAGE_ERROR(
        "We can't continue with the current options, exiting instead");
    return false;
  }

  write_header(cli_options);

  if (!check_existing_results(cli_options)) {
    if (!cli_options.redo.value_or(false)) {
      MESSAGE_ERROR("Refusing to run with existing results. Please specify the "
                    "--redo option if you want to overwrite existing results");
      return false;
    }
  }
  return true;
}

/**
 * Checks that the config file and the cli options are compatible
 */
[[nodiscard]] bool config_compatible(const cli_options_t &cli_options) {
  bool ok = true;

  if (cli_options.config_filename.has_value()
      && cli_options.cli_arg_specified()) {
    MESSAGE_ERROR("Both config and CLI arguments have been specified. Please "
                  "only specify one");
    ok = false;
  }
  if (cli_options.config_filename.has_value()) {
    auto &config_filename = cli_options.config_filename.value();
    if (!verify_config_file(config_filename)) {
      LOG_ERROR("There was an issue with the config file %s",
                config_filename.c_str());
      ok = false;
    }
  }

  return ok;
}
