#include "io.hpp"

#include "clioptions.hpp"
#include "logger.hpp"
#include "model.hpp"

#include <cmath>
#include <csv.hpp>
#include <expected>
#include <filesystem>
#include <optional>
#include <string_view>

using namespace std::string_view_literals; // for the 'sv' suffix

constexpr size_t MAX_REGIONS = 64;

void print_periods(const std::vector<bigrig::period_params_t> &periods) {
  LOG_INFO("   Running with %lu periods:", periods.size());
  for (const auto &p : periods) {
    LOG_INFO("      - Start time: %.2f", p.start);
    MESSAGE_INFO("        Rate parameters:");
    LOG_INFO("          Dispersion(d): %.2f, Extinction(e): %.2f",
             p.rates.dis,
             p.rates.ext);
    MESSAGE_INFO("        Cladogenesis parameters:");
    LOG_INFO("          Allopatry(v): %.2f, Sympatry(s): %.2f, Copy(y): %.2f, "
             "Jump(j): %.2f",
             p.clado.allopatry,
             p.clado.sympatry,
             p.clado.copy,
             p.clado.jump);
  }
}

void print_model_parameters(const bigrig::period_params_t &period) {
  MESSAGE_INFO("   Model Parameters:");
  MESSAGE_INFO("     Rate parameters:");
  LOG_INFO("       Dispersion(d): %.2f, Extinction(e): %.2f",
           period.rates.dis,
           period.rates.ext);
  MESSAGE_INFO("     Cladogenesis parameters:");
  LOG_INFO("       Allopatry(v): %.2f, Sympatry(s): %.2f, Copy(y): %.2f, "
           "Jump(j): %.2f",
           period.clado.allopatry,
           period.clado.sympatry,
           period.clado.copy,
           period.clado.jump);
}

/**
 * Print the message at the start of the run.
 */
void write_header(const cli_options_t &cli_options) {
  MESSAGE_INFO("Running simulation with the following options:");
  if (cli_options.tree_filename.has_value()) {
    LOG_INFO("   Tree file: %s", cli_options.tree_filename.value().c_str());
  } else {
    MESSAGE_INFO("   Tree: Simulate");
  }
  LOG_INFO("   Prefix: %s", cli_options.prefix.value().c_str());
  LOG_INFO("   Root range: %s",
           cli_options.root_range.value().to_str().c_str());
  LOG_INFO("   Region count: %u", cli_options.root_range->regions());
  if (cli_options.periods.size() == 1) {
    print_model_parameters(cli_options.periods.front());
  } else {
    print_periods(cli_options.periods);
  }
  if (cli_options.rng_seed.has_value()) {
    LOG_INFO("   Seed: %lu", cli_options.rng_seed.value());
  }
  if (cli_options.mode.has_value()
      && cli_options.mode.value() == bigrig::operation_mode_e::SIM) {
    MESSAGE_WARNING(
        "Setting the operation mode to simulation, results will be slow");
  }
}

/**
 * Produce a phylip file as a string.
 */
std::string to_phylip(const bigrig::tree_t &tree) {
  std::ostringstream oss;
  oss << std::to_string(tree.leaf_count()) << " " << tree.region_count()
      << "\n";

  tree.to_phylip_body(oss);

  return oss.str();
}

/**
 * Produce a phylip file as a string, including inner nodes.
 */
std::string to_phylip_all_nodes(const bigrig::tree_t &tree) {
  std::ostringstream oss;
  oss << std::to_string(tree.node_count()) << " " << tree.region_count()
      << "\n";

  tree.to_phylip_body(oss, true);

  return oss.str();
}

/**
 * Check that the tree file path is ok to use.
 */
[[nodiscard]] bool validate_tree_filename(
    const std::optional<std::filesystem::path> &tree_filename_option) {
  if (!tree_filename_option.has_value()) {
    MESSAGE_WARNING("No tree file was provided");
    return true;
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

/**
 * Checks that the prefix is alright, and then make the directories.
 *
 * When the user passes a prefix, we need to check that its all kosher, and
 * then _actually_ make the directories. We make them here.
 */
[[nodiscard]] bool validate_and_make_prefix(
    const std::optional<std::filesystem::path> &prefix_option) {
  if (!prefix_option.has_value()) {
    MESSAGE_ERROR("No prefix was provided");
    return false;
  }
  auto prefix = prefix_option.value();
  bool ok     = true;

  if (prefix.has_parent_path()
      && !std::filesystem::exists(prefix.parent_path())) {
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

[[nodiscard]] bool validate_model_parameter(const std::optional<double> &param,
                                            const char                  *name) {
  bool ok = true;
  if (!param.has_value()) {
    LOG_ERROR("The model parameter '%s' was not set. Please provide a "
              "value for this parameter",
              name);
    ok = false;
  } else if (param.value() < 0) {
    LOG_ERROR("Simulating with '%s' = %f is not valid, please pick "
              "a positive number",
              name,
              param.value());
    ok = false;
  }
  return ok;
}

[[nodiscard]] bool
validate_root_region(const std::optional<bigrig::dist_t> &root_range,
                     const std::optional<size_t>         &region_count) {
  bool ok = true;
  if (!root_range.has_value() && !region_count.has_value()) {
    MESSAGE_ERROR("The root range was not provided. Please provide a value for "
                  "the root range");
    return false;
  }
  if (root_range.has_value() && region_count.has_value()
      && root_range.value().regions() != region_count.value()) {
    MESSAGE_ERROR("Both a root range and region count was provided, but they "
                  "differ in size");
    ok = false;
  }
  if (root_range.has_value()) {
    if (root_range.value().regions() >= MAX_REGIONS) {
      LOG_ERROR("Simulating with %u regions is unsupported. Please choose a "
                "number less than %lu regions",
                root_range.value().regions(),
                MAX_REGIONS);
      ok = false;
    }
    if (root_range.value().empty()) {
      MESSAGE_ERROR("Cannot simulate with an empty root range. Please provide "
                    "a range with at least one region set");
      ok = false;
    }
  }
  if (region_count.has_value() && region_count.value() >= MAX_REGIONS) {
    LOG_ERROR("region-count is set to %lu, but that number of regions is "
              "unsupported. Please a choose a number less than %lu "
              "regions",
              region_count.value(),
              MAX_REGIONS);
    ok = false;
  }
  return ok;
}

[[nodiscard]] bigrig::adjustment_matrix_symmetry
determine_matrix_symmetry(const std::vector<bigrig::adjacency_arc_t> &matrix,
                          size_t region_count) {
  size_t symmetric_size = (region_count * (region_count + 1)) / 2;

  if (matrix.size() == symmetric_size) {
    return bigrig::adjustment_matrix_symmetry::symmetric;
  }

  size_t nonsymmetric_size = (region_count - 1) * region_count;

  if (matrix.size() == nonsymmetric_size) {
    return bigrig::adjustment_matrix_symmetry::nonsymmetric;
  }
  return bigrig::adjustment_matrix_symmetry::unknown;
}

[[nodiscard]] bool
validate_matrix_symmetry(const std::vector<bigrig::adjacency_arc_t> &matrix) {
  for (auto a : matrix) {
    bool row_symmetry = false;
    for (auto b : matrix) { row_symmetry |= a.reverse(b); }
    if (!row_symmetry) { return false; }
  }
  return true;
}

[[nodiscard]] bool
validiate_adjustment_matrix(const std::vector<bigrig::adjacency_arc_t> &matrix,
                            size_t region_count) {
  /*
   * there are two cases:
   *  - Symmetric adjustment matrix
   *  - Nonsymmetric adjustment matrix
   * In both of these cases, the diagonal does nothing, so we should not include
   * it. Therfore, the two formulas are
   *  - n(n-1)
   *  - n(n+1)/2
   * So, we need to determine which case it is, and then check the matrix size.
   *
   * This data structure at this moment is a list of tuples:
   *  (from, to, value)
   * So, to check if the matrix is intended to be symmetric, we need to check if
   * (a, b, _) and (b, a, _) are in the list.
   */

  auto symmetry = determine_matrix_symmetry(matrix, region_count);

  switch (symmetry) {
  case bigrig::adjustment_matrix_symmetry::symmetric:
    if (!validate_matrix_symmetry(matrix)) {
      MESSAGE_ERROR(
          "A matrix is not fully symmetric, despite being the correct size");
      return false;
    }
    return true;
  case bigrig::adjustment_matrix_symmetry::nonsymmetric:
    return true;

  default:
  case bigrig::adjustment_matrix_symmetry::unknown:
    return false;
  }
}

[[nodiscard]] bool validate_adjustment_matrix_params(
    const std::optional<bigrig::adjustment_matrix_params_t>
        &adjustment_params) {
  if (!adjustment_params.has_value()) { return true; }
  bool ok = true;

  auto params = adjustment_params.value();

  if (params.simulate.has_value() && params.simulate.value()
      && params.adjustments.has_value()) {
    MESSAGE_ERROR("Both an adjustment matrix and the simulate option were set. "
                  "These are incompatible");
    ok &= false;
  }

  if (params.exponent.has_value()) {
    double &exponent = params.exponent.value();
    if (!std::isfinite(exponent)) {
      MESSAGE_ERROR("There is an issue with the adjustment matrix exponent");
      ok &= false;
    }
  }

  return ok;
}

/**
 * Check that the program options are valid
 *
 * Check that the program options work. Note, it tries to be thorough when
 * checking, as the loop of "change a thing, find something else wrong" is
 * annoying. So, this function tries to check as much as it can, and not to
 * bail out at the first error.
 */
[[nodiscard]] bool validate_cli_options(const cli_options_t &cli_options) {
  bool ok = true;

  ok &= validate_tree_filename(cli_options.tree_filename);
  ok &= validate_and_make_prefix(cli_options.prefix);
  ok &= validate_root_region(cli_options.root_range, cli_options.region_count);

  for (const auto &p : cli_options.periods) {
    ok &= validate_model_parameter(p.rates.dis, "dispersion");
    ok &= validate_model_parameter(p.rates.ext, "extinction");

    ok &= validate_model_parameter(p.clado.allopatry, "allopatry");
    ok &= validate_model_parameter(p.clado.sympatry, "sympatry");
    ok &= validate_model_parameter(p.clado.copy, "copy");
    ok &= validate_model_parameter(p.clado.jump, "jump");

    ok &= validate_adjustment_matrix_params(p.adjustment_matrix);
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
  } catch (const std::bad_optional_access &err) { return false; }
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

std::expected<std::vector<bigrig::adjacency_arc_t>, bigrig::io_err>
read_adjustment_matrix(const bigrig::adjustment_matrix_params_t &params) {
  auto &filename = *params.matrix_filename;

  std::vector<bigrig::adjacency_arc_t> rows;

  if (!std::filesystem::exists(filename)) {
    LOG_ERROR("The matrix file '%s' does not exist", filename.c_str());
    return std::unexpected{bigrig::io_err::ReadError};
  }
  if (!verify_path_is_readable(filename)) {
    return std::unexpected{bigrig::io_err::ReadError};
  }

  csv::CSVReader reader{filename.string()};

  for (auto &cur_row : reader) {
    rows.push_back(bigrig::adjacency_arc_t{
        .from  = cur_row[0].get(),
        .to    = cur_row[1].get(),
        .value = cur_row[2].get<double>(),
    });
  }

  std::sort(rows.begin(), rows.end(), [](const auto &a, const auto &b) {
    return a.from < b.from && a.to < b.to;
  });

  return rows;
}

/**
 * Check if results files exist already.
 *
 * Again, we use the philosophy that we should check as much as possible. So,
 * we check all the possible outputs, as long as they are set.
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
  if (cli_options.csv_file_set()) {
    if (std::filesystem::exists(cli_options.csv_splits_filename())) {
      LOG_WARNING("Results file %s exists already",
                  cli_options.csv_splits_filename().c_str());
      ok = false;
    }
    if (std::filesystem::exists(cli_options.csv_events_filename())) {
      LOG_WARNING("Results file %s exists already",
                  cli_options.csv_events_filename().c_str());
      ok = false;
    }
    if (std::filesystem::exists(cli_options.csv_periods_filename())) {
      LOG_WARNING("Results file %s exists already",
                  cli_options.csv_periods_filename().c_str());
      ok = false;
    }
    if (std::filesystem::exists(cli_options.csv_program_stats_filename())) {
      LOG_WARNING("Results file %s exists already",
                  cli_options.csv_program_stats_filename().c_str());
      ok = false;
    }
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

template <typename T, typename U>
void write_yaml_value(YAML::Emitter &yaml, const T &label, const U &value) {
  yaml << YAML::Key << label << YAML::Value << value;
}

template <typename T>
void write_yaml_value(YAML::Emitter        &yaml,
                      const T              &label,
                      const bigrig::dist_t &value) {
  yaml << YAML::Key << label << YAML::Value << value.to_str();
}

void write_yaml_tree(YAML::Emitter &yaml, const bigrig::tree_t &tree) {
  write_yaml_value(yaml, "tree", tree.to_newick());
}

void write_yaml_taxa(YAML::Emitter &yaml, const bigrig::tree_t &tree) {
  write_yaml_value(yaml, "taxa", tree.leaf_count());
}

void write_yaml_regions(YAML::Emitter &yaml, const size_t &regions) {
  write_yaml_value(yaml, "region-count", regions);
}

void write_yaml_root_range(YAML::Emitter        &yaml,
                           const bigrig::dist_t &root_dist) {
  write_yaml_value(yaml, "root-range", root_dist);
}

void write_yaml_alignment(YAML::Emitter &yaml, const bigrig::tree_t &tree) {
  yaml << YAML::Key << "align";
  yaml << YAML::BeginMap;

  for (const auto &n : tree) {
    yaml << YAML::Key << n->string_id() << YAML::Value
         << n->final_state().to_str();
  }
  yaml << YAML::EndMap;
}

void write_yaml_splits(YAML::Emitter &yaml, const bigrig::tree_t &tree) {
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
    yaml << YAML::Key << "period" << YAML::Value
         << n->node_split().period_index;
    yaml << YAML::EndMap;
  }
  yaml << YAML::EndMap;
}

void write_yaml_events(YAML::Emitter &yaml, const bigrig::tree_t &tree) {
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
        yaml << YAML::Key << "period" << YAML::Value << t.period_index;
        yaml << YAML::EndMap;
      }
      yaml << YAML::EndSeq;
    }
  }
  yaml << YAML::EndMap;
}

void write_yaml_period(YAML::Emitter &yaml, const bigrig::period_t &period) {
  auto &model = period.model();
  yaml << YAML::BeginMap;

  yaml << YAML::Key << "start";
  yaml << YAML::Value << period.start();

  auto [dis, ext] = model.rates();

  yaml << YAML::Key << "rates";
  yaml << YAML::Value;

  yaml << YAML::BeginMap;

  yaml << YAML::Key << "dispersion";
  yaml << YAML::Value << dis;

  yaml << YAML::Key << "extinction";
  yaml << YAML::Value << ext;

  yaml << YAML::EndMap;

  auto [allo, symp, copy, jump] = model.cladogenesis_params();

  yaml << YAML::Key << "cladogenesis";
  yaml << YAML::Value;

  yaml << YAML::BeginMap;
  yaml << YAML::Key << "allopatry";
  yaml << YAML::Value << allo;

  yaml << YAML::Key << "sympatry";
  yaml << YAML::Value << symp;

  yaml << YAML::Key << "copy";
  yaml << YAML::Value << copy;

  yaml << YAML::Key << "jump";
  yaml << YAML::Value << jump;
  yaml << YAML::EndMap;

  yaml << YAML::EndSeq;
}

void write_yaml_period_list(YAML::Emitter               &yaml,
                            const bigrig::period_list_t &periods) {
  yaml << YAML::Key << "periods";
  yaml << YAML::Value;
  yaml << YAML::BeginSeq;

  for (const auto &p : periods) { write_yaml_period(yaml, p); }

  yaml << YAML::EndMap;
}

void write_yaml_program_stats(YAML::Emitter         &yaml,
                              const program_stats_t &program_stats) {
  yaml << YAML::Key << "stats";
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "time";
  yaml << YAML::Value;
  yaml << program_stats.execution_time_in_seconds();
  yaml << YAML::EndMap;
}

/**
 * Write the output as a YAML file.
 */
void write_yaml_file(std::ostream                &os,
                     const bigrig::tree_t        &tree,
                     const bigrig::period_list_t &periods,
                     const program_stats_t       &program_stats) {
  YAML::Emitter yaml;
  yaml << YAML::BeginMap;

  write_yaml_tree(yaml, tree);
  write_yaml_regions(yaml, tree.region_count());
  write_yaml_root_range(yaml, tree.get_root_range());
  write_yaml_alignment(yaml, tree);
  write_yaml_splits(yaml, tree);
  write_yaml_events(yaml, tree);
  write_yaml_period_list(yaml, periods);
  write_yaml_program_stats(yaml, program_stats);

  yaml << YAML::EndMap;
  os << yaml.c_str() << std::endl;
}

void write_json_file(std::ostream                &os,
                     const bigrig::tree_t        &tree,
                     const bigrig::period_list_t &periods,
                     const program_stats_t       &program_stats) {
  nlohmann::json j;

  j["tree"]          = tree.to_newick();
  j["taxa"]          = tree.leaf_count();
  j["regions"]       = tree.region_count();
  j["root-range"]    = tree.get_root_range().to_str();
  j["stats"]["time"] = program_stats.execution_time_in_seconds();

  for (const auto &n : tree) {
    j["align"][n->string_id()] = n->final_state().to_str();
  }

  for (const auto &n : tree) {
    if (n->is_leaf()) { continue; }
    auto split                  = n->node_split();
    j["splits"][n->string_id()] = {
        {"left", split.left.to_str()},
        {"right", split.right.to_str()},
        {"type", split.to_type_string()},
        {"period", split.period_index},
    };
  }

  for (const auto &n : tree) {
    if (n->is_leaf()) { continue; }
    for (const auto &c : n->children()) {
      auto   node_key = std::format("{} -> {}", n->string_id(), c->string_id());
      double total_time = 0;
      for (const auto &t : c->transitions()) {
        total_time += t.waiting_time;
        j["events"][node_key].push_back({
            {"abs-time", n->abs_time() + total_time},
            {"waiting_time", t.waiting_time},
            {"initial-state", t.initial_state.to_str()},
            {"final-state", t.final_state.to_str()},
            {"period", t.period_index},
        });
      }
    }
  }

  for (const auto &p : periods) {
    auto &model                   = p.model();
    auto [dis, ext]               = model.rates();
    auto [allo, symp, copy, jump] = model.cladogenesis_params();

    j["periods"].push_back({{"start", p.start()},
                            {"rates",
                             {
                                 {"dispersion", dis},
                                 {"extinction", ext},
                             }},
                            {"cladogenesis",
                             {
                                 {"allopatry", allo},
                                 {"sympatry", symp},
                                 {"copy", copy},
                                 {"jump", jump},
                             }}});
  }

  os << j.dump() << std::endl;
}

template <typename T, size_t N>
inline std::string make_csv_row(const std::array<T, N> &entries) {
  return std::accumulate(std::next(entries.begin()),
                         entries.end(),
                         std::string(*entries.begin()),
                         [](const T &acc, const T &entry) -> std::string {
                           return acc + ", " + entry;
                         })
       + "\n";
}

template <size_t N>
inline std::string make_csv_row(const std::array<std::string_view, N> &fields) {
  return std::accumulate(
             std::next(fields.begin()),
             fields.end(),
             std::string(*fields.begin()),
             [](std::string acc, const std::string_view &entry) -> std::string {
               acc += ", ";
               acc += entry;
               return acc;
             })
       + "\n";
}

template <size_t N>
inline std::string make_csv_row(const std::array<double, N> &fields) {
  return std::accumulate(
             std::next(fields.begin()),
             fields.end(),
             std::to_string(*fields.begin()),
             [](std::string acc, const double &entry) -> std::string {
               acc += ", ";
               acc += std::to_string(entry);
               return acc;
             })
       + "\n";
}

template <size_t N>
inline std::ofstream init_csv(const std::filesystem::path           &filename,
                              const std::array<std::string_view, N> &fields) {
  std::ofstream csv_file(filename);
  csv_file << make_csv_row(fields);
  return csv_file;
}

void write_split_csv_file(const cli_options_t  &cli_options,
                          const bigrig::tree_t &tree) {
  auto                 output_filename = cli_options.csv_splits_filename();
  constexpr std::array fields{
      "node"sv, "left"sv, "right"sv, "type"sv, "period"sv};
  auto output_file = init_csv(output_filename, fields);

  for (const auto &n : tree) {
    if (n->is_leaf()) { continue; }
    auto split = n->node_split();

    output_file << make_csv_row(std::array{n->string_id(),
                                           split.left.to_str(),
                                           split.right.to_str(),
                                           split.to_type_string(),
                                           std::to_string(split.period_index)});
  };
}

void write_events_csv_file(const cli_options_t  &cli_options,
                           const bigrig::tree_t &tree) {
  auto                 output_filename = cli_options.csv_events_filename();
  constexpr std::array fields{"node"sv,
                              "waiting-time"sv,
                              "initial-state"sv,
                              "final-state"sv,
                              "period"sv};
  auto                 output_file = init_csv(output_filename, fields);

  for (const auto &n : tree) {
    for (const auto &t : n->transitions()) {
      output_file << make_csv_row(std::array{n->string_id(),
                                             std::to_string(t.waiting_time),
                                             t.initial_state.to_str(),
                                             t.final_state.to_str(),
                                             std::to_string(t.period_index)});
    }
  }
}

void write_periods_csv_file(const cli_options_t         &cli_options,
                            const bigrig::period_list_t &periods) {
  auto                 output_filename = cli_options.csv_periods_filename();
  constexpr std::array fields{"index"sv,
                              "start"sv,
                              "dispersion"sv,
                              "extinction"sv,
                              "allopatry"sv,
                              "sympatry"sv,
                              "copy"sv,
                              "jump"sv};
  auto                 output_file = init_csv(output_filename, fields);

  for (const auto &p : periods) {
    auto &model                   = p.model();
    auto [dis, ext]               = model.rates();
    auto [allo, symp, copy, jump] = model.cladogenesis_params();

    output_file << make_csv_row(std::array{static_cast<double>(p.index()),
                                           p.start(),
                                           dis,
                                           ext,
                                           allo,
                                           symp,
                                           copy,
                                           jump});
  }
}

void write_program_stats_csv_file(const cli_options_t   &cli_options,
                                  const program_stats_t &program_stats) {
  auto output_filename = cli_options.csv_program_stats_filename();
  constexpr std::array fields{"stat"sv, "value"sv};

  auto output_file = init_csv(output_filename, fields);
  output_file << make_csv_row(std::array<std::string, 2>{
      "execution-time",
      std::to_string(program_stats.execution_time_in_seconds())});
}

void write_csv_files(const cli_options_t         &cli_options,
                     const bigrig::tree_t        &tree,
                     const bigrig::period_list_t &periods,
                     const program_stats_t       &program_stats) {
  write_split_csv_file(cli_options, tree);
  write_events_csv_file(cli_options, tree);
  write_periods_csv_file(cli_options, periods);
  write_program_stats_csv_file(cli_options, program_stats);
}

void write_clean_tree_file(const cli_options_t  &cli_options,
                           const bigrig::tree_t &tree) {
  auto clean_cb = [](std::ostream &os, bigrig::node_t n) {
    os << n.string_id();
    os << ":" << n.brlen();
  };

  auto clean_tree_filename  = cli_options.prefix.value();
  clean_tree_filename      += ".clean";
  clean_tree_filename      += bigrig::util::NEWICK_EXT;
  std::ofstream clean_tree_file(clean_tree_filename);
  clean_tree_file << tree.to_newick(clean_cb) << std::endl;
}

void write_annotated_tree_file(const cli_options_t  &cli_options,
                               const bigrig::tree_t &tree) {
  auto annotated_cb = [](std::ostream &os, bigrig::node_t n) {
    os << n.string_id();
    os << ":" << n.brlen();
    os << "[&&NHX:";
    if (n.is_leaf()) {
      os << "dist=" << n.final_state().to_str();
    } else {
      os << n.node_split().to_nhx_string();
    }
    os << "]";
  };

  auto annotated_tree_filename  = cli_options.prefix.value();
  annotated_tree_filename      += ".annotated";
  annotated_tree_filename      += bigrig::util::NEWICK_EXT;
  std::ofstream annotated_tree_file(annotated_tree_filename);
  annotated_tree_file << tree.to_newick(annotated_cb) << std::endl;
}

/**
 * Write the output files given a sampled tree and model.
 *
 * Automatically selects which outputs need to be created based on
 * `cli_options_t`.
 */
void write_output_files(const cli_options_t         &cli_options,
                        const bigrig::tree_t        &tree,
                        const bigrig::period_list_t &periods,
                        const program_stats_t       &program_stats) {
  auto          phylip_filename = cli_options.phylip_filename();
  std::ofstream phylip_file(phylip_filename);
  phylip_file << to_phylip(tree);

  auto phylip_all_filename  = cli_options.prefix.value();
  phylip_all_filename      += ".all.phy";
  std::ofstream phylip_all_file(phylip_all_filename);
  phylip_all_file << to_phylip_all_nodes(tree);

  write_annotated_tree_file(cli_options, tree);

  if (cli_options.simulate_tree.value_or(false)) {
    write_clean_tree_file(cli_options, tree);
  }

  if (cli_options.yaml_file_set()) {
    auto          output_yaml_filename = cli_options.yaml_filename();
    std::ofstream output_yaml_file(output_yaml_filename);
    write_yaml_file(output_yaml_file, tree, periods, program_stats);
  }
  if (cli_options.json_file_set()) {
    auto          output_json_filename = cli_options.json_filename();
    std::ofstream output_json_file(output_json_filename);
    write_json_file(output_json_file, tree, periods, program_stats);
  }
  if (cli_options.csv_file_set()) {
    write_csv_files(cli_options, tree, periods, program_stats);
  }
}

[[nodiscard]] bool finalize_options(cli_options_t &cli_options) {
  bool ok = true;

  if (cli_options.rng_seed.has_value()) {
    cli_options.get_rng_wrapper().seed(cli_options.rng_seed.value());
  } else {
    cli_options.get_rng_wrapper().seed();
  }
  if (!cli_options.root_range.has_value()) {
    cli_options.root_range = bigrig::make_random_dist(
        cli_options.region_count.value(), cli_options.get_rng());
  }

  for (auto [index, p] : std::views::enumerate(cli_options.periods)) {
    if (p.adjustment_matrix.has_value()) {
      auto res = read_adjustment_matrix(p.adjustment_matrix.value());
      if (!res) {
        LOG_ERROR("Could not read the matrix file for period %lu", index);
        ok &= false;
      } else {
        auto &matrix = *res;
        if (!validiate_adjustment_matrix(matrix,
                                         cli_options.root_range->regions())) {
          LOG_ERROR("The matrix was malformed for period %lu", index);
          ok &= false;
          continue;
        }
        p.adjustment_matrix->adjustments = *res;
      }
    }
  }

  return ok;
}

/**
 * Validate the CLI options, and merge the config file if passed
 *
 * This does much of the setup for the runtime of the program, including:
 * - Merging the config file with the current `cli_options_t`.
 * - Making directories
 * - Normalizing paths
 * - Checking existing results
 *   Loading the adjustment matrices
 */
bool validate_and_finalize_options(cli_options_t &cli_options) {
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

  if (!finalize_options(cli_options)) {
    MESSAGE_ERROR("Failed to finalize the setup exiting");
    return false;
  };

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
