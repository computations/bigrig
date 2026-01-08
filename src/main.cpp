#include "clioptions.hpp"
#include "dist.hpp"
#include "io.hpp"
#include "model.hpp"
#include "pcg_random.hpp"

#include <chrono>
#include <corax/corax.hpp>
#include <logger.hpp>

inline bigrig::tree_t get_tree(const cli_options_t &cli_options) {
  if (!cli_options.simulate_tree.value_or(false)) {
    return bigrig::tree_t(cli_options.tree_filename.value());
  }
  return bigrig::tree_t{};
}

int main(int argc, char **argv) {
  const auto start_time{std::chrono::high_resolution_clock::now()};

  logger::get_log_states().add_stream(stdout, logger::info | logger::defaults);

  CLI::App app{"A tool to simulate (ancestal) range distributions under the "
               "DEC[+J] model."};

  cli_options_t         cli_options;
  std::optional<double> dispersion, extinction, allopatry, sympatry, copy, jump;

  app.add_option("--config",
                 cli_options.config_filename,
                 "YAML file containing the program configuration. Information "
                 "about the config file can be found in the readme.");
  app.add_option("--tree",
                 cli_options.tree_filename,
                 "A file containing a newick encoded tree which "
                 "will be used to perform the simulation.");
  app.add_option(
      "--prefix", cli_options.prefix, "Prefix for the output files.");
  app.add_option(
      "--root-range",
      cli_options.root_range,
      "Range for the species at the root for the start of the "
      "simulation. Should be given as a binary string (e.g. 01010). Required "
      "if region-count is not specified.");
  app.add_option("--region-count",
                 cli_options.region_count,
                 "Specify the number of regions to simulate. If "
                 "given, and root-dist is not given, then a random root "
                 "distribution is generated and used");
  app.add_option(
      "-d,--dispersion", dispersion, "The dispersion rate for the simulation.");
  app.add_option(
      "-e,--extinction", extinction, "The extinction rate for the simulation.");

  app.add_option("-v,--allopatry",
                 allopatry,
                 "The allopatry/vicariance rate for cladogenesis "
                 "for the simulation.");
  app.add_option("-s,--sympatry",
                 sympatry,
                 "The sympatry rate for cladogenesis for the simulation.");
  app.add_option(
      "-y,--copy", copy, "The copy rate for cladogenesis for the simulation.");
  app.add_option(
      "-j,--jump", jump, "The jump rate for cladogenesis for the simulation.");
  app.add_option("--seed", cli_options.rng_seed, "Seed for the RNG");

  app.add_flag("--redo", cli_options.redo, "Ignore existing result files");
  app.add_flag("--debug-log",
               cli_options.debug_log,
               "Create a file in the prefix that contains the debug "
               "log. Don't enable this without a good reason.");
  app.add_flag(
      "--json",
      [&cli_options](std::int64_t count) {
        (void)(count); // Silence a warning
        cli_options.output_format_type = output_format_type_e::JSON;
      },
      "Output results in a JSON file.");
  app.add_flag(
      "--yaml",
      [&cli_options](std::int64_t count) {
        (void)(count); // Silence a warning
        cli_options.output_format_type = output_format_type_e::YAML;
      },
      "Output results in a YAML file.");
  app.add_flag(
      "--csv",
      [&cli_options](std::int64_t count) {
        (void)(count); // Silence a warning
        cli_options.output_format_type = output_format_type_e::CSV;
      },
      "Output results in a YAML file.");
  app.add_flag("--two-region-duplicity",
               cli_options.two_region_duplicity,
               "[Optional] Allow for outcome duplicity in the case of 2 region "
               "splits. See the README.md for more information.")
      ->group("");
  app.add_flag(
      "--sim",
      [&cli_options](std::int64_t count) {
        (void)(count);
        cli_options.mode = bigrig::operation_mode_e::SIM;
      },
      "Run in simulation mode (warning: slow).");
  app.add_flag(
      "--fast",
      [&cli_options](std::int64_t count) {
        (void)(count);
        cli_options.mode = bigrig::operation_mode_e::FAST;
      },
      "Run in fast mode (default on).");

  CLI11_PARSE(app, argc, argv);

  if (!cli_options.convert_cli_parameters(
          dispersion, extinction, allopatry, sympatry, copy, jump)) {
    LOG_ERROR("If model parameters are passed on the command line, all of "
                  "the parameters must be provided");
    return 1;
  }

  if (!validate_and_finalize_options(cli_options)) {
    LOG_ERROR("Use --help to get a list of all options");
    return 1;
  }

  if (cli_options.debug_log) {
    std::filesystem::path debug_filename  = cli_options.prefix.value();
    debug_filename                       += ".debug.log";
    LOG_INFO("Logging debug information to {}", debug_filename.c_str());
    logger::get_log_states().add_file_stream(
        debug_filename.c_str(),
        logger::info | logger::debug | logger::defaults);
  }

  auto gen = cli_options.get_rng();

  if (cli_options.simulate_tree.value_or(false)) {
    LOG_INFO("Parsing tree");
  }

  auto tree = get_tree(cli_options);
  tree.set_mode(cli_options.mode.value_or(bigrig::operation_mode_e::FAST));

  bool ok = true;

  auto periods = cli_options.make_periods(gen);
  if (!periods.validate(cli_options.compute_region_count())) {
    LOG_ERROR("There was an issue with the specified periods");
    ok = false;
  }

  if (periods.empty()) { return 1; }
  tree.set_periods(periods);

  if (!tree.is_ready(cli_options.simulate_tree.value_or(false))) {
    LOG_ERROR("Could not use the tree provided");
    ok = false;
  }

  if (!ok) { return 1; }

  const auto config_time{std::chrono::high_resolution_clock::now()};

  if (!cli_options.simulate_tree.value_or(false)) {
    LOG_INFO("Tree has {} taxa", tree.leaf_count());
    LOG_INFO("Simulating ranges on the tree");
    tree.simulate(cli_options.root_range.value(), gen);
  } else {
    periods.set_extinction(true);
    LOG_INFO("Simulating ranges and tree");
    tree.simulate_tree(cli_options.root_range.value(),
                       periods,
                       cli_options.tree_height.value_or(1.0),
                       true /*prune tree after simulation*/,
                       gen);
    LOG_INFO("Simulated tree with {} taxa", tree.leaf_count());
  }
  const auto      end_time{std::chrono::high_resolution_clock::now()};
  program_stats_t program_stats{.start_time  = start_time,
                                .config_time = config_time,
                                .end_time    = end_time};

  LOG_INFO("Writing results to files");
  write_output_files(cli_options, tree, periods, program_stats);

  LOG_INFO("Done!");
  return 0;
}
