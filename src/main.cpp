#include "clioptions.hpp"
#include "dist.hpp"
#include "io.hpp"
#include "model.hpp"
#include "pcg_extras.hpp"
#include "pcg_random.hpp"

#include <corax/corax.hpp>
#include <logger.hpp>

int main() {
  logger::get_log_states().add_stream(
      stdout,
      logger::log_level::info | logger::log_level::warning
          | logger::log_level::important | logger::log_level::error
          | logger::log_level::progress);

  CLI::App app{"A tool to simulate (ancestal) range distributions under the "
               "DEC[+J] model."};

  cli_options_t         cli_options;
  std::optional<double> dispersion, extinction, allopatry, sympatry, copy, jump;

  app.add_option("--config",
                 cli_options.config_filename,
                 "[Optional] YAML file containing the program configuration.");
  app.add_option("--tree",
                 cli_options.tree_filename,
                 "[Required] A file containing a newick encoded tree which "
                 "will be used to perform the simulation.");
  app.add_option("--prefix",
                 cli_options.prefix,
                 "[Optional] Prefix for the output files.");
  app.add_option(
      "--root-dist",
      cli_options.root_distribution,
      "[Required] Range for the species at the root for the start of the "
      "simulation. Should be given as a binary string (e.g. 01010).");
  app.add_option("-d,--dispersion",
                 dispersion,
                 "[Required] The dispersion rate for the simulation.");
  app.add_option("-e,--extinction",
                 extinction,
                 "[Required] The extinction rate for the simulation.");

  app.add_option("-v,--allopatry",
                 allopatry,
                 "[Required] The allopatry/vicariance rate for cladogenesis "
                 "for the simulation.");
  app.add_option(
      "-s,--sympatry",
      sympatry,
      "[Required] The sympatry rate for cladogenesis for the simulation.");
  app.add_option(
      "-y,--copy",
      copy,
      "[Required] The copy rate for cladogenesis for the simulation.");
  app.add_option(
      "-j,--jump",
      jump,
      "[Required] The jump rate for cladogenesis for the simulation.");
  app.add_option("--seed", cli_options.rng_seed, "[Optional] Seed for the RNG");

  app.add_flag(
      "--redo", cli_options.redo, "[Optional] Ignore existing result files");
  app.add_flag("--debug-log",
               cli_options.debug_log,
               "[Optional] Create a file in the prefix that contains the debug "
               "log. Don't enable this without a good reason.");
  app.add_flag(
      "--json",
      [&cli_options](std::int64_t count) {
        (void)(count); // Silence a warning
        cli_options.output_format_type = output_format_type_e::JSON;
      },
      "[Optional] Output results in a JSON file.");
  app.add_flag(
      "--yaml",
      [&cli_options](std::int64_t count) {
        (void)(count); // Silence a warning
        cli_options.output_format_type = output_format_type_e::YAML;
      },
      "[Optional] Output results in a YAML file.");
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
      "[Optional] Run in simulation mode (warning: slow).");
  app.add_flag(
      "--fast",
      [&cli_options](std::int64_t count) {
        (void)(count);
        cli_options.mode = bigrig::operation_mode_e::FAST;
      },
      "[Optional] Run in fast mode.");

  CLI11_PARSE(app);

  if (!cli_options.convert_cli_parameters(
          dispersion, extinction, allopatry, sympatry, copy, jump)) {
    MESSAGE_ERROR("If model parameters are passed on the command line, all of "
                  "the parameters must be provided");
    return 1;
  }

  if (!validate_options(cli_options)) { return 1; }

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

  MESSAGE_INFO("Parsing tree");
  auto tree = bigrig::tree_t(cli_options.tree_filename.value());
  tree.set_mode(cli_options.mode.value_or(bigrig::operation_mode_e::FAST));

  auto periods = cli_options.make_periods();
  if (periods.empty()) { return 1; }
  tree.set_periods(periods);

  if (!tree.is_ready()) {
    MESSAGE_ERROR("Could not use the tree provided, exiting");
    return 1;
  }

  LOG_INFO("Tree has %lu taxa", tree.leaf_count());

  pcg64_fast gen;
  if (cli_options.rng_seed.has_value()) {
    gen.seed(cli_options.rng_seed.value());
  } else {
    gen.seed(pcg_extras::seed_seq_from<std::random_device>{});
  }

  MESSAGE_INFO("Simulating ranges on the tree")
  tree.simulate(cli_options.root_distribution.value(), gen);

  MESSAGE_INFO("Writing results to files")
  write_output_files(cli_options, tree, periods);

  MESSAGE_INFO("Done!")
  return 0;
}
