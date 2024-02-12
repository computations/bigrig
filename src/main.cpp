#include "clioptions.hpp"
#include "io.hpp"
#include "model.hpp"
#include "node.hpp"
#include "pcg_random.hpp"
#include "tree.hpp"

#include <corax/corax.hpp>
#include <fstream>
#include <iostream>
#include <logger.hpp>
#include <sstream>
#include <string>

int main() {
  logger::get_log_states().add_stream(
      stdout,
      logger::log_level::info | logger::log_level::warning
          | logger::log_level::important | logger::log_level::error
          | logger::log_level::progress);

  CLI::App app{"A tool to simulate (ancestal) range distribution under the "
               "DEC[+J] model."};

  cli_options_t cli_options;

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
      "simulation. Should be given as a binary string (e.g. 1010).");
  app.add_option("--d",
                 cli_options.dispersion_rate,
                 "[Required] The dispersion rate for the simulation.");
  app.add_option("--e",
                 cli_options.extinction_rate,
                 "[Required] The extinction rate for the simulation.");
  app.add_flag(
         "--redo", cli_options.redo, "[Optional] Ignore existing result files")
      ->default_val("yes");
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
               "splits. See the README.md for more information.");

  CLI11_PARSE(app);

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

  if (!tree.is_valid()) {
    MESSAGE_ERROR("The tree provided is not valid, exiting");
    return 1;
  }

  if (!tree.is_binary()) {
    MESSAGE_ERROR("The tree provided is not a binary tree, refusing to run");
    return 1;
  }

  pcg_extras::seed_seq_from<std::random_device> seed_source;
  pcg64_fast                                    gen(seed_source);

  bigrig::biogeo_model_t model(cli_options.dispersion_rate.value(),
                               cli_options.extinction_rate.value(),
                               cli_options.root_distribution.value().regions(),
                               cli_options.two_region_duplicity.value_or(true));

  MESSAGE_INFO("Sampling from tree")
  tree.simulate(cli_options.root_distribution.value(), model, gen);

  MESSAGE_INFO("Writing results to files")
  write_output_files(cli_options, tree, model);

  MESSAGE_INFO("Done!")
  return 0;
}
