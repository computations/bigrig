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
  app.add_flag(
      "--json",
      [&cli_options](std::int64_t count) {
        cli_options.output_format_type = output_format_type_e::JSON;
      },
      "Output results in JSON, where possible");
  app.add_flag(
      "--yaml",
      [&cli_options](std::int64_t count) {
        cli_options.output_format_type = output_format_type_e::YAML;
      },
      "Output results in YAML, where possible");

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

  auto tree = biogeosim::tree_t(cli_options.tree_filename.value());

  pcg_extras::seed_seq_from<std::random_device> seed_source;
  pcg64_fast                                    gen(seed_source);

  biogeosim::substitution_model_t model(
      cli_options.dispersion_rate,
      cli_options.extinction_rate,
      cli_options.root_distribution.regions());

  tree.sample(cli_options.root_distribution, model, gen);

  write_output_files(cli_options, tree, model);

  return 0;
}
