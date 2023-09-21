#include "node.hpp"
#include "tree.hpp"
#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <corax/corax.hpp>
#include <filesystem>
#include <logger.hpp>
#include <optional>

struct cli_options_t {
  std::filesystem::path tree_filename;
  std::optional<std::filesystem::path> debug_filename;
};

int main(int argc, char **argv) {
  logger::get_log_states().add_stream(
      stdout, logger::log_level::info | logger::log_level::warning |
                  logger::log_level::important | logger::log_level::error |
                  logger::log_level::progress);

  CLI::App app{"Biogeosim"};

  cli_options_t cli_options;

  app.add_option("--tree", cli_options.tree_filename,
                 "A file containing a newick encoded tree which will be used "
                 "to perform the simulation")
      ->required();
  app.add_option("--debug-file", cli_options.debug_filename,
                 "A file to write the full debug log to");

  CLI11_PARSE(app);

  if (cli_options.debug_filename.has_value()) {
    logger::get_log_states().add_file_stream(
        cli_options.debug_filename.value().c_str(),
        logger::log_level::info | logger::log_level::warning |
            logger::log_level::important | logger::log_level::error |
            logger::log_level::debug);
  }

  MESSAGE_DEBUG("this is the debug message");

  auto tree = biogeosim::tree_t(cli_options.tree_filename.c_str());

  std::random_device rd;
  std::minstd_rand gen(rd());
  biogeosim::substitution_model_t model(1.0, 1.0, 6);

  tree.sample(0110011, model, gen);
  std::cout << tree.to_newick() << std::endl;

  return 0;
}
