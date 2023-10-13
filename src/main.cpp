#include "model.hpp"
#include "node.hpp"
#include "pcg_random.hpp"
#include "tree.hpp"

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <corax/corax.hpp>
#include <filesystem>
#include <logger.hpp>
#include <optional>
#include <sstream>
#include <string>

struct cli_options_t {
  std::filesystem::path                tree_filename;
  std::filesystem::path                prefix;
  std::optional<std::filesystem::path> debug_filename;
};

std::string to_phylip(const biogeosim::tree_t              &tree,
                      const biogeosim::substitution_model_t model) {
  std::ostringstream oss;
  oss << std::to_string(tree.leaf_count()) << " " << model.region_count()
      << "\n";

  tree.to_phylip_body(oss, model.region_count());

  return oss.str();
}

int main(int argc, char **argv) {
  logger::get_log_states().add_stream(
      stdout,
      logger::log_level::info | logger::log_level::warning
          | logger::log_level::important | logger::log_level::error
          | logger::log_level::progress);

  CLI::App app{"Biogeosim"};

  cli_options_t cli_options;

  app.add_option("--tree",
                 cli_options.tree_filename,
                 "A file containing a newick encoded tree which will be used "
                 "to perform the simulation")
      ->required();
  app.add_option("--prefix", cli_options.prefix, "prefix for the output files");
  app.add_option("--debug-file",
                 cli_options.debug_filename,
                 "A file to write the full debug log to");

  CLI11_PARSE(app);

  if (cli_options.debug_filename.has_value()) {
    logger::get_log_states().add_file_stream(
        cli_options.debug_filename.value().c_str(),
        logger::log_level::info | logger::log_level::warning
            | logger::log_level::important | logger::log_level::error
            | logger::log_level::debug);
  }

  MESSAGE_DEBUG("this is the debug message");

  auto tree = biogeosim::tree_t(cli_options.tree_filename.c_str());

  pcg_extras::seed_seq_from<std::random_device> seed_source;
  pcg64_fast                                    gen(seed_source);

  biogeosim::substitution_model_t model(1.0, 1.0, 6);

  tree.sample({0b11'0011, 6}, model, gen);
  std::cout << tree.to_newick() << std::endl;
  std::cout << to_phylip(tree, model) << std::endl;

  return 0;
}
