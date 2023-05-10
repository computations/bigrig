#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>
#include <corax/corax.hpp>
#include <filesystem>
#include <logger.hpp>

struct cli_options_t {
  std::filesystem::path tree_filename;
};

int main(int argc, char **argv) {
  logger::get_log_states().add_stream(
      stdout,
      logger::log_level::info | logger::log_level::warning |
          logger::log_level::important | logger::log_level::error |
          logger::log_level::progress);

  CLI::App app{"Biogeosim"};

  cli_options_t cli_options;

  app.add_option("--tree",
                 cli_options.tree_filename,
                 "A file containing a newick encoded tree which will be used "
                 "to perform the simulation");

  CLI11_PARSE(app);

  auto tree =
      corax_utree_parse_newick_rooted(cli_options.tree_filename.c_str());

  std::vector<double> rates;

  return 0;
}
