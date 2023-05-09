#include "cli.hpp"
#include <corax/corax.hpp>
#include <logger.hpp>

int main(int argc, char **argv) {

  logger::get_log_states().add_stream(
      stdout,
      logger::log_level::info | logger::log_level::warning |
          logger::log_level::important | logger::log_level::error |
          logger::log_level::progress);

  auto cli_options = cli_options_t(argc, argv);
  MESSAGE_INFO("Hello World");
  LOG_INFO("You passed '%s' for --tree",
           cli_options["tree"].value<std::filesystem::path>().c_str());

  return 0;
}
