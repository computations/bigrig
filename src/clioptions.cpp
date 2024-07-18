#include "clioptions.hpp"

#include "dist.hpp"
#include "logger.hpp"
#include "period.hpp"
#include "rng.hpp"
#include "util.hpp"

#include <algorithm>

pcg64_fast &cli_options_t::get_rng() { return bigrig::rng_wrapper_t::rng(); }
bigrig::rng_wrapper_t &cli_options_t::get_rng_wrapper() {
  return bigrig::rng_wrapper_t::get_instance();
}

std::filesystem::path cli_options_t::phylip_filename() const {
  auto tmp  = prefix.value();
  tmp      += bigrig::util::PHYILP_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::yaml_filename() const {
  auto tmp  = prefix.value();
  tmp      += bigrig::util::YAML_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::json_filename() const {
  auto tmp  = prefix.value();
  tmp      += bigrig::util::JSON_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::csv_splits_filename() const {
  constexpr auto splits_subprefix  = ".splits";
  auto           tmp               = prefix.value();
  tmp                             += splits_subprefix;
  tmp                             += bigrig::util::CSV_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::csv_events_filename() const {
  constexpr auto events_subprefix  = ".events";
  auto           tmp               = prefix.value();
  tmp                             += events_subprefix;
  tmp                             += bigrig::util::CSV_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::csv_periods_filename() const {
  constexpr auto state_subprefix  = ".periods";
  auto           tmp              = prefix.value();
  tmp                            += state_subprefix;
  tmp                            += bigrig::util::CSV_EXT;
  return tmp;
}

std::filesystem::path cli_options_t::csv_program_stats_filename() const {
  constexpr auto state_subprefix  = ".program-stats";
  auto           tmp              = prefix.value();
  tmp                            += state_subprefix;
  tmp                            += bigrig::util::CSV_EXT;
  return tmp;
}

/**
 * Checks if all the required args for the CLI have been specified.
 */
bool cli_options_t::cli_arg_specified() const {
  return tree_filename.has_value() || prefix.has_value()
      || debug_log.has_value() || tree_filename.has_value();
}

/**
 * Checks if the output format is a YAML file.
 */
bool cli_options_t::yaml_file_set() const {
  return output_format_type.has_value()
      && output_format_type.value() == output_format_type_e::YAML;
}

/**
 * Checks if the output format is a JSON file.
 */
bool cli_options_t::json_file_set() const {
  return output_format_type.has_value()
      && output_format_type.value() == output_format_type_e::JSON;
}

bool cli_options_t::csv_file_set() const {
  return output_format_type.has_value()
      && output_format_type.value() == output_format_type_e::CSV;
}

/**
 * Merges a `cli_options_t` with the current value. Specifically, it overwrites
 * the current values with the values from the passed `cli_options_t`. Values
 * affected are:
 *  - `config_filename`
 *  - `tree_filename`
 *  - `debug_log`
 *  - `output_format_type`
 *  - `root_distribution`
 *  - `dispersion_rate`
 *  - `extinction_rate`
 *  - `redo`
 *  - `two_region_duplicity`
 */

void print_config_cli_warning(const char *option_name) {
  LOG_WARNING("The '%s' option is specified in both the config file and "
              "the command line. Using the value from the command line",
              option_name);
}

template <typename T>
void merge_variable(std::optional<T>       &ours,
                    const std::optional<T> &theirs,
                    const char             *name) {
  if (theirs.has_value()) {
    if (ours.has_value()) {
      print_config_cli_warning(name);
    } else {
      ours = theirs;
    }
  }
}

void cli_options_t::merge(const cli_options_t &other) {
  merge_variable(tree_filename, other.tree_filename, "tree");
  merge_variable(prefix, other.prefix, "prefix");
  merge_variable(debug_log, other.debug_log, "debug-log");
  merge_variable(output_format_type, other.output_format_type, "output-format");
  merge_variable(root_range, other.root_range, "root-range");
  merge_variable(region_count, other.region_count, "region-count");

  /* periods are a bit different, so we handle them differently */
  if (!other.periods.empty()) {
    if (!periods.empty()) {
      print_config_cli_warning("periods");
    } else {
      periods = other.periods;
    }
  }

  merge_variable(redo, other.redo, "redo");
  merge_variable(
      two_region_duplicity, other.two_region_duplicity, "two-region-duplicity");
  merge_variable(mode, other.mode, "mode");
  merge_variable(rng_seed, other.rng_seed, "seed");
}

std::filesystem::path cli_options_t::get_tree_filename(const YAML::Node &yaml) {
  constexpr auto TREE_KEY = "tree";
  return yaml[TREE_KEY].as<std::string>();
}

std::optional<std::filesystem::path>
cli_options_t::get_prefix(const YAML::Node &yaml) {
  constexpr auto PREFIX_KEY = "prefix";
  if (yaml[PREFIX_KEY]) { return yaml[PREFIX_KEY].as<std::string>(); }
  return {};
}

std::optional<bool> cli_options_t::get_debug_log(const YAML::Node &yaml) {
  constexpr auto DEBUG_LOG_KEY = "debug-log";
  if (yaml[DEBUG_LOG_KEY]) { return yaml[DEBUG_LOG_KEY].as<bool>(); }
  return {};
}

std::optional<output_format_type_e>
cli_options_t::get_output_format(const YAML::Node &yaml) {
  constexpr auto OUPUT_FORMAT_KEY = "output-format";
  if (yaml[OUPUT_FORMAT_KEY]) {
    auto value = yaml[OUPUT_FORMAT_KEY].as<std::string>();
    std::transform(value.begin(),
                   value.end(),
                   value.begin(),
                   [](char c) -> char { return std::tolower(c); });

    if (value == "json") { return output_format_type_e::JSON; }
    if (value == "yaml") { return output_format_type_e::YAML; }
    if (value == "csv") { return output_format_type_e::CSV; }
  }
  return {};
}

std::optional<bool>
cli_options_t::get_two_region_duplicity(const YAML::Node &yaml) {
  constexpr auto DUPLICITY_KEY = "two-region-duplicity";
  if (yaml[DUPLICITY_KEY]) { return yaml[DUPLICITY_KEY].as<bool>(); }
  return {};
}

std::optional<bigrig::dist_t>
cli_options_t::get_root_range(const YAML::Node &yaml) {
  constexpr auto ROOT_DIST_KEY = "root-range";
  if (yaml[ROOT_DIST_KEY]) { return yaml[ROOT_DIST_KEY].as<std::string>(); }
  return {};
}

std::optional<size_t> cli_options_t::get_region_count(const YAML::Node &yaml) {
  constexpr auto REGION_COUNT_KEY = "region-count";
  if (yaml[REGION_COUNT_KEY]) { return yaml[REGION_COUNT_KEY].as<size_t>(); }
  return {};
}

std::optional<bigrig::rate_params_t>
cli_options_t::get_rates(const YAML::Node &yaml) {
  double         dis = 0.0, ext = 0.0;
  constexpr auto RATE_DICT_KEY = "rates";
  if (yaml[RATE_DICT_KEY]) {
    auto rates = yaml[RATE_DICT_KEY];

    constexpr auto DISPERSION_KEY = "dispersion";
    if (rates[DISPERSION_KEY]) {
      dis = rates[DISPERSION_KEY].as<double>();
    } else {
      MESSAGE_ERROR("No dispersion parameter provided for a period");
    }

    constexpr auto EXTINCTION_KEY = "extinction";
    if (rates[EXTINCTION_KEY]) {
      ext = rates[EXTINCTION_KEY].as<double>();
    } else {
      MESSAGE_ERROR("No extinction parameter provided for a period");
    }
    return {{dis, ext}};
  }
  return {};
}

std::optional<bigrig::cladogenesis_params_t>
cli_options_t::get_cladogenesis(const YAML::Node &yaml) {
  double allopatry_rate = 0.0, sympatry_rate = 0.0, copy_rate = 0.0,
         jump_rate              = 0.0;
  bool           ok             = true;
  constexpr auto CLADO_DICT_KEY = "cladogenesis";
  if (yaml[CLADO_DICT_KEY]) {
    auto clado = yaml[CLADO_DICT_KEY];

    constexpr auto ALLOPATRY_KEY = "allopatry";
    if (clado[ALLOPATRY_KEY]) {
      allopatry_rate = clado[ALLOPATRY_KEY].as<double>();
    } else {
      MESSAGE_ERROR("No allopatry parameter provided for a period");
      ok = false;
    }

    constexpr auto SYMPATRY_KEY = "sympatry";
    if (clado[SYMPATRY_KEY]) {
      sympatry_rate = clado[SYMPATRY_KEY].as<double>();
    } else {
      MESSAGE_ERROR("No sympatry parameter provided for a period");
      ok = false;
    }

    constexpr auto COPY_KEY = "copy";
    if (clado[COPY_KEY]) {
      copy_rate = clado[COPY_KEY].as<double>();
    } else {
      MESSAGE_ERROR("No copy parameter provided for a period");
      ok = false;
    }

    constexpr auto JUMP_KEY = "jump";
    if (clado[JUMP_KEY]) {
      jump_rate = clado[JUMP_KEY].as<double>();
    } else {
      MESSAGE_ERROR("No jump parameter provided for a period");
      ok = false;
    }

    if (ok) { return {{allopatry_rate, sympatry_rate, copy_rate, jump_rate}}; }
  }
  return {};
}

std::optional<period_params_t>
cli_options_t::get_period(const YAML::Node &yaml) {
  period_params_t period_params;
  bool            ok        = true;
  constexpr auto  START_KEY = "start";
  if (yaml[START_KEY]) {
    period_params.start = yaml[START_KEY].as<double>();
  } else {
    MESSAGE_ERROR("No start time provided for a period");
    ok = false;
  }

  auto rates = get_rates(yaml);
  if (rates.has_value()) {
    period_params.rates = rates.value();
  } else {
    MESSAGE_ERROR("Rates for a period are malformed");
    ok = false;
  }

  auto clado_params = get_cladogenesis(yaml);
  if (clado_params.has_value()) {
    period_params.clado = clado_params.value();
  } else {
    MESSAGE_ERROR("Cladogenesis parameters for a period are malformed");
    ok = false;
  }

  if (ok) { return period_params; }
  return {};
}

std::vector<period_params_t>
cli_options_t::get_periods(const YAML::Node &yaml) {
  std::vector<period_params_t> ret;
  constexpr auto               PERIOD_KEY = "periods";
  if (yaml[PERIOD_KEY]) {
    auto   list  = yaml[PERIOD_KEY];
    size_t index = 0;
    for (auto y : list) {
      auto period = get_period(y);
      if (!period.has_value()) {
        LOG_ERROR("Period %lu is malformed", index);
        return {};
      }
      period.value().index = index++;
      ret.push_back(period.value());
    }
  } else {
    auto rates        = get_rates(yaml);
    auto clado_params = get_cladogenesis(yaml);
    bool ok           = true;

    if (!rates.has_value()) {
      ok = false;
      MESSAGE_ERROR("Failed to find rates in the config file");
    }
    if (!clado_params.has_value()) {
      ok = false;
      MESSAGE_ERROR("Failed to find cladogensis parameters in the config file");
    }

    if (ok) {
      return {period_params_t{.rates = rates.value(),
                              .clado = clado_params.value(),
                              .start = 0.0,
                              .index = 0}};
    }
  }
  return ret;
}

std::optional<bool> cli_options_t::get_redo(const YAML::Node &yaml) {
  constexpr auto REDO_KEY = "redo";
  if (yaml[REDO_KEY]) { return yaml[REDO_KEY].as<bool>(); }
  return {};
}

std::optional<bigrig::operation_mode_e>
cli_options_t::get_mode(const YAML::Node &yaml) {
  constexpr auto MODE_KEY = "mode";
  if (yaml[MODE_KEY]) {
    std::string value = yaml[MODE_KEY].as<std::string>();
    std::transform(value.begin(),
                   value.end(),
                   value.begin(),
                   [](char c) -> char { return std::tolower(c); });
    if (value == "fast") {
      return bigrig::operation_mode_e::FAST;
    } else if (value == "sim") {
      return bigrig::operation_mode_e::SIM;
    } else {
      throw cli_option_invalid_parameter{
          "Failed to recognize the run mode in the config file"};
    }
  }
  return {};
}

std::optional<uint64_t> cli_options_t::get_seed(const YAML::Node &yaml) {
  constexpr auto SEED_KEY = "seed";
  if (yaml[SEED_KEY]) { return yaml[SEED_KEY].as<uint64_t>(); }
  return {};
}

template <typename T>
[[nodiscard]] bool check_passed_cli_parameter(const std::optional<T> &o,
                                              const char             *name) {
  if (!o.has_value()) {
    LOG_ERROR("No value provided for parameter '%s' when other values were "
              "provided on the command line",
              name);
    return false;
  }
  return true;
}

[[nodiscard]] bool
cli_options_t::convert_cli_parameters(std::optional<double> dis,
                                      std::optional<double> ext,
                                      std::optional<double> allo,
                                      std::optional<double> symp,
                                      std::optional<double> copy,
                                      std::optional<double> jump) {
  if (!(dis.has_value() || ext.has_value() || allo.has_value()
        || symp.has_value() || copy.has_value() || jump.has_value())) {
    return true;
  }
  bool ok  = true;
  ok      &= check_passed_cli_parameter(dis, "dispersion");
  ok      &= check_passed_cli_parameter(ext, "extinction");
  ok      &= check_passed_cli_parameter(allo, "allopatry");
  ok      &= check_passed_cli_parameter(symp, "sympatry");
  ok      &= check_passed_cli_parameter(copy, "copy");
  ok      &= check_passed_cli_parameter(jump, "jump");

  if (ok) {
    period_params_t period;
    period.start = 0.0;
    period.rates = {.dis = dis.value(), .ext = ext.value()};
    period.clado = {.allopatry = allo.value(),
                    .sympatry  = ext.value(),
                    .copy      = copy.value(),
                    .jump      = jump.value()};
    periods.push_back(period);
    return true;
  }
  return false;
}

std::vector<bigrig::period_t> cli_options_t::make_periods() const {
  std::vector<bigrig::period_t> ret;
  for (const auto &cli_period : periods) {
    bigrig::period_t period{
        cli_period.start,
        0.0,
        {.dis = cli_period.rates.dis, .ext = cli_period.rates.ext},
        {.allopatry = cli_period.clado.allopatry,
         .sympatry  = cli_period.clado.sympatry,
         .copy      = cli_period.clado.copy,
         .jump      = cli_period.clado.jump},
        two_region_duplicity.value_or(true),
        cli_period.index};
    if (!period.model().check_ok(root_range.value().regions())) {
      LOG_ERROR("There is an issue with the model for period '%lu', we can't "
                "continue",
                cli_period.index);
      return {};
    }
    ret.push_back(period);
  }

  for (size_t i = 0; i < ret.size() - 1; ++i) {
    auto &p1 = ret[i];
    auto &p2 = ret[i + 1];
    p1.set_length(p2.start() - p1.start());
  }
  ret.back().set_length(std::numeric_limits<double>::infinity());
  return ret;
}
