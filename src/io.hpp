#pragma once

#include "clioptions.hpp"
#include "period.hpp"
#include "tree.hpp"

std::string to_phylip(const bigrig::tree_t &tree);

std::string to_phylip_all_nodes(const bigrig::tree_t &tree);

[[nodiscard]] bool config_compatible(const cli_options_t &cli_options);

bool validate_and_finalize_options(cli_options_t &cli_options);

void write_output_files(const cli_options_t         &cli_options,
                        const bigrig::tree_t        &tree,
                        const bigrig::period_list_t &period,
                        const program_stats_t       &program_stats);

void write_header(const cli_options_t &cli_options);
