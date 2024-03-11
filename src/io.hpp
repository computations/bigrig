#pragma once

#include "clioptions.hpp"
#include "tree.hpp"

std::string to_phylip(const bigrig::tree_t &tree);

std::string to_phylip_all_nodes(const bigrig::tree_t &tree);

[[nodiscard]] bool config_compatible(const cli_options_t &cli_options);

bool validate_options(cli_options_t &cli_options);

void write_output_files(const cli_options_t  &cli_options,
                        const bigrig::tree_t &tree);

void write_header(const cli_options_t &cli_options);
