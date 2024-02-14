#pragma once

#include "clioptions.hpp"
#include "model.hpp"
#include "tree.hpp"

std::string to_phylip(const bigrig::tree_t        &tree,
                      const bigrig::biogeo_model_t model);

std::string to_phylip_all_nodes(const bigrig::tree_t         &tree,
                                const bigrig::biogeo_model_t &model);

[[nodiscard]] bool config_compatible(const cli_options_t &cli_options);

bool validate_options(cli_options_t &cli_options);

void write_output_files(const cli_options_t          &cli_options,
                        const bigrig::tree_t         &tree,
                        const bigrig::biogeo_model_t &model);

void write_header(const cli_options_t &cli_options);
