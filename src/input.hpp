#ifndef __INPUT_HPP
#define __INPUT_HPP

#include <string>
#include <fstream>
#include <vector>

#include "peptide.hpp"

void parse_raw_input( std::vector<peptide>& peplist, std::string const& in_file, bool trimming = true );
void parse_airr_input( std::vector<peptide>& peplist, std::string const& in_file, bool trimming = true );
void parse_trust4_input( std::vector<peptide>& peplist, std::string const& in_file, std::vector<std::string>& inputlines, bool trimming = true );
#endif