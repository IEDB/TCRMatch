#ifndef __DATABASE_HPP
#define __DATABASE_HPP


#include <string>
#include <vector>
#include <map>

// This should become DEPRECATED in the near future and replaced by a map<cdr_db, rest_of_row>
struct IEDB_data_row { 
  std::string amino_acid_sequence;
  std::string original_sequence;
  std::string receptor_group;
  std::string epitope;
  std::string organism;
  std::string antigen;
};


// TODO: check if all residues of the database are valid when parsing
std::vector<std::string> parse_database( std::string database_file ); // this is redundant; info also in iedb_map
std::map<std::string, std::vector<IEDB_data_row>> map_database( std::string IEDB_data_file );


#endif
