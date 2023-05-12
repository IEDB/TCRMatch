#ifndef __KMER_ALGORITHM_HPP
#define __KMER_ALGORITHM_HPP

#include <array>
#include <vector>
#include <string>
#include <map>
#include <array>

#include "peptide.hpp"
#include "database.hpp"
#include "valid_residues.hpp"

float k3_sum( const peptide& pep1, 
              const peptide& pep2, 
			  const std::array<std::array<float, 20>, 20>& k1, 
			  int p_kmin = 1, 
			  int p_kmax = 31 );


std::vector<std::string> multi_calc_k3( const std::vector<peptide>& peplist1, 
                                        const std::vector<peptide>& peplist2,
                                        float threshold,
							            std::map<std::string, std::vector<IEDB_data_row>>& iedb_map,
										std::array<std::array<float, 20>, 20>& k1 ); 

void calc_selfmatch_score( std::vector<peptide>& peplist, 
						   const std::array<std::array<float, 20>, 20>& smatrix );

#endif
