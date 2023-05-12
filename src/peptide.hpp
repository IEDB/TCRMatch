#ifndef __PEPTIDE_HPP
#define __PEPTIDE_HPP

#include <string>
#include <vector>

struct peptide 
{
  std::string seq;
  int len;
  float aff;
  std::vector<int> i; // same as string seq, but with integers indicating the position in the matrix
};

void trim_flanking_residues( peptide& pep );

#endif