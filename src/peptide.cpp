#include "peptide.hpp"
// #include <iostream>

using namespace std;

// -----------------------------------------
void trim_flanking_residues( peptide& pep )
{
	int last_idx = pep.seq.size() - 1;
  
	if ( pep.seq[0] == 'C' and ( pep.seq[last_idx] == 'F' or pep.seq[last_idx] == 'W' ) )
		pep.seq = pep.seq.substr( 1, last_idx-1 );
}

// -----------------------------------------




// int main()
// {
// 	peptide p;
// 	p.seq = "AQBLJADFPJ";
// 	cout <<p.seq <<endl;
// 	p.seq = "CAQBLJADFPJF";
// 	cout <<p.seq <<endl;
// 	trim_flanking_residues(p);
// 	cout <<p.seq <<endl;
// }