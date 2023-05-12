#include "valid_residues.hpp"

using namespace std;

string valid_residues = "ARNDCQEGHILKMFPSTWYV"; // Should be read from a file (submatrix?)

//  ------------------------------------------
char contain_invalid_residue( string peptide )
{
	for (const auto& residue : peptide) 
		if( valid_residues.find( residue ) == -1 )
			return residue;
	return 0;
}

//  ------------------------------------------

// int main()
// {
// 	// cout <<contain_invalid_residue("QWPKIMPLSARNDCQEGHILKMFPSTWYB") <<endl;
// 	// cout <<contain_invalid_residue("#") <<endl;
// 	char test = contain_invalid_residue("QWPKIMPLSARNDCQEGHILKMFPSTWY");
// 	if( test )
// 		cout <<"Contains invalid residue " <<test <<endl;
	
// 	test = contain_invalid_residue("QWP+P");
// 	if( test )
// 		cout <<"Contains invalid residue " <<test <<endl;
// }