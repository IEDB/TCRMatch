#include "substitution_matrix.hpp"
#include <iostream>
using namespace std;

// ---------------------------------------------------
array<array<float, 20>, 20> parse_matrix( const string& inputfilename )
{

	array<array<float, 20>, 20> submatrix;

	ifstream f( inputfilename );

	for (auto& row : submatrix) 
		for (auto& entry : row) 
			f >> entry;


	return submatrix;
}

// ---------------------------------------------------
array<array<float, 20>, 20> normalize( const array<array<float, 20>, 20>& matrix, float norm_factor )
{
	// Sum values of the row of aminoacids
	array<float,20> row_sum;
	for (int j = 0; j < 20; j++) {
		float sum = 0;
		for (int k = 0; k < 20; k++)
			sum += matrix[j][k];
		row_sum[j] = sum;
	}

	// Normalize matrix by marginal frequencies
	array<array<float, 20>, 20> norm_matrix; // Called k1 on the paper
	for (int j = 0; j < 20; j++) 
		for (int k = 0; k < 20; k++) 
			norm_matrix[j][k] = pow(matrix[j][k] / (row_sum[j] * row_sum[k]), norm_factor );

	return norm_matrix;
}

// ---------------------------------------------------
void print_matrix( const array<array<float, 20>, 20> &submatrix  )
{
	for (const auto& row : submatrix) 
	{
		for (const auto& entry : row) 
			cout << entry <<" ";
		cout << endl;
	}
}

// ---------------------------------------------------
// int main() 
// {
// 	array<array<float, 20>, 20> mat = parse_matrix( "data/blosum62.ssv" );
// 	print( mat );
// 	exit(0);
// 	cout <<mat[0].size();
// 	mat = normalize( mat, 0.11387);
// 	print( mat );
	
// }

