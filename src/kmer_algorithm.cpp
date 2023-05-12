#include <omp.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "substitution_matrix.hpp"
#include "kmer_algorithm.hpp"

using namespace std;


// --------------------------------------------------------
void calc_selfmatch_score( vector<peptide>& peplist, const array<array<float, 20>, 20>& smatrix )
{

	#pragma omp parallel for

	for( auto& p : peplist )
	{
		// this should be automatically done when initializing the peptide
		for( int x = 0; x < p.seq.size(); x++ ) 
			p.i.push_back( valid_residues.find( p.seq[x]));

		p.aff = k3_sum( p, p, smatrix );

	}

}


// --------------------------------------------------------
float k3_sum( const peptide& pep1, const peptide& pep2, const array<array<float, 20>, 20>& k1, int p_kmin, int p_kmax ) 
{
	// Calculate Kernel 3 using Kernel 1 lookups
	float k2, term, k3 = 0.0;
	int start1, start2;
	int k, j1, j2;

	float k2_prod_save[50][50][50]; // 31 is a problem for ATILYEILLGKATLYAVLVSALVLMATNEKLF, which is a 32-mer. 
	// This should be allocated dinamically, based on the dataset and query
	
	p_kmax = min( pep1.seq.size(), pep2.seq.size() );
	
	for( k = p_kmin; k <= p_kmax; k++ ) 
		for( start1 = 0; start1 <= pep1.seq.size() - k; start1++ ) 	
			for( start2 = 0; start2 <= pep2.seq.size() - k; start2++ ) 
			{
				j1 = pep1.i[start1 + k - 1];
				j2 = pep2.i[start2 + k - 1];
				term = k1[j1][j2];

				if (k == 1) 
					k2 = term;
				else 
					k2 = k2_prod_save[k - 1][start1][start2] * term;

				k2_prod_save[k][start1][start2] = k2;
				k3 += k2;
			}

	return( k3 );

}

// --------------------------------------------------------
vector<string> multi_calc_k3( const vector<peptide>& peplist1, 
                              const vector<peptide>& peplist2,
                              float threshold,
                              map<string, vector<IEDB_data_row>>& iedb_map, 
							  std::array<std::array<float, 20>, 20>& k1) 
{

	// cout <<peplist2[14612].seq <<" " <<peplist2[14612].aff <<std::endl;
	// cout <<peplist2[14613].seq <<" " <<peplist2[14613].aff <<endl;
	// for( auto& p : peplist2)
	// 	std::cout <<p.seq <<" " <<p.aff <<std::endl;
	// exit(0);
	// Simple method to calculate pairwise TCRMatch scores using two peptide
	// vectors
	vector<tuple<string, string, float, int>>::iterator it2[omp_get_max_threads()];
	vector<tuple<string, string, float, int>>
	results[omp_get_max_threads()];

	#pragma omp parallel for
	for (int i = 0; i < peplist1.size(); i++) 
    	for (int j = 0; j < peplist2.size(); j++) 
		{
      		peptide pep1 = peplist1[i];
      		peptide pep2 = peplist2[j];
      		float score = 0.0;
		
      		score = k3_sum(pep1, pep2, k1) / sqrt(pep1.aff * pep2.aff);
      		if (score > threshold) 
			{
        		int tid = omp_get_thread_num();
				// If input seq-match seq is unique, add tuple to results
				// Repeat IEDB matches are skipped (prevents duplicate rows in results) but duplicate input
				// sequences are permitted (e.g. for repertoires with identical TCRs)
				it2[tid] = find(results[tid].begin(), results[tid].end(), make_tuple(pep1.seq, pep2.seq, score, i));
				if (it2[tid] == results[tid].end()) 
					results[tid].push_back(make_tuple(pep1.seq, pep2.seq, score, i));
					
			}
	    }
  
//   cout << "input_sequence\tmatch_"
//                "sequence\tscore\treceptor_"
//                "group\tepitope\tantigen\torganism\t"
//             << endl;

	vector<string> tcrmatch_output;

	for (int k = 0; k < omp_get_max_threads(); k++) 
		for (auto &tuple : results[k]) 
		{
			string match_peptide = get<1>(tuple);
			map<string, vector<IEDB_data_row>>::iterator it = iedb_map.find(match_peptide);
			vector<IEDB_data_row> match_row_vec = it->second;
  
			for (int l = 0; l < match_row_vec.size(); l++) 
			{
				IEDB_data_row match_row = match_row_vec[l];
				//   cout << fixed << setprecision(4) << get<0>(tuple)
				//             << "\t" << get<1>(tuple) << "\t" << get<2>(tuple)
				//             << "\t" << match_row.receptor_group << "\t"
				//             << match_row.epitope << "\t" << match_row.antigen << "\t"
				//             << match_row.organism << endl;
				string score = to_string(get<2>(tuple));
				string output_string = get<0>(tuple) + "\t" + get<1>(tuple) + "\t" + score
							+ "\t" + match_row.receptor_group + "\t"
							+ match_row.epitope + "\t" + match_row.antigen + "\t"
							+ match_row.organism + "\n";

				tcrmatch_output.push_back( output_string );

      		}
    	}

  return tcrmatch_output;
}
