#define VERSION "1.0.3"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <omp.h>
#include <string>
#include <tuple>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <sstream>

#include "input.hpp"
#include "peptide.hpp"
#include "valid_residues.hpp"
#include "substitution_matrix.hpp" 
#include "kmer_algorithm.hpp"
#include "database.hpp"


void trim( std::string& seq )
{
	int last_idx = seq.size() - 1;
  
	if ( seq[0] == 'C' and ( seq[last_idx] == 'F' or seq[last_idx] == 'W' ) )
		seq = seq.substr( 1, last_idx-1 );
}

// -----------------------------------------------------------------
std::vector<peptide> iedb_data_to_peplist( std::vector<std::string>& iedb_data )
{
  // Comment 1. Change to IEDB data: Can't we get this from IEDB map????
  // Comment 2. Also, why do we even need 2 lists? Why not use the map? check if it impacts performance
  // Comment 3: Checking should happen as database is parsed
  
  std::vector<peptide> peplist;
  for (std::vector<std::string>::iterator it = iedb_data.begin(); it != iedb_data.end(); it++) 
  {
    std::vector<int> int_vec;
    for (int i = 0; i < (*it).length(); i++) {
      if (valid_residues.find((*it)[i]) == -1) {
        std::cerr << "Invalid amino acid found in " << *it << " at position " << i + 1 << std::endl;
        exit(1);
      }
    }
    peplist.push_back({*it, int((*it).length()), -99.9, int_vec});
  }
  return peplist;
}

// -----------------------------------------------------------------
std::string get_nth_field( std::string& str, int n ){
	std::istringstream iss(str);
	std::string field;
	for (int i = 0; i <= n; i++)
		iss >>field;
	return field;
}

// -----------------------------------------------------------------
// std::array<std::array<float, 20>, 20> blm_qij; 
// int p_kmin = 1;   // turn into parameters
// int p_kmax = 30;  // turn into parameters 


// std::array<std::array<float, 20>, 20> k1;// It's global because when it isn't, it gives a diff result
// Move this to outside -> import everything you need
int main(int argc, char *argv[]) {
  std::array<std::array<float, 20>, 20> k1; // Must be global otherwise the result changes (don't know why)
  float p_beta = 0.11387;
  int opt;
  int n_threads;
  float threshold;
  std::string in_file;
  std::string iedb_file = "data/IEDB_data.tsv";
  int i_flag = -1;
  int t_flag = -1;
  int thresh_flag = -1;
  int airr_flag = -1;
  bool trimming = true;
  bool trust4flag = false;
  // Command line argument parsing
  while ((opt = getopt(argc, argv, "rakvt:i:s:d:")) != -1) {
    switch (opt) {
      case 't':
        n_threads = std::stoi(optarg);
        t_flag = 1;
        break;
      case 'a':
        airr_flag = 1;
        break;
      case 'i':
        in_file = optarg;
        i_flag = 1;
        break;
      case 's':
        threshold = std::stof(optarg);
        thresh_flag = 1;
        break;
      case 'd':
        iedb_file = optarg;
        break;
      case 'k':
        trimming = false;
        break;
      case 'r':
        trust4flag = true;
        break;
      case 'v':
        std::cout << VERSION << std::endl;
        return 0;
      default:
        std::cerr << "Usage: ./tcrmatch -i infile_name.txt -a -t num_threads -s "
                    "score_threshold -d /path/to/database -k -r"
                  << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Check that required parameters are there + error correcting
  if (i_flag == -1 || t_flag == -1) {
    std::cerr << "Missing mandatory parameters" << std::endl
              << "Usage: ./tcrmatch -i infile_name.txt -t num_threads"
              << std::endl;
    return EXIT_FAILURE;
  }
  if (n_threads < 1) {
    std::cerr << "Number of threads cannot be less than one." << std::endl;
    return EXIT_FAILURE;
  }
  if (thresh_flag == -1) {
    threshold = .97;
  }
  if (t_flag == -1) {
    n_threads = 1;
  }
  if (threshold < 0 || threshold > 1) {
    std::cerr << "Threshold must be between 0 and 1" << std::endl;
    return EXIT_FAILURE;
  }

  std::map<std::string, std::vector<IEDB_data_row>> iedb_map = map_database(iedb_file);

  std::vector<std::string> iedb_data = parse_database(iedb_file);
  
  // Check if we didn't read any data into IEDB data // RT: remove this: check if file is valid inside
  if (iedb_data.size() == 0) {
    std::cout
        << "Failed to read database file - please specify the path to database "
           "file if not ./data/IEDB_data.tsv"
        << std::endl;

    return EXIT_FAILURE;
  }


  
  omp_set_num_threads(n_threads);

  // Parses and normalizes matrix
  k1 = parse_matrix( "data/blosum62.ssv" );
  k1 = normalize( k1, p_beta);
  // std::cout <<std::endl;
  // print_matrix(k1);
  
  std::vector<std::string> inputlines;  

  // Parse input file
  std::vector<peptide> peplist1;
  if (airr_flag == 1) 
    parse_airr_input( peplist1, in_file );  

  else if( trust4flag ) 
  {
    parse_trust4_input( peplist1, in_file, inputlines, trimming );
  }

  else 
  	parse_raw_input( peplist1, in_file, trimming );    
  // std::cout <<std::endl;
  // print_matrix(k1);
  // Parse database file
  std::vector<peptide> peplist2 = iedb_data_to_peplist( iedb_data );
  
  // std::cout <<peplist1[0].aff <<std::endl;  
  calc_selfmatch_score( peplist1, k1 );

  // for( auto& p : peplist1 )
  //   std::cout << p.aff <<std::endl;
  // exit(0);
  calc_selfmatch_score( peplist2, k1);
  // std::cout <<std::endl;
  // print_matrix(k1);
  // Calculate input data vs database with multi-threading
  std::vector<std::string> extended_output;

  // std::cout <<peplist2[14612].seq <<" " <<peplist2[14612].aff <<std::endl;
  // std::cout <<peplist2[14613].seq <<" " <<peplist2[14613].aff <<std::endl;
  // exit(0);
  
  std::vector<std::string> tcrmatch_output = multi_calc_k3( peplist1, peplist2, threshold, iedb_map, k1 );

  // Output: this is ugly as sin and must be refactored
  std::string tcrmatch_output_header;
  std::string was_trimmed;
  if( trimming )
    was_trimmed = "trimmed_";
  else 
    was_trimmed = "";

  tcrmatch_output_header = was_trimmed + "input_sequence\tmatch_sequence\tscore\treceptor_group\tepitope\tantigen\torganism" ;

  if( trust4flag == true )
  {
    std::ifstream file1(in_file);
    std::string trust4header;
    getline(file1, trust4header); // parses header
    file1.close();

    std::cout << trust4header <<"\t";
    std::cout << tcrmatch_output_header << std::endl; 

    // Appending TCRmatch to Trust4 output
    for( int i = 0; i < inputlines.size(); i++)
      for( int j = 0; j < tcrmatch_output.size(); j++) {
        std::string final_line = inputlines[i] + "\t" + tcrmatch_output[j];
        std::string AAseq_trust  = get_nth_field(inputlines[i], 3);
        std::string seq_tcrmatch = get_nth_field(tcrmatch_output[j], 0);
        // std::cout << (trim(AAseq_trust) == seq_tcrmatch) <<std::endl;
        if( trimming )
          trim(AAseq_trust);

        if( AAseq_trust == seq_tcrmatch )
	        extended_output.push_back( final_line );
    }

    for( int i = 0; i < extended_output.size(); i++)
      std::cout <<extended_output[i];
  }

  // #Change from 'input_sequence' to 'trimmed_input_sequence' if trimmed

  else {
    std::cout << tcrmatch_output_header << std::endl;
    
    for( int j = 0; j < tcrmatch_output.size(); j++) 
      std::cout<< tcrmatch_output[j];
  }
  return 0;
}
