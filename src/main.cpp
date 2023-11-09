#define VERSION "1.1.2"

#include "tcrmatch.cpp"

int main(int argc, char *argv[]) {
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

  std::map<std::string, std::vector<IEDB_data_row>> iedb_map =
      create_IEDB_map(iedb_file);

  std::vector<std::string> iedb_data = read_IEDB_data(iedb_file);
  // Check if we didn't read any data into IEDB data
  if (iedb_data.size() == 0) {
    std::cout
        << "Failed to read database file - please specify the path to database "
           "file if not ./data/IEDB_data.tsv"
        << std::endl;

    return EXIT_FAILURE;
  }
  std::string line;
  std::string seq;
  std::string alphabet;
  std::vector<peptide> peplist1;
  std::vector<peptide> peplist2;
  std::vector<std::string> inputlines;  
  std::string trust4header;
  omp_set_num_threads(n_threads);

  alphabet = "ARNDCQEGHILKMFPSTWYV";
  k1 = fmatrix_k1();
  if (airr_flag == 1) {
    // AIRR input format
    std::vector<std::string> airr_data = read_AIRR_data(in_file);
    if (airr_data.empty())
      return EXIT_FAILURE;
    for (std::vector<std::string>::iterator it = airr_data.begin();
         it != airr_data.end(); it++) {
      std::vector<int> int_vec;
      bool invalid_residue = false;

      for (int i = 0; i < (*it).length(); i++) {
         if (alphabet.find((*it)[i]) == -1) {
          std::cerr << "Invalid amino acid found in " << *it << " at position "
                    << i + 1 <<". Skipping." << std::endl;
            invalid_residue = true;
            break;
        //   return EXIT_FAILURE;
        }
      }
      if(invalid_residue)
         continue;
      peplist1.push_back({*it, int((*it).length()), -99.9, int_vec});
    }
  } 
  
  else if( trust4flag ){
	
    // text file input
    std::ifstream file1(in_file);

	
    getline(file1, trust4header); // skips header 
    erase_newline( trust4header );

    std::string valid_header = "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\tcid_full_length";
    
    // Checking if TRUST4 file is correctly formatted
    if ( trust4header != valid_header )
    {
      std::cerr <<"ERROR: Invalid TRUST4 header. Please check input file format." <<std::endl ;
      std::cerr <<"Input file: " + in_file <<std::endl ;
      exit(1);
    }

    while (getline(file1, line)) {
  		
      erase_newline( line );
      
      if( is_TCR_gene(line) )
      {
      
        inputlines.push_back( line );

        std::istringstream is( line );
        std::string skip;
        is >>skip >>skip >>skip >>seq;  // 4th field contains peptide seq
        std::vector<int> int_vec;
        bool invalid_char_found = false;

        for (int i = 0; i < seq.length(); i++) {
          if (alphabet.find(seq[i]) == -1) {
            // return EXIT_FAILURE; // TRUST4 contains many illegal chars. Skipping them for now.
            invalid_char_found = true;
            break;
          }
        }

        if( invalid_char_found ) {
          std::cerr << "Invalid amino acid found in " << seq <<std::endl;
          std::cerr <<"Skipping" <<std::endl;
          continue;
        }

        if( trimming ) // removes flanking residues (C and F or W). Works for text input, not AIRR.
          seq = trim( seq );

        peplist1.push_back({seq, int(seq.length()), -99.9, int_vec});
      }
    }
    file1.close();
  }

  else {
    // text file input
    std::ifstream file1(in_file);
    while (getline(file1, line)) {
      std::vector<int> int_vec;
      bool invalid_residue = false;
      for (int i = 0; i < line.length(); i++) {
         if (alphabet.find(line[i]) == -1) {
            std::cerr << "Invalid amino acid found in " << line << " at position "
                    << i + 1 <<". Skipping." << std::endl;
            invalid_residue = true;
            break;   
        //   return EXIT_FAILURE;

        }
      }

      if(invalid_residue)
         continue;

      if( trimming ) 
        // std::cout <<line <<'\t';
        line = trim( line ); // this will only work for text file input, not for AIRR data. Must be processed afterwards.
        // std::cout <<line <<std::endl;

      peplist1.push_back({line, int(line.length()), -99.9, int_vec});
    }
    file1.close();
  }

// Calculate the normalization score (aff) (kernel 3 self vs self) list 1
#pragma omp parallel for
  for (int i = 0; i < peplist1.size(); i++) {
    peptide *pep_ptr = &peplist1[i];
    for (int x = 0; x < pep_ptr->len; x++) {
      pep_ptr->i.push_back(alphabet.find(pep_ptr->seq[x]));
    }
    pep_ptr->aff = k3_sum(*pep_ptr, *pep_ptr);
  }

  // change to IEDB data
  for (std::vector<std::string>::iterator it = iedb_data.begin();
       it != iedb_data.end(); it++) {
    std::vector<int> int_vec;
    for (int i = 0; i < (*it).length(); i++) {
      if (alphabet.find((*it)[i]) == -1) {
        std::cerr << "Invalid amino acid found in " << *it << " at position "
                  << i + 1 << std::endl;
        return EXIT_FAILURE;
      }
    }
    peplist2.push_back({*it, int((*it).length()), -99.9, int_vec});
  }

// Calculate the normalization score (aff) (kernel 3 self vs self) for list 2
#pragma omp parallel for
  for (int i = 0; i < peplist2.size(); i++) {
    peptide *pep_ptr = &peplist2[i];
    for (int x = 0; x < pep_ptr->len; x++) {
      pep_ptr->i.push_back(alphabet.find(pep_ptr->seq[x]));
    }
    pep_ptr->aff = k3_sum(*pep_ptr, *pep_ptr);
  }

  // Calculate input data vs database with multi-threading
  std::vector<std::string> extended_output;
  std::vector<std::string> tcrmatch_output = multi_calc_k3(peplist1, peplist2, threshold, iedb_map);

  std::string tcrmatch_output_header;
  std::string was_trimmed;
  if( trimming )
    was_trimmed = "trimmed_";
  else 
    was_trimmed = "";

  tcrmatch_output_header = was_trimmed + "input_sequence\tmatch_sequence\tscore\treceptor_group\tepitope\tantigen\torganism" ;

  if( trust4flag == true ){
    std::cout << trust4header <<"\t";
    std::cout << tcrmatch_output_header << std::endl; 

    // Appending TCRmatch to Trust4 output
    for( int i = 0; i < inputlines.size(); i++)
      for( int j = 0; j < tcrmatch_output.size(); j++) {
        std::string final_line = inputlines[i] + "\t" + tcrmatch_output[j];
        std::string AAseq_trust  = get_nth_field(inputlines[i], 3);
        std::string seq_tcrmatch = get_nth_field(tcrmatch_output[j], 0);
        // std::cout << (trim(AAseq_trust) == seq_tcrmatch) <<std::endl;
        if( trim(AAseq_trust) == seq_tcrmatch )
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
