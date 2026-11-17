#define VERSION "1.4"

#include "tcrmatch.cpp"

const char *USAGE =
    "Usage: ./tcrmatch -i infile_name.txt -t num_threads [-s score_threshold] "
    "[-d /path/to/database] [-a | -r] [-k] [-m MAX_RAM_GB] "
    "[-p per_receptor_group.tsv] [-c collapsed_by_epitope.tsv] "
    "[--species \"Homo sapiens,...\"] "
    "[--antigen-type neoantigen,viral,germline,other]\n"
    "       ./tcrmatch -l [-d /path/to/database]   (list the receptor species)\n"
    "       ./tcrmatch -v                          (print the version)";

int main(int argc, char *argv[]) {
  int opt;
  int n_threads;
  float threshold;
  std::string in_file;
  std::string iedb_file = "data/IEDB_data.schema2.tsv";
  std::vector<std::string> species;  // --species: keep only these receptors
  std::vector<std::string> antigen_types;  // --antigen-type: keep only these epitopes
  bool list_species_flag = false;
  std::string per_group_file;  // -p: also write one line per receptor group
  std::string collapsed_file;  // -c: also write matches collapsed by epitope
  int i_flag = -1;
  int t_flag = -1;
  int thresh_flag = -1;
  int airr_flag = -1;
  bool trimming = true;
  bool trust4flag = false;
  // Command line argument parsing
  enum { OPT_SPECIES = 1000, OPT_ANTIGEN_TYPE };
  static struct option long_options[] = {
      {"list-species", no_argument, nullptr, 'l'},
      {"species", required_argument, nullptr, OPT_SPECIES},
      {"antigen-type", required_argument, nullptr, OPT_ANTIGEN_TYPE},
      {"version", no_argument, nullptr, 'v'},
      {nullptr, 0, nullptr, 0}};
  while ((opt = getopt_long(argc, argv, "rakvlt:i:s:d:m:c:p:", long_options,
                            nullptr)) != -1) {
    switch (opt) {
      case 'l':
        list_species_flag = true;
        break;
      case OPT_SPECIES:
        species = parse_name_list(optarg);
        if (species.empty()) {
          std::cerr << "--species needs at least one species name, e.g. "
                       "--species \"Homo sapiens\""
                    << std::endl;
          return EXIT_FAILURE;
        }
        break;
      case OPT_ANTIGEN_TYPE:
        antigen_types = parse_antigen_types(optarg);
        break;
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
      case 'c':
        collapsed_file = optarg;
        break;
      case 'p':
        per_group_file = optarg;
        break;
      case 'k':
        trimming = false;
        break;
      case 'r':
        trust4flag = true;
        break;
      case 'm':
        MEMORY_LIMIT = static_cast<size_t>(std::stol(optarg)) * 1024 * 1024 * 1024; // Convert GB to bytes
        break;
      case 'v':
        std::cout << VERSION << std::endl;
        return 0;
      default:
        std::cerr << USAGE << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (list_species_flag)
    return list_species(iedb_file);

  // Check that required parameters are there + error correcting
  if (i_flag == -1 || t_flag == -1) {
    std::cerr << "Missing mandatory parameters" << std::endl
              << USAGE << std::endl;
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

  // The distinct database sequences to score, in file order.
  std::vector<std::string> iedb_data;
  std::map<std::string, std::vector<IEDB_data_row>> iedb_map =
      create_IEDB_map(iedb_file, &iedb_data, species, antigen_types);
  if (iedb_data.empty()) {
    if (!antigen_types.empty())
      std::cerr << "No receptors in " << iedb_file
                << " match the requested antigen types"
                << (species.empty() ? "" : " and species") << "." << std::endl;
    else if (!species.empty())
      std::cerr << "No receptors of the requested species in " << iedb_file
                << ". Run with --list-species to see the species it has."
                << std::endl;
    else
      std::cerr << "No sequences found in the database file " << iedb_file
                << std::endl;
    return EXIT_FAILURE;
  }

  // Check -p and -c before scoring, so a bad path fails fast.
  std::ofstream per_group_out, collapsed_out;
  if (!open_summary_file(per_group_file, per_group_out) ||
      !open_summary_file(collapsed_file, collapsed_out))
    return EXIT_FAILURE;
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

  // k3_sum can't score sequences longer than MAX_SEQ_LEN.
  {
    std::vector<peptide> kept;
    for (const peptide &pep : peplist1) {
      if (pep.len > MAX_SEQ_LEN)
        std::cerr << "Input sequence " << pep.seq << " is longer than "
                  << MAX_SEQ_LEN << " residues. Skipping." << std::endl;
      else
        kept.push_back(pep);
    }
    peplist1.swap(kept);
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
    if ((int)(*it).length() > MAX_SEQ_LEN) {
      std::cerr << "Database sequence " << *it << " is longer than "
                << MAX_SEQ_LEN << " residues" << std::endl;
      return EXIT_FAILURE;
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
  // Standard output: one line per input, match and untrimmed CDR3b, with
  // ';' lists by receptor group (see group_by_cdr3).
  std::vector<std::string> per_group_output, collapsed_output;
  std::vector<std::string> tcrmatch_output = multi_calc_k3(
      peplist1, peplist2, threshold, iedb_map,
      per_group_file.empty() ? nullptr : &per_group_output,
      collapsed_file.empty() ? nullptr : &collapsed_output);

  std::string tcrmatch_output_header;
  std::string was_trimmed;
  if( trimming )
    was_trimmed = "trimmed_";
  else
    was_trimmed = "";

  // Extra views. Per receptor group (-p): one line per input, match and
  // receptor group. Collapsed (-c): one line per input, match and epitope.
  // In TRUST4 mode both are keyed on the CDR3 sequence, without the TRUST4
  // columns.
  if (!write_summary_file(
          per_group_out, per_group_file,
          was_trimmed + "input_sequence\tmatch_sequence\tscore\t"
                        "untrimmed_match_sequence\treceptor_group\t"
                        "n_epitopes\tepitopes\tn_antigens\tantigens\t"
                        "n_organisms\torganisms\treceptor_species",
          per_group_output) ||
      !write_summary_file(
          collapsed_out, collapsed_file,
          was_trimmed + "input_sequence\tmatch_sequence\tscore\tepitope\t"
                        "antigen\torganism\tn_receptor_groups\t"
                        "receptor_groups\tn_receptor_species\t"
                        "receptor_species",
          collapsed_output))
    return EXIT_FAILURE;

  tcrmatch_output_header =
      was_trimmed + "input_sequence\tmatch_sequence\tscore\t"
                    "untrimmed_match_sequence\tn_receptor_groups\t"
                    "receptor_groups\tn_epitopes\tepitopes\tn_antigens\t"
                    "antigens\tn_organisms\torganisms\tn_receptor_species\t"
                    "receptor_species";

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
