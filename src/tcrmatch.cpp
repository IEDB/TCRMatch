#include <algorithm>
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

bool is_TCR_gene( std::string& str )
{
  /**
   * Checks if any of the VDJ fields
   * of the TRUST4 input file is a
   * gene for TCR (TR)
   *
   * @param str: line from input file
   * @return bool: true if is a gene for immunoglobulin
   */  

  std::istringstream iss( str );
  std::string field;
  
  for (int i = 0; i <= 6; i++)
  {
    iss >> field;

    if( i == 4 or i == 5 or i == 6)
      if( field.substr(0, 3) == "TRB" )
        return true;
  }

  return false;
}



void erase_newline( std::string& str ) 
{
  /**
   * Removes newline characters (\rn \n)
   * to avoid conflicts of inputs from
   * windows vs unix systems
   *
   * @param str: line from input file
   */  
  
	if (str.back() == '\n' || str.back() == '\r') 
		str.erase(str.find_last_of("\r\n"));
}




std::string get_nth_field(std::string &str, int n) {
  /**
   * Get the nth field of a string
   *
   * @param str: string to parse
   * @param n: field number to return
   * @return: nth field of the string
   */
  std::istringstream iss(str);
  std::string field;
  for (int i = 0; i <= n; i++)
    iss >> field;
  return field;

}

std::array<std::array<float, 20>, 20> k1;
int P_KMIN = 1;
int P_KMAX = 30;
float P_BETA = 0.11387;

// Hardcoded because parsing + computing matrix is annoying
float blm_qij[20][20] = {
    {0.0215, 0.0023, 0.0019, 0.0022, 0.0016, 0.0019, 0.003,
     0.0058, 0.0011, 0.0032, 0.0044, 0.0033, 0.0013, 0.0016,
     0.0022, 0.0063, 0.0037, 0.0004, 0.0013, 0.0051},
    {0.0023, 0.0178, 0.002,  0.0016, 0.0004, 0.0025, 0.0027,
     0.0017, 0.0012, 0.0012, 0.0024, 0.0062, 0.0008, 0.0009,
     0.001,  0.0023, 0.0018, 0.0003, 0.0009, 0.0016},
    {0.0019, 0.002,  0.0141, 0.0037, 0.0004, 0.0015, 0.0022,
     0.0029, 0.0014, 0.001,  0.0014, 0.0024, 0.0005, 0.0008,
     0.0009, 0.0031, 0.0022, 0.0002, 0.0007, 0.0012},
    {0.0022, 0.0016, 0.0037, 0.0213, 0.0004, 0.0016, 0.0049,
     0.0025, 0.001,  0.0012, 0.0015, 0.0024, 0.0005, 0.0008,
     0.0012, 0.0028, 0.0019, 0.0002, 0.0006, 0.0013},
    {0.0016, 0.0004, 0.0004, 0.0004, 0.0119, 0.0003, 0.0004,
     0.0008, 0.0002, 0.0011, 0.0016, 0.0005, 0.0004, 0.0005,
     0.0004, 0.001,  0.0009, 0.0001, 0.0003, 0.0014},
    {0.0019, 0.0025, 0.0015, 0.0016, 0.0003, 0.0073, 0.0035,
     0.0014, 0.001,  0.0009, 0.0016, 0.0031, 0.0007, 0.0005,
     0.0008, 0.0019, 0.0014, 0.0002, 0.0007, 0.0012},
    {0.003,  0.0027, 0.0022, 0.0049, 0.0004, 0.0035, 0.0161,
     0.0019, 0.0014, 0.0012, 0.002,  0.0041, 0.0007, 0.0009,
     0.0014, 0.003,  0.002,  0.0003, 0.0009, 0.0017},
    {0.0058, 0.0017, 0.0029, 0.0025, 0.0008, 0.0014, 0.0019,
     0.0378, 0.001,  0.0014, 0.0021, 0.0025, 0.0007, 0.0012,
     0.0014, 0.0038, 0.0022, 0.0004, 0.0008, 0.0018},
    {0.0011, 0.0012, 0.0014, 0.001,  0.0002, 0.001,  0.0014,
     0.001,  0.0093, 0.0006, 0.001,  0.0012, 0.0004, 0.0008,
     0.0005, 0.0011, 0.0007, 0.0002, 0.0015, 0.0006},
    {0.0032, 0.0012, 0.001,  0.0012, 0.0011, 0.0009, 0.0012,
     0.0014, 0.0006, 0.0184, 0.0114, 0.0016, 0.0025, 0.003,
     0.001,  0.0017, 0.0027, 0.0004, 0.0014, 0.012},
    {0.0044, 0.0024, 0.0014, 0.0015, 0.0016, 0.0016, 0.002,
     0.0021, 0.001,  0.0114, 0.0371, 0.0025, 0.0049, 0.0054,
     0.0014, 0.0024, 0.0033, 0.0007, 0.0022, 0.0095},
    {0.0033, 0.0062, 0.0024, 0.0024, 0.0005, 0.0031, 0.0041,
     0.0025, 0.0012, 0.0016, 0.0025, 0.0161, 0.0009, 0.0009,
     0.0016, 0.0031, 0.0023, 0.0003, 0.001,  0.0019},
    {0.0013, 0.0008, 0.0005, 0.0005, 0.0004, 0.0007, 0.0007,
     0.0007, 0.0004, 0.0025, 0.0049, 0.0009, 0.004,  0.0012,
     0.0004, 0.0009, 0.001,  0.0002, 0.0006, 0.0023},
    {0.0016, 0.0009, 0.0008, 0.0008, 0.0005, 0.0005, 0.0009,
     0.0012, 0.0008, 0.003,  0.0054, 0.0009, 0.0012, 0.0183,
     0.0005, 0.0012, 0.0012, 0.0008, 0.0042, 0.0026},
    {0.0022, 0.001,  0.0009, 0.0012, 0.0004, 0.0008, 0.0014,
     0.0014, 0.0005, 0.001,  0.0014, 0.0016, 0.0004, 0.0005,
     0.0191, 0.0017, 0.0014, 0.0001, 0.0005, 0.0012},
    {0.0063, 0.0023, 0.0031, 0.0028, 0.001,  0.0019, 0.003,
     0.0038, 0.0011, 0.0017, 0.0024, 0.0031, 0.0009, 0.0012,
     0.0017, 0.0126, 0.0047, 0.0003, 0.001,  0.0024},
    {0.0037, 0.0018, 0.0022, 0.0019, 0.0009, 0.0014, 0.002,
     0.0022, 0.0007, 0.0027, 0.0033, 0.0023, 0.001,  0.0012,
     0.0014, 0.0047, 0.0125, 0.0003, 0.0009, 0.0036},
    {0.0004, 0.0003, 0.0002, 0.0002, 0.0001, 0.0002, 0.0003,
     0.0004, 0.0002, 0.0004, 0.0007, 0.0003, 0.0002, 0.0008,
     0.0001, 0.0003, 0.0003, 0.0065, 0.0009, 0.0004},
    {0.0013, 0.0009, 0.0007, 0.0006, 0.0003, 0.0007, 0.0009,
     0.0008, 0.0015, 0.0014, 0.0022, 0.001,  0.0006, 0.0042,
     0.0005, 0.001,  0.0009, 0.0009, 0.0102, 0.0015},
    {0.0051, 0.0016, 0.0012, 0.0013, 0.0014, 0.0012, 0.0017,
     0.0018, 0.0006, 0.012,  0.0095, 0.0019, 0.0023, 0.0026,
     0.0012, 0.0024, 0.0036, 0.0004, 0.0015, 0.0196}};

struct peptide {
  std::string seq;
  int len;
  float aff;
  std::vector<int> i;
};

std::string trim(std::string line) {
  /**
   * Trim the sequence from the IEDB data file
   *
   * @param line: line from IEDB data file
   * @return: trimmed sequence
   */
  int last_idx = line.size() - 1;

  if (line[0] == 'C' and (line[last_idx] == 'F' or line[last_idx] == 'W'))
    return line.substr(1, last_idx - 1);
  else
    return line;
}

std::vector<std::string> read_IEDB_data(std::string IEDB_data_file) {
  /**
   * Read the IEDB data file
   *
   * @param IEDB_data_file: path to IEDB data file
   * @return: vector of trimmed sequences
   */
  std::vector<std::string> iedb_data;
  std::ifstream iedb_file(IEDB_data_file);
  std::string line;
  while (getline(iedb_file, line)) {
    std::stringstream ss(line);
    std::string sequence;
    ss >> sequence;
    if (sequence != "trimmed_seq") {
      iedb_data.push_back(sequence);
    }
  }

  return iedb_data;
}

struct IEDB_data_row {
  std::string amino_acid_sequence;
  std::string original_sequence;
  std::string receptor_group;
  std::string epitope;
  std::string organism;
  std::string antigen;
};

std::map<std::string, std::vector<IEDB_data_row>>
create_IEDB_map(std::string IEDB_data_file) {
  /**
   * Create a map of IEDB data
   *
   * @param IEDB_data_file: path to IEDB data file
   * @return: map of IEDB data
   */
  std::map<std::string, std::vector<IEDB_data_row>> iedb_map;
  std::map<std::string, std::vector<IEDB_data_row>>::iterator it;
  std::vector<IEDB_data_row> iedb_data;
  std::ifstream iedb_file(IEDB_data_file);
  IEDB_data_row input_row;
  std::string line;
  std::string temp;
  int idx;

  // read values
  while (getline(iedb_file, line)) {
    idx = 0;
    std::stringstream buffer(line);
    std::string values[6];
    while (getline(buffer, temp, '\t')) {
      values[idx] = temp;
      idx++;
    }
    it = iedb_map.find(values[0]);
    if (it != iedb_map.end()) {
      iedb_map[values[0]].push_back(
          {values[0], values[1], values[2], values[3], values[4], values[5]});
    } else {
      iedb_map[values[0]] = {
          {values[0], values[1], values[2], values[3], values[4], values[5]}};
    }
  }

  return iedb_map;
}

std::vector<std::string> read_AIRR_data(std::string AIRR_data_file) {
  /**
   * Read the AIRR data file
   *
   * @param AIRR_data_file: path to AIRR data file
   * @return: vector of CDR3 sequences
   */
  std::vector<std::string> airr_data;
  std::ifstream airr_file(AIRR_data_file);
  std::string line;
  std::string temp;
  int column = -1, idx;

  // get cdr3_aa column
  idx = 0;
  std::getline(airr_file, line);
  std::stringstream buffer(line);
  while (getline(buffer, temp, '\t')) {
    if (temp == "cdr3_aa") {
      column = idx;
      break;
    }
    idx++;
  }

  if (column < 0) {
    std::cerr << "Cannot find cdr3_aa column in AIRR TSV input file"
              << std::endl;
    return airr_data;
  }

  // read values
  while (getline(airr_file, line)) {
    std::stringstream buffer(line);
    std::vector<std::string> values;
    while (getline(buffer, temp, '\t')) {
      values.push_back(temp);
    }
    airr_data.push_back(values[column]);
  }
  return airr_data;
}

std::array<std::array<float, 20>, 20> fmatrix_k1() {
  /**
   * Calculate K1 matrix
   *
   * @return: K1 matrix
   */
  int k, j;
  float marg[20];
  float sum;

  // initialize margin array
  for (int i = 0; i < 20; i++) {
    marg[i] = 0.0;
  }
  // initialize k1
  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 20; j++) {
      k1[i][j] = 0.0;
    }
  }
  // normalize matrix by marginal frequencies
  for (j = 0; j < 20; j++) {
    sum = 0;
    for (k = 0; k < 20; k++)
      sum += blm_qij[j][k];
    marg[j] = sum;
  }
  // calculate K1
  for (j = 0; j < 20; j++) {
    for (k = 0; k < 20; k++) {
      k1[j][k] = blm_qij[j][k] / (marg[j] * marg[k]);
      k1[j][k] = pow(k1[j][k], P_BETA);
    }
  }

  return (k1);
}

float k3_sum(peptide pep1, peptide pep2) {
  /**
   * Recursively calculate K3 sum
   *
   * @param pep1: peptide 1
   * @param pep2: peptide 2
   * @return: K3 sum
   */
  float k2, term, k3 = 0.0;
  int start1, start2;
  int k, j1, j2;

  float k2_prod_save[31][31][31];

  for (k = P_KMIN; k <= P_KMAX; k++) {
    for (start1 = 0; start1 <= pep1.len - k; start1++) {
      for (start2 = 0; start2 <= pep2.len - k; start2++) {
        j1 = pep1.i[start1 + k - 1];
        j2 = pep2.i[start2 + k - 1];
        term = k1[j1][j2];

        if (k == 1) {
          k2 = term;
        } else {
          k2 = k2_prod_save[k - 1][start1][start2] * term;
        }

        k2_prod_save[k][start1][start2] = k2;
        k3 += k2;
      }
    }
  }
  return (k3);
}

std::vector<std::string>
multi_calc_k3(std::vector<peptide> peplist1, std::vector<peptide> peplist2,
              float threshold,
              std::map<std::string, std::vector<IEDB_data_row>> iedb_map) {
  /**
   * Calculate K3 for all peptides in peplist1 and peplist2
   *
   * @param peplist1: peptide list 1
   * @param peplist2: peptide list 2
   * @param threshold: K3 threshold
   * @param iedb_map: IEDB map
   * @return: vector of matching CDR3 sequences
   */
  std::vector<std::tuple<std::string, std::string, float, int>>::iterator
      it2[omp_get_max_threads()];
  std::vector<std::tuple<std::string, std::string, float, int>>
      results[omp_get_max_threads()];

#pragma omp parallel for
  for (int i = 0; i < peplist1.size(); i++) {
    for (int j = 0; j < peplist2.size(); j++) {
      peptide pep1 = peplist1[i];
      peptide pep2 = peplist2[j];
      float score = 0.0;
      score = k3_sum(pep1, pep2) / sqrt(pep1.aff * pep2.aff);
      if (score > threshold) {
        int tid = omp_get_thread_num();
        // If input seq-match seq is unique, add tuple to results
        // Repeat IEDB matches are skipped (prevents duplicate rows in results)
        // but duplicate input sequences are permitted (e.g. for repertoires
        // with identical TCRs)
        it2[tid] = find(results[tid].begin(), results[tid].end(),
                        make_tuple(pep1.seq, pep2.seq, score, i));
        if (it2[tid] == results[tid].end()) {
          results[tid].push_back(make_tuple(pep1.seq, pep2.seq, score, i));
        }
      }
    }
  }
  //   std::cout << "input_sequence\tmatch_"
  //                "sequence\tscore\treceptor_"
  //                "group\tepitope\tantigen\torganism\t"
  //             << std::endl;

  std::vector<std::string> tcrmatch_output;

  for (int k = 0; k < omp_get_max_threads(); k++) {
    for (auto &tuple : results[k]) {
      std::string match_peptide = std::get<1>(tuple);
      std::map<std::string, std::vector<IEDB_data_row>>::iterator it =
          iedb_map.find(match_peptide);
      std::vector<IEDB_data_row> match_row_vec = it->second;

      for (int l = 0; l < match_row_vec.size(); l++) {
        IEDB_data_row match_row = match_row_vec[l];
        //   std::cout << std::fixed << std::setprecision(4) <<
        //   std::get<0>(tuple)
        //             << "\t" << std::get<1>(tuple) << "\t" <<
        //             std::get<2>(tuple)
        //             << "\t" << match_row.receptor_group << "\t"
        //             << match_row.epitope << "\t" << match_row.antigen << "\t"
        //             << match_row.organism << std::endl;
        std::string score = std::to_string(std::get<2>(tuple));
        std::string output_string =
            std::get<0>(tuple) + "\t" + std::get<1>(tuple) + "\t" + score +
            "\t" + match_row.receptor_group + "\t" + match_row.epitope + "\t" +
            match_row.antigen + "\t" + match_row.organism + "\n";

        tcrmatch_output.push_back(output_string);
      }
    }
  }
  return tcrmatch_output;
}
