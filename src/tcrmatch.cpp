#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <math.h>
#include <omp.h>
#include <string>
#include <tuple>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <unordered_set>
#include <atomic>

std::atomic<size_t> current_memory_usage(0); // Global memory tracker
size_t MEMORY_LIMIT = 2L * 1024 * 1024 * 1024; // max: 2GB

struct peptide {
  std::string seq;
  int len;
  float aff;
  std::vector<int> i;
};


struct TupleHash {
  template <typename T>
  std::size_t operator()(const T& tuple) const {
    auto hash1 = std::hash<std::string>()(std::get<0>(tuple));
    auto hash2 = std::hash<std::string>()(std::get<1>(tuple));
    auto hash3 = std::hash<float>()(std::get<2>(tuple));
    auto hash4 = std::hash<int>()(std::get<3>(tuple));
    return hash1 ^ hash2 ^ hash3 ^ hash4;
  }
};


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
  
	if (!str.empty() && (str.back() == '\n' || str.back() == '\r'))
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
   * A database can have several rows for the same trimmed sequence (one per
   * receptor group in the schema2 file); each sequence is returned once, in
   * order of first appearance, so it is only scored once.
   *
   * @param IEDB_data_file: path to IEDB data file
   * @return: vector of distinct trimmed sequences
   */
  std::vector<std::string> iedb_data;
  std::unordered_set<std::string> seen;
  std::ifstream iedb_file(IEDB_data_file);
  std::string line;
  while (getline(iedb_file, line)) {
    std::stringstream ss(line);
    std::string sequence;
    ss >> sequence;
    if (!sequence.empty() && sequence != "trimmed_seq" &&
        seen.insert(sequence).second) {
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
  std::string receptor_species;  // empty for databases without this column
};

std::vector<std::string> split_tabs(const std::string &line) {
  /**
   * Split a line on tabs. A trailing empty field (line ending in a tab) is
   * not returned; callers treat missing fields as empty.
   */
  std::vector<std::string> fields;
  std::stringstream buffer(line);
  std::string field;
  while (getline(buffer, field, '\t'))
    fields.push_back(field);
  return fields;
}

std::string field_or_empty(const std::vector<std::string> &fields, int idx) {
  if (idx < 0 || idx >= (int)fields.size())
    return "";
  return fields[idx];
}

std::map<std::string, std::vector<IEDB_data_row>>
create_IEDB_map(std::string IEDB_data_file, bool *is_schema2 = nullptr) {
  /**
   * Create a map of IEDB data, keyed by trimmed sequence
   *
   * If the first line starts with "trimmed_seq" it is a header and columns
   * are found by name, so both the legacy IEDB_data.tsv (6 columns) and
   * IEDB_data.schema2.tsv (adds receptor_species) work, and unknown extra
   * columns are ignored. Without a header (e.g. a plain list of CDR3b
   * sequences), columns are read positionally in the legacy order and missing
   * ones are left empty.
   *
   * @param IEDB_data_file: path to IEDB data file
   * @param is_schema2: if given, set to whether the header has the schema2
   *                    receptor_species column ('|'-separated, aligned lists)
   * @return: map of IEDB data
   */
  std::map<std::string, std::vector<IEDB_data_row>> iedb_map;
  if (is_schema2)
    *is_schema2 = false;
  std::ifstream iedb_file(IEDB_data_file);
  std::string line;

  // Column indices in the legacy (headerless) order.
  int c_trimmed = 0, c_original = 1, c_group = 2, c_epitopes = 3,
      c_organisms = 4, c_antigens = 5, c_species = -1;

  bool first_line = true;
  while (getline(iedb_file, line)) {
    erase_newline(line);
    std::vector<std::string> fields = split_tabs(line);

    if (first_line) {
      first_line = false;
      if (!fields.empty() && fields[0] == "trimmed_seq") {
        std::map<std::string, int> col;
        for (int i = 0; i < (int)fields.size(); i++)
          col[fields[i]] = i;
        const char *required[] = {"trimmed_seq",      "original_seq",
                                  "receptor_group",   "epitopes",
                                  "source_organisms", "source_antigens"};
        for (const char *name : required) {
          if (col.find(name) == col.end()) {
            std::cerr << "Database file " << IEDB_data_file
                      << " is missing the required column '" << name << "'"
                      << std::endl;
            exit(EXIT_FAILURE);
          }
        }
        c_trimmed = col["trimmed_seq"];
        c_original = col["original_seq"];
        c_group = col["receptor_group"];
        c_epitopes = col["epitopes"];
        c_organisms = col["source_organisms"];
        c_antigens = col["source_antigens"];
        c_species =
            col.count("receptor_species") ? col["receptor_species"] : -1;
        if (is_schema2)
          *is_schema2 = c_species >= 0;
        continue;
      }
    }

    std::string key = field_or_empty(fields, c_trimmed);
    if (key.empty())
      continue;
    iedb_map[key].push_back({key, field_or_empty(fields, c_original),
                             field_or_empty(fields, c_group),
                             field_or_empty(fields, c_epitopes),
                             field_or_empty(fields, c_organisms),
                             field_or_empty(fields, c_antigens),
                             field_or_empty(fields, c_species)});
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


std::vector<std::tuple<std::string, std::string, float, int>> remove_duplicate_matches(
    std::vector<std::tuple<std::string, std::string, float, int>>* results) 
    
{
  std::unordered_set<std::tuple<std::string, std::string, float, int>, TupleHash> seen; 
  std::vector<std::tuple<std::string, std::string, float, int>> unique_vec;
  for (int tid = 0; tid < omp_get_max_threads(); ++tid) 
    for (const auto& item : results[tid]) {
      if (seen.insert(item).second) { // Only insert if it's unique
        unique_vec.push_back(item);
    }
  }
  delete[] results; // Ensure to free memory of results

  return unique_vec;
}




std::vector<std::tuple<std::string, std::string, float, int>>* 
find_matches(std::vector<peptide>& peplist1, std::vector<peptide>& peplist2, float threshold)
{
  auto *results = new std::vector<std::tuple<std::string, std::string, float, int>>[omp_get_max_threads()];
  int max_matches = 1000000;
  int nresults = 0;

  #pragma omp parallel for
    for (int i = 0; i < peplist1.size(); i++) {
      for (int j = 0; j < peplist2.size(); j++) {
        peptide pep1 = peplist1[i];
        peptide pep2 = peplist2[j];
        float score = 0.0;
        score = k3_sum(pep1, pep2) / sqrt(pep1.aff * pep2.aff);
        if (score > threshold) {
          int tid = omp_get_thread_num();
          auto item = make_tuple(pep1.seq, pep2.seq, score, i);
          results[tid].push_back(item);
          nresults += 1;
          // Increment memory usage
          current_memory_usage += sizeof(item);
          
          if ( current_memory_usage > MEMORY_LIMIT){
            std::cerr << nresults << " matches processed so far. TCRMatch exceeded the " 
                      << MEMORY_LIMIT / (1024 * 1024 * 1024) << "GB limit for RAM usage." << std::endl;
            std::cerr << "Consider partitioning your input or increasing the memory limit using the -m flag." << std::endl;
            std::cerr << "Aborting." << std::endl;
			exit(0);
          }
        }
      }
    }
    return results;
}



std::vector<std::string> split_list(const std::string &list, char sep) {
  /**
   * Split a database list field, keeping empty entries (schema2 uses empty
   * placeholders so that entry n of each list describes the same epitope).
   */
  std::vector<std::string> items;
  size_t start = 0, pos;
  while ((pos = list.find(sep, start)) != std::string::npos) {
    items.push_back(list.substr(start, pos - start));
    start = pos + 1;
  }
  items.push_back(list.substr(start));
  return items;
}

std::string join_list(const std::vector<std::string> &items, char sep) {
  std::string out;
  for (size_t i = 0; i < items.size(); i++) {
    if (i)
      out += sep;
    out += items[i];
  }
  return out;
}

struct epitope_summary {
  std::string epitope;
  std::vector<std::string> receptor_groups;  // distinct, in database order
  std::unordered_set<std::string> seen_groups;
  // (antigen, organism) labels and how many receptor groups gave each
  std::vector<std::pair<std::pair<std::string, std::string>, int>> labels;
  std::vector<std::string> species;  // distinct non-empty, in database order
};

std::vector<std::string> collapse_by_epitope(
    const std::vector<std::tuple<std::string, std::string, float, int>> &matches,
    std::map<std::string, std::vector<IEDB_data_row>> &iedb_map) {
  /**
   * Collapse matches to one line per (input sequence, match sequence,
   * epitope), for the schema2 database. Each receptor group of the matched
   * sequence counts once for every epitope in its list.
   *
   * Output columns: input_sequence, match_sequence, score, epitope, antigen,
   * organism, n_receptor_groups, receptor_groups, n_receptor_species,
   * receptor_species.
   *   - antigen/organism: the label most receptor groups give the epitope,
   *     preferring a non-empty one (ties: first seen)
   *   - receptor_groups, receptor_species: '|'-separated, distinct
   * Within a match, epitopes are ordered by n_receptor_groups (descending),
   * then by epitope. Receptor groups without epitopes are left out.
   */
  std::vector<std::string> out;
  std::set<std::pair<std::string, std::string>> done;
  for (const auto &match : matches) {
    const std::string &input = std::get<0>(match);
    const std::string &match_seq = std::get<1>(match);
    // Repeated input lines give repeated matches; collapse them too.
    if (!done.insert(std::make_pair(input, match_seq)).second)
      continue;
    auto it = iedb_map.find(match_seq);
    if (it == iedb_map.end())
      continue;

    std::vector<epitope_summary> summaries;
    std::map<std::string, size_t> index;
    for (const IEDB_data_row &row : it->second) {
      std::vector<std::string> epitopes = split_list(row.epitope, '|');
      std::vector<std::string> antigens = split_list(row.antigen, '|');
      std::vector<std::string> organisms = split_list(row.organism, '|');
      for (size_t e = 0; e < epitopes.size(); e++) {
        if (epitopes[e].empty())
          continue;
        auto found = index.find(epitopes[e]);
        if (found == index.end()) {
          found = index.insert(std::make_pair(epitopes[e], summaries.size())).first;
          summaries.push_back(epitope_summary());
          summaries.back().epitope = epitopes[e];
        }
        epitope_summary &s = summaries[found->second];
        // The same group can appear on several rows (one per original
        // sequence); count it once.
        if (!s.seen_groups.insert(row.receptor_group).second)
          continue;
        s.receptor_groups.push_back(row.receptor_group);
        std::pair<std::string, std::string> label(
            field_or_empty(antigens, e), field_or_empty(organisms, e));
        bool counted = false;
        for (auto &l : s.labels)
          if (l.first == label) {
            l.second++;
            counted = true;
          }
        if (!counted)
          s.labels.push_back(std::make_pair(label, 1));
        if (!row.receptor_species.empty() &&
            std::find(s.species.begin(), s.species.end(),
                      row.receptor_species) == s.species.end())
          s.species.push_back(row.receptor_species);
      }
    }

    std::stable_sort(summaries.begin(), summaries.end(),
                     [](const epitope_summary &a, const epitope_summary &b) {
                       if (a.receptor_groups.size() != b.receptor_groups.size())
                         return a.receptor_groups.size() > b.receptor_groups.size();
                       return a.epitope < b.epitope;
                     });

    std::string prefix = input + "\t" + match_seq + "\t" +
                         std::to_string(std::get<2>(match)) + "\t";
    for (const epitope_summary &s : summaries) {
      const std::pair<std::string, std::string> *best = nullptr;
      int best_count = 0;
      bool best_empty = true;
      for (const auto &l : s.labels) {
        bool empty = l.first.first.empty() && l.first.second.empty();
        if (!best || (best_empty && !empty) ||
            (best_empty == empty && l.second > best_count)) {
          best = &l.first;
          best_count = l.second;
          best_empty = empty;
        }
      }
      out.push_back(prefix + s.epitope + "\t" + best->first + "\t" +
                    best->second + "\t" +
                    std::to_string(s.receptor_groups.size()) + "\t" +
                    join_list(s.receptor_groups, '|') + "\t" +
                    std::to_string(s.species.size()) + "\t" +
                    join_list(s.species, '|') + "\n");
    }
  }
  return out;
}

bool open_summary_file(char flag, const std::string &path, bool is_schema2,
                       std::ofstream &out) {
  /**
   * Open the file for an extra output view (-p, -c), if one was requested.
   * Both views split the database's lists by receptor group, so they need a
   * schema2 database. Called before scoring, so problems are reported early.
   *
   * @return: false (after printing an error) if the run should stop
   */
  if (path.empty())
    return true;
  if (!is_schema2) {
    std::cerr << "-" << flag
              << " needs a schema2 database (with a receptor_species column), "
                 "e.g. -d data/IEDB_data.schema2.tsv. The legacy database "
                 "pools epitopes across receptor groups, so they cannot be "
                 "split by receptor group."
              << std::endl;
    return false;
  }
  out.open(path);
  if (!out) {
    std::cerr << "Cannot write to " << path << std::endl;
    return false;
  }
  return true;
}

bool write_summary_file(std::ofstream &out, const std::string &path,
                        const std::string &header,
                        const std::vector<std::string> &lines) {
  /** Write an extra output view opened by open_summary_file, if any. */
  if (!out.is_open())
    return true;
  out << header << "\n";
  for (const std::string &line : lines)
    out << line;
  out.close();
  if (!out) {
    std::cerr << "Error writing " << path << std::endl;
    return false;
  }
  return true;
}

std::string without_semicolons(std::string text) {
  // ';' separates receptor groups in the standard output; a few names
  // contain one (e.g. "UL123; IE1"), so write it as ',' there.
  std::replace(text.begin(), text.end(), ';', ',');
  return text;
}

std::string count_distinct(const std::vector<std::string> &values) {
  // The n_* columns: how many distinct non-empty values a list holds.
  std::unordered_set<std::string> seen;
  for (const std::string &value : values)
    if (!value.empty())
      seen.insert(value);
  return std::to_string(seen.size());
}

void append_split(std::vector<std::string> &values, const std::string &list) {
  // Add the entries of a schema2 '|' list to values.
  for (const std::string &item : split_list(list, '|'))
    values.push_back(item);
}

typedef std::tuple<std::string, std::string, float, int> match_tuple;

std::vector<const match_tuple *>
distinct_matches(const std::vector<match_tuple> &matches) {
  // Each (input sequence, match sequence) pair once, in order of first
  // appearance: an input sequence given more than once is reported once.
  std::vector<const match_tuple *> out;
  std::set<std::pair<std::string, std::string>> done;
  for (const match_tuple &match : matches)
    if (done.insert(std::make_pair(std::get<0>(match), std::get<1>(match))).second)
      out.push_back(&match);
  return out;
}

std::vector<IEDB_data_row>
rows_for_match(const std::string &match_seq,
               const std::map<std::string, std::vector<IEDB_data_row>> &iedb_map) {
  // The database rows for a matched sequence; a sequence without metadata
  // still gets one (unannotated) row.
  auto it = iedb_map.find(match_seq);
  if (it != iedb_map.end())
    return it->second;
  return {{match_seq, "", "", "", "", "", ""}};
}

std::string match_prefix(const match_tuple &match) {
  return std::get<0>(match) + "\t" + std::get<1>(match) + "\t" +
         std::to_string(std::get<2>(match)) + "\t";
}

std::vector<std::string>
group_by_cdr3(const std::vector<match_tuple> &matches,
              const std::map<std::string, std::vector<IEDB_data_row>> &iedb_map,
              bool is_schema2) {
  /**
   * The standard output: one line per (input sequence, match sequence,
   * untrimmed CDR3b as recorded in IEDB), the same lines as the legacy
   * database's rows.
   *
   * Columns: input_sequence, match_sequence, score, untrimmed_match_sequence,
   * then n_receptor_groups, receptor_groups, n_epitopes, epitopes,
   * n_antigens, antigens, n_organisms, organisms, n_receptor_species,
   * receptor_species. Each list column has one ';'-separated entry per
   * receptor group, in the same order; with schema2, each entry is that
   * group's '|'-separated list, unchanged. Each n_* column counts the
   * distinct non-empty values in its list.
   *
   * A legacy or headerless database has one receptor group per line and its
   * lists cannot be trusted, so the n_* columns are left empty.
   */
  std::vector<std::string> out;
  for (const match_tuple *match : distinct_matches(matches)) {
    std::vector<IEDB_data_row> rows = rows_for_match(std::get<1>(*match), iedb_map);

    // One block of rows per output line, in database order: schema2 rows
    // are grouped by untrimmed sequence. Legacy and headerless databases
    // already have one row per untrimmed sequence (or none recorded), so
    // each of their rows is its own line, as in earlier versions.
    std::vector<std::vector<const IEDB_data_row *>> blocks;
    std::map<std::string, size_t> block_of_original;
    for (const IEDB_data_row &row : rows) {
      auto found = block_of_original.find(row.original_sequence);
      if (is_schema2 && found != block_of_original.end()) {
        blocks[found->second].push_back(&row);
      } else {
        block_of_original[row.original_sequence] = blocks.size();
        blocks.push_back({&row});
      }
    }

    for (const auto &block : blocks) {
      const std::string &original = block.front()->original_sequence;
      std::vector<std::string> groups, epitopes, antigens, organisms, species;
      std::vector<std::string> all_epitopes, all_antigens, all_organisms;
      for (const IEDB_data_row *row : block) {
        groups.push_back(row->receptor_group);
        epitopes.push_back(without_semicolons(row->epitope));
        antigens.push_back(without_semicolons(row->antigen));
        organisms.push_back(without_semicolons(row->organism));
        species.push_back(without_semicolons(row->receptor_species));
        append_split(all_epitopes, row->epitope);
        append_split(all_antigens, row->antigen);
        append_split(all_organisms, row->organism);
      }
      auto count = [is_schema2](const std::vector<std::string> &values) {
        return is_schema2 ? count_distinct(values) : std::string();
      };
      out.push_back(match_prefix(*match) + original + "\t" +
                    count(groups) + "\t" + join_list(groups, ';') + "\t" +
                    count(all_epitopes) + "\t" + join_list(epitopes, ';') + "\t" +
                    count(all_antigens) + "\t" + join_list(antigens, ';') + "\t" +
                    count(all_organisms) + "\t" + join_list(organisms, ';') + "\t" +
                    count(species) + "\t" + join_list(species, ';') + "\n");
    }
  }
  return out;
}

std::vector<std::string> per_receptor_group(
    const std::vector<match_tuple> &matches,
    const std::map<std::string, std::vector<IEDB_data_row>> &iedb_map) {
  /**
   * The -p view, for the schema2 database: one line per (input sequence,
   * match sequence, receptor group).
   *
   * Columns: input_sequence, match_sequence, score, untrimmed_match_sequence,
   * receptor_group, n_epitopes, epitopes, n_antigens, antigens, n_organisms,
   * organisms, receptor_species. The lists are the group's '|'-separated
   * lists from the database; each n_* column counts the distinct non-empty
   * values in its list.
   */
  std::vector<std::string> out;
  for (const match_tuple *match : distinct_matches(matches)) {
    for (const IEDB_data_row &row : rows_for_match(std::get<1>(*match), iedb_map)) {
      out.push_back(match_prefix(*match) + row.original_sequence + "\t" +
                    row.receptor_group + "\t" +
                    count_distinct(split_list(row.epitope, '|')) + "\t" + row.epitope + "\t" +
                    count_distinct(split_list(row.antigen, '|')) + "\t" + row.antigen + "\t" +
                    count_distinct(split_list(row.organism, '|')) + "\t" + row.organism + "\t" +
                    row.receptor_species + "\n");
    }
  }
  return out;
}

std::vector<std::string>
multi_calc_k3(std::vector<peptide> peplist1, std::vector<peptide> peplist2,
              float threshold,
              std::map<std::string, std::vector<IEDB_data_row>> &iedb_map,
              bool is_schema2,
              std::vector<std::string> *per_group_output = nullptr,
              std::vector<std::string> *collapsed_output = nullptr) {
  /**
   * Score all peptides in peplist1 against peplist2 and build the output
   *
   * @param peplist1: peptide list 1 (input)
   * @param peplist2: peptide list 2 (database)
   * @param threshold: K3 threshold
   * @param iedb_map: IEDB map
   * @param is_schema2: whether the database is schema2 (see create_IEDB_map)
   * @param per_group_output: if given, also filled with one line per
   *                          receptor group (see per_receptor_group)
   * @param collapsed_output: if given, also filled with the matches
   *                          collapsed by epitope (see collapse_by_epitope)
   * @return: the standard output lines (see group_by_cdr3)
   */
  auto results = find_matches(peplist1, peplist2, threshold);
  auto unique_vec = remove_duplicate_matches(results);
  if (per_group_output)
    *per_group_output = per_receptor_group(unique_vec, iedb_map);
  if (collapsed_output)
    *collapsed_output = collapse_by_epitope(unique_vec, iedb_map);
  return group_by_cdr3(unique_vec, iedb_map, is_schema2);
}
