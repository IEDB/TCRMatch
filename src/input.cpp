#include <stdexcept>

#include "input.hpp"
#include "valid_residues.hpp"
#include <sstream>
#include <iostream>

using namespace std;

// --------------------------------------------------------------------
void erase_newline( string& str )
{
	// to be used with trust4 inputs, airr inputs, etc
	if (str.back() == '\n' || str.back() == '\r')
		str.erase(str.find_last_of("\r\n"));
}


// --------------------------------------------------------------------
void _assert_valid_aminoacids( string line, int line_count, string input_file_name )
{
	char check_invalid = contain_invalid_residue( line );

	if( check_invalid )
		throw invalid_argument( string("Non-standard aminoacid '") + check_invalid + string("' found in line ") + to_string(line_count) + " of '" + string(input_file_name) +"' input file. Aborting.");

}

// --------------------------------------------------------------------
void parse_raw_input( vector<peptide>& peplist, string const& in_file, bool trimming )
{
	string line;
	int line_count = 1;

	ifstream input_file( in_file );

	while( getline(input_file, line) ) 
	{

		vector<int> int_vec; // this is not doing anything here

		_assert_valid_aminoacids( line, line_count++, in_file );

		peptide p_aux = {line, int(line.length()), -99.9, int_vec};
		
		if( trimming ) // maybe take trimming out of here?
			trim_flanking_residues(p_aux);

		peplist.push_back( p_aux );

	}
	
	input_file.close();

}

// ----------------------------------
vector<string> read_AIRR_data(string AIRR_data_file) 
{
  vector<string> airr_data;
  ifstream airr_file(AIRR_data_file);
  string line;
  string temp;
  int column = -1, idx;

  // get cdr3_aa column
  idx = 0;
  getline(airr_file, line);
  stringstream buffer(line);
  while (getline(buffer, temp, '\t')) {
    if (temp == "cdr3_aa") {
      column = idx;
      break;
    }
    idx++;
  }

  if (column < 0) 
		throw invalid_argument( "Cannot find cdr3_aa column in AIRR TSV input file. Aborting.");


  // read values
  while (getline(airr_file, line)) {
    stringstream buffer(line);
    vector<string> values;
    while (getline(buffer, temp, '\t')) {
      values.push_back(temp);
    }
    airr_data.push_back(values[column]);
  }
  return airr_data;
}
// ---------------------------------------------------------------------
void parse_airr_input( vector<peptide>& peplist, string const& in_file, bool trimming )
{
	// TODO: it's not trimming yet!
	vector<string> airr_data = read_AIRR_data(in_file);

	vector<int> int_vec;
	for (vector<string>::iterator it = airr_data.begin(); it != airr_data.end(); it++) 
		peplist.push_back({*it, int((*it).length()), -99.9, int_vec});
	
} 

// ---------------------------------------------------------------------
void parse_trust4_input( vector<peptide>& peplist, string const& in_file, vector<string>& inputlines, bool trimming )
{
	std::string trust4header;
    std::ifstream file1(in_file);

    // Checking if TRUST4 file is correctly formatted
    getline(file1, trust4header); // parses header
    std::string valid_header = "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\tcid_full_length";
    if ( trust4header != valid_header )
		throw invalid_argument( "ERROR: Invalid TRUST4 header. Please check input file format. Aborting.");
	
	std::string line, seq;
    
	while ( getline(file1, line) ) 
	{
		inputlines.push_back( line );

		std::istringstream is( line );
		std::string skip;
		is >>skip >>skip >>skip >>seq;  // 4th field contains peptide seq
		std::vector<int> int_vec;
		bool invalid_char_found = false;

		for (int i = 0; i < seq.length(); i++) 
		{
			if (valid_residues.find(seq[i]) == -1) {
			// return EXIT_FAILURE; // TRUST4 contains many illegal chars. Skipping them for now.
			invalid_char_found = true; 
			break;
			}
		}

		if( invalid_char_found ) 
		{
			cerr << "Invalid amino acid found in " << seq <<endl;
			cerr << "Skipping" <<endl;
			continue;
		}

		peptide p_aux = {seq, int(seq.length()), -99.9, int_vec};
		if( trimming ) // removes flanking residues (C and F or W). Works for text input, not AIRR.
			trim_flanking_residues( p_aux );

		peplist.push_back( p_aux );
    }

    file1.close();
}

// -------------------------------------------------------
// int main()
// {
// 	parse_raw_input( "testfiles/iedb/10seqs.test" );

// }