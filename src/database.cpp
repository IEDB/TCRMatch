#include <fstream>
#include <sstream>

#include "database.hpp"

using namespace std;


// -------------------------------------------
vector<string> parse_database(string db_file) 
{
	vector<string> iedb_data;
	ifstream iedb_file(db_file);
	string line;

	getline( iedb_file, line ); // skip header

	while( getline(iedb_file, line) ) 
	{
		
		stringstream ss(line);
		string sequence;
		ss >> sequence;
		
		iedb_data.push_back(sequence);
		
	}

	return iedb_data;
}

// -------------------------------------------
map<string, vector<IEDB_data_row>> map_database( string IEDB_data_file ) 
{
	map<string, vector<IEDB_data_row>> iedb_map;
	map<string, vector<IEDB_data_row>>::iterator it;
	vector<IEDB_data_row> iedb_data;
	ifstream iedb_file(IEDB_data_file);
	IEDB_data_row input_row;
	string line;
	string temp;
	int idx;

	// read values
	while (getline(iedb_file, line)) 
	{
		idx = 0;
		stringstream buffer(line);
		string values[6];

		while (getline(buffer, temp, '\t')) 
		{
			values[idx] = temp;
			idx++;
		}

		it = iedb_map.find(values[0]);

		if (it != iedb_map.end()) 
			iedb_map[values[0]].push_back({values[0], values[1], values[2], values[3], values[4], values[5]});
		else 
			iedb_map[values[0]] = {{values[0], values[1], values[2], values[3], values[4], values[5]}};
		
	}

	return iedb_map;
}
