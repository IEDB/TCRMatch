#include <array>
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>

std::array<std::array<float, 20>, 20> parse_matrix( const std::string& inputfilename );
std::array<std::array<float, 20>, 20> normalize( const std::array<std::array<float, 20>, 20>& matrix, float norm_factor );
void print_matrix( const std::array<std::array<float, 20>, 20> &submatrix  );
