#include "../src/tcrmatch.cpp"
#include "catch_amalgamated.hpp"

TEST_CASE("Test int_vec", "[k3_sum]") {
    
}

TEST_CASE("Test k3_sum", "[k3_sum]") {
    std::string alphabet = "ARNDCQEGHILKMFPSTWYV";
    std::vector<int> int_vec;
    std::string test_seq_one = "AAAAA";
    
    for (int x = 0; x < test_seq_one.length(); x++) {
      int_vec.push_back(alphabet.find(test_seq_one[x]));
    }
    
    peptide test_peptide_one = {test_seq_one, int(test_seq_one.length()), -99.9, int_vec};
    
    REQUIRE(k3_sum(test_peptide_one, test_peptide_one) == 1);
}