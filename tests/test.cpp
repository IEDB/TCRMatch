#include "../src/tcrmatch.cpp"
#include "catch_amalgamated.hpp"

TEST_CASE("Test k3_sum", "[k3_sum]") {
  std::string alphabet = "ARNDCQEGHILKMFPSTWYV";
  std::vector<int> int_vec_one;
  std::vector<int> int_vec_two;
  std::string test_seq_one = "ASSSRSSYEQY";
  std::string test_seq_two = "ASSNRSSYEQY";

  fmatrix_k1();

  for (int x = 0; x < test_seq_one.length(); x++) {
    int_vec_one.push_back(alphabet.find(test_seq_one[x]));
  }

  for (int x = 0; x < test_seq_two.length(); x++) {
    int_vec_two.push_back(alphabet.find(test_seq_two[x]));
  }

  peptide test_peptide_one = {test_seq_one, int(test_seq_one.length()), -99.9,
                              int_vec_one};
  peptide test_peptide_two = {test_seq_two, int(test_seq_two.length()), -99.9,
                              int_vec_two};

  REQUIRE_THAT(k3_sum(test_peptide_one, test_peptide_two), Catch::Matchers::WithinRelMatcher(620.97955f, 1e-5));
}

TEST_CASE("Test Trim", "[trim]") {
  std::string sequence = "CASSSRSSYEQF";

  REQUIRE(trim(sequence) == "ASSSRSSYEQ");
  REQUIRE(trim(sequence) != "CASSSRSSYEQF");
}

TEST_CASE("Test Integration", "[integration]") {
  std::string alphabet = "ARNDCQEGHILKMFPSTWYV";
  std::vector<int> int_vec_one;
  std::vector<int> int_vec_two;
  std::string test_seq_one = "ASSSRSSYEQY";
  std::string test_seq_two = "ASSNRSSYEQY";

  fmatrix_k1();

  for (int x = 0; x < test_seq_one.length(); x++) {
    int_vec_one.push_back(alphabet.find(test_seq_one[x]));
  }

  for (int x = 0; x < test_seq_two.length(); x++) {
    int_vec_two.push_back(alphabet.find(test_seq_two[x]));
  }

  peptide test_peptide_one = {test_seq_one, int(test_seq_one.length()), -99.9,
                              int_vec_one};
  peptide test_peptide_two = {test_seq_two, int(test_seq_two.length()), -99.9,
                              int_vec_two};
  test_peptide_one.aff = k3_sum(test_peptide_one, test_peptide_one);
  test_peptide_two.aff = k3_sum(test_peptide_two, test_peptide_two);

  float score = k3_sum(test_peptide_one, test_peptide_two) /
                sqrt(test_peptide_one.aff * test_peptide_two.aff);

  REQUIRE_THAT(score, Catch::Matchers::WithinRelMatcher(0.97332f, 1e-5));
}