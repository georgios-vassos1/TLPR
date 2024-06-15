#include "utils.hpp" // Include the header file to ensure consistency

namespace utils {

  // Function to flatten a vector of vectors of equal lengths
  std::vector<double> flatten(const std::vector<std::vector<double>>& vecOfVecs) {
    std::vector<double> flattened;
    for (const auto& innerVec : vecOfVecs) {
      flattened.insert(flattened.end(), innerVec.begin(), innerVec.end());
    }
    return flattened;
  }
  
  // Function to generate random integers
  std::vector<int> generateRandomIntegers(int n, int lowerBound, int upperBound) {
    std::vector<int> randomIntegers;
    std::random_device rd;  // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(lowerBound, upperBound); // define the range

    for (size_t i = 0; i < n; ++i) {
      randomIntegers.push_back(distr(gen));
    }

    return randomIntegers;
  }

} // namespace utils
