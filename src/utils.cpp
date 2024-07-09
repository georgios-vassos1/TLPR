#ifdef UTILS_HPP 
  #include "utils.hpp"
#else
  #include "../inst/include/utils.hpp"
#endif

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

  // Helper function to create a vector with elements from 0 to n-1
  std::vector<int> createIndexVector(int n) {
    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i) {
      idx[i] = i;
    }
    return idx;
  }

  // Helper function to append copies of an index vector to a target vector
  void appendIndexVectors(std::vector<std::vector<int>>& target, const std::vector<int>& idx, int count) {
    for (int i = 0; i < count; ++i) {
      target.push_back(idx);
    }
  }

} // namespace utils
