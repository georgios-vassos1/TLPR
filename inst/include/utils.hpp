#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <random>
#include <vector>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <functional>
#include <type_traits>

// Namespace for utility functions
namespace utils {

  // Helper function to print tuples
  template<class... Args>
  void printTuple(const std::tuple<Args...>& t) {
      std::apply([](auto&&... args){((std::cout << args << ", "), ...);}, t);
      std::cout << std::endl;
  }

  // Function to concatenate n vectors
  template <typename T>
  std::vector<T> concatenateVectors(const std::vector<std::vector<T>>& vectors) {
    // Determine the total size needed
    size_t total_size = 0;
    for (const auto& v : vectors) {
      total_size += v.size();
    }

    // Create the concatenated vector with reserved capacity
    std::vector<T> concatenated_vector;
    concatenated_vector.reserve(total_size);

    // Insert elements from each vector
    for (const auto& v : vectors) {
      concatenated_vector.insert(concatenated_vector.end(), v.begin(), v.end());
    }

    return concatenated_vector;
  }

  template <typename T>
  std::vector<T> mirrorAndNegateVector(const std::vector<T>& original) {
      static_assert(std::is_arithmetic<T>::value, "T must be a numeric type");

      std::vector<T> extended;
      extended.reserve(2 * original.size() - 1);

      // Prepend negatives in reverse order
      for (auto it = original.rbegin(); it != original.rend(); ++it) {
          extended.push_back(-*it);
      }

      // Append the original vector
      for (const auto& value : original) {
          extended.push_back(value);
      }

      return extended;
  }

  // Helper function to create a vector with elements from 0 to n-1
  std::vector<int> createIndexVector(int n);

  // Helper function to append copies of an index vector to a target vector
  void appendIndexVectors(std::vector<std::vector<int>>& target, const std::vector<int>& idx, int count);

  // Function to flatten a vector of vectors of equal lengths
  std::vector<double> flatten(const std::vector<std::vector<double>>& vecOfVecs);

  // Function to generate random integers
  std::vector<int> generateRandomIntegers(int n, int lowerBound, int upperBound);

} // namespace utils

#endif // UTILS_HPP
