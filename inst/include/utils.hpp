#ifndef UTILS_HPP
#define UTILS_HPP

#include <Rcpp.h>
#include <vector>
#include <tuple>
#include <map>
#include <algorithm>
#include <numeric>
#include <functional>
#include <random>
#include <iostream>
#include <type_traits>

// Namespace for utility functions
namespace utils {

// Helper function to create a vector with elements from 0 to n-1
std::vector<int> createIndexVector(int n);

// Helper function to append copies of an index vector to a target vector
void appendIndexVectors(std::vector<std::vector<int>>& target, const std::vector<int>& idx,
                        int count);

// Function to flatten a vector of vectors of equal lengths
std::vector<double> flatten(const std::vector<std::vector<double>>& vecOfVecs);

// Function to generate random integers
std::vector<int> generateRandomIntegers(int n, int lowerBound, int upperBound);

// Helper function to print tuples
template <class... Args> void printTuple(const std::tuple<Args...>& t) {
  std::apply([](auto&&... args) { ((std::cout << args << ", "), ...); }, t);
  std::cout << '\n';
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

template <typename T> std::vector<T> mirrorAndNegateVector(const std::vector<T>& original) {
  static_assert(std::is_arithmetic_v<T>, "T must be a numeric type");

  std::vector<T> extended;
  extended.reserve(2 * original.size() - 1);

  // Prepend negatives in reverse order
  for (auto it = original.rbegin(); it != original.rend() - 1; ++it) {
    extended.push_back(-*it);
  }

  // Append the original vector
  for (const auto& value : original) {
    extended.push_back(value);
  }

  return extended;
}

// Decompose a destination extended-state index into surplus (S+) and deficit (S-)
// components.  The extended index idx ∈ {0,...,2R} encodes signed inventory
// val = idx - R via mirrorAndNegateVector.  Only one of (sp, sm) is nonzero.
//
//   sp = max(val,  0)  — surplus at terminal  (val > 0)
//   sm = max(-val, 0)  — deficit at terminal  (val < 0)
//
inline std::pair<int, int> splitSignedIdx(int idx, int R) {
  const int val = idx - R;
  return {std::max(val, 0), std::max(-val, 0)};
}

// Inverse of splitSignedIdx: reconstruct the extended-state index from (sp, sm).
// Caller must guarantee sp >= 0, sm >= 0, and sp * sm == 0.
inline int joinSignedIdx(int sp, int sm, int R) {
  return R + sp - sm;
}

// Helper function to append copies of a vector to a target vector
template <typename T>
void appendVectors(std::vector<std::vector<T>>& target, const std::vector<T>& idx, int count) {
  for (int i = 0; i < count; ++i) {
    target.push_back(idx);
  }
}

// Function to convert R list of vectors to std::map<KeyType, std::vector<ValueType>>
template <typename KeyType, typename ValueType>
std::map<KeyType, std::vector<ValueType>> convertListToMap(SEXP inputList) {
  std::map<KeyType, std::vector<ValueType>> result;

  // Ensure the input is a list
  if (!Rf_isNewList(inputList)) {
    Rcpp::stop("Input must be a list");
  }

  Rcpp::List list(inputList);
  const bool useIndices = (list.names() == R_NilValue);
  Rcpp::CharacterVector names{};
  if (!useIndices) {
    names = list.names();
  }

  for (int i = 0; i < list.size(); i++) {
    KeyType key;
    if constexpr (std::is_same_v<KeyType, int>) {
      key = i;
    } else {
      if (useIndices) {
        Rcpp::stop("Named list required for non-integer key type");
      }
      key = Rcpp::as<KeyType>(names[i]);
    }

    const SEXP valuesSEXP = list[i];
    Rcpp::NumericVector values(valuesSEXP);
    std::vector<ValueType> vec(values.begin(), values.end());

    result[key] = vec;
  }

  return result;
}

// Function to print a map with generic key and value types
template <typename KeyType, typename ValueType>
void printMap(const std::map<KeyType, std::vector<ValueType>>& map) {
  for (const auto& pair : map) {
    std::cout << pair.first << ": ";
    for (const auto& val : pair.second) {
      std::cout << val << " ";
    }
    std::cout << '\n';
  }
}

} // namespace utils

#endif // UTILS_HPP
