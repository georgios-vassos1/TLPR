#include <iostream>
#include <future>
#include <thread>
#include <Rcpp.h>

#include "cartesian.hpp"

using namespace std;
using namespace Rcpp;

// Function using only STL containers
std::vector<std::vector<int>> CartesianProductIntSTL(const std::vector<std::vector<int>>& vectors) {
  std::vector<std::vector<int>> results;

  // Compute the total number of combinations and reserve memory
  size_t totalCombinations = 1;
  for (const auto& vector : vectors) {
    if (vector.empty()) return results;
    totalCombinations *= vector.size();
  }
  results.reserve(totalCombinations);

  // Main loop to generate all combinations
  for (size_t count = 0; count < totalCombinations; ++count) {
    std::vector<int> current;
    size_t temp = count;
    for (size_t i = 0; i < vectors.size(); ++i) {
      current.push_back(vectors[i][temp % vectors[i].size()]);
      temp /= vectors[i].size();
    }
    results.push_back(current);
  }

  return results;
}

// Function returning an Eigen matrix to avoid type conversion to Rcpp supported container
Eigen::MatrixXi CartesianProductInt(const std::vector<std::vector<int>>& vectors) {
  if (vectors.empty()) return Eigen::MatrixXi(0, 0);

  // Compute the total number of combinations
  size_t totalCombinations = 1;
  for (const auto& vec : vectors) {
    if (vec.empty()) return Eigen::MatrixXi(0, vectors.size()); // Early return if any input vector is empty
    totalCombinations *= vec.size();
  }

  // Create an Eigen matrix to hold the results
  Eigen::MatrixXi results(totalCombinations, vectors.size());

  // Generate all combinations
  for (size_t count = 0; count < totalCombinations; ++count) {
    size_t temp = count;
    for (size_t i = 0; i < vectors.size(); ++i) {
      results(count, i) = vectors[i][temp % vectors[i].size()];
      temp /= vectors[i].size();
    }
  }

  return results;
}

// Function to generate a range of combinations in parallel
void generateCombinations(const std::vector<std::vector<int>>& vectors, Eigen::MatrixXi& results, size_t start, size_t end) {

  for (size_t count = start; count < end; ++count) {
    size_t temp = count;

    for (size_t i = 0; i < vectors.size(); ++i) {
      results(count, i) = vectors[i][temp % vectors[i].size()];
      temp /= vectors[i].size();
    }
  }
}

// Function to generate all combinations in parallel
Eigen::MatrixXi CartesianProductIntParallel(const std::vector<std::vector<int>>& vectors, const size_t numThreads) {
  if (vectors.empty()) return Eigen::MatrixXi(0, 0);

  // Compute the total number of combinations
  size_t totalCombinations = 1;
  for (const auto& vec : vectors) {
    if (vec.empty()) return Eigen::MatrixXi(0, vectors.size()); // Early return if any input vector is empty
    totalCombinations *= vec.size();
  }

  // Create an Eigen matrix to hold the results
  Eigen::MatrixXi results(totalCombinations, vectors.size());

  // Determine the range of combinations each thread will handle
  size_t combinationsPerThread = totalCombinations / numThreads;
  size_t remainder = totalCombinations % numThreads;

  // Launch threads
  std::vector<std::thread> threads;
  size_t start = 0;
  for (size_t i = 0; i < numThreads; ++i) {
    size_t end = start + combinationsPerThread + (i < remainder ? 1 : 0);
    threads.emplace_back(generateCombinations, std::cref(vectors), std::ref(results), start, end);
    start = end;
  }

  // Join threads
  for (auto& thread : threads) {
    thread.join();
  }

  return results;
}

// Function to generate all combinations in parallel using load balancing between threads
Eigen::MatrixXi CartesianProductIntParallelxLB(const std::vector<std::vector<int>>& vectors, const size_t numThreads) {
  if (vectors.empty()) return Eigen::MatrixXi(0, 0);

  // Compute the total number of combinations
  size_t totalCombinations = 1;
  for (const auto& vec : vectors) {
    if (vec.empty()) return Eigen::MatrixXi(0, vectors.size()); // Early return if any input vector is empty
    totalCombinations *= vec.size();
  }

  // Create an Eigen matrix to hold the results
  Eigen::MatrixXi results(totalCombinations, vectors.size());
  std::vector<std::future<void>> futures;

  // Determine the range of combinations each thread will handle
  size_t combinationsPerThread = totalCombinations / numThreads;
  size_t remainder = totalCombinations % numThreads;

  // Launch threads
  size_t start = 0;
  for (size_t i = 0; i < numThreads; ++i) {
    size_t end = start + combinationsPerThread + (i < remainder ? 1 : 0);
    futures.push_back(
      std::async(std::launch::async, generateCombinations, std::cref(vectors), std::ref(results), start, end)
    );
    start = end;
  }

  // Join threads
  for (auto& fut : futures) {
    fut.get();
  }

  return results;
}

// Wrapper function to call CartesianProductInt from R
// [[Rcpp::export]]
Eigen::MatrixXi CartesianProductRcpp(List vectors) {
    std::vector<std::vector<int>> cpp_vectors;
    for (int i = 0; i < vectors.size(); ++i) {
        cpp_vectors.push_back(as<std::vector<int>>(vectors[i]));
    }

    return CartesianProductInt(cpp_vectors);
}

// Wrapper function to call CartesianProductIntParallel from R
// [[Rcpp::export]]
Eigen::MatrixXi CartesianProductRcppParallel(List vectors, int numThreads) {
    std::vector<std::vector<int>> cpp_vectors;
    for (int i = 0; i < vectors.size(); ++i) {
        cpp_vectors.push_back(as<std::vector<int>>(vectors[i]));
    }

    return CartesianProductIntParallel(cpp_vectors, numThreads);
}

// Wrapper function to call CartesianProductIntParallel from R
// [[Rcpp::export]]
Eigen::MatrixXi CartesianProductRcppParallelxLB(List vectors, int numThreads) {
    std::vector<std::vector<int>> cpp_vectors;
    for (int i = 0; i < vectors.size(); ++i) {
        cpp_vectors.push_back(as<std::vector<int>>(vectors[i]));
    }

    return CartesianProductIntParallelxLB(cpp_vectors, numThreads);
}

// int main(int argc, char** argv) {
//   // Example usage
//   std::vector<std::vector<int>> vectors = {{1, 2, 3}, {4, 5}, {6, 7, 8}};
// 
//   // Compute the Cartesian product using STL containers
//   std::vector<std::vector<int>> resultsSTL = CartesianProductIntSTL(vectors);
//   // Print the results
//   std::cout << "STL results:" << std::endl;
//   for (const auto& v : resultsSTL) {
//     for (const auto& elem : v) {
//       std::cout << elem << " ";
//     }
//     std::cout << std::endl;
//   }
//   std::cout << std::endl;
// 
//   // Compute the Cartesian product using Eigen
//   Eigen::MatrixXi resultsEigen = CartesianProductInt(vectors);
//   // Print the results
//   std::cout << "Eigen results:" << std::endl;
//   std::cout << resultsEigen << std::endl;
//   std::cout << std::endl;
// 
//   // Compute the Cartesian product using Eigen in parallel
//   Eigen::MatrixXi resultsEigenParallel = CartesianProductIntParallel(vectors, 4);
//   // Print the results
//   std::cout << "Eigen parallel results:" << std::endl;
//   std::cout << resultsEigenParallel << std::endl;
//   std::cout << std::endl;
// 
//   return 0;
// }
