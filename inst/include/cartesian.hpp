#ifndef CARTESIAN_HPP
#define CARTESIAN_HPP

#include <Eigen/Dense>
#include <vector>

// Function using only STL containers
std::vector<std::vector<int>> CartesianProductIntSTL(const std::vector<std::vector<int>>& vectors);

// Function returning an Eigen matrix to avoid type conversion to Rcpp supported container
Eigen::MatrixXi CartesianProductInt(const std::vector<std::vector<int>>& vectors);

// Function to generate a range of combinations in parallel
void generateCombinations(const std::vector<std::vector<int>>& vectors, Eigen::MatrixXi& results, size_t start, size_t end);

// Function to generate all combinations in parallel
Eigen::MatrixXi CartesianProductIntParallel(const std::vector<std::vector<int>>& vectors, const size_t numThreads);

// Function to generate all combinations in parallel with load balancing
Eigen::MatrixXi CartesianProductIntParallelxLB(const std::vector<std::vector<int>>& vectors, const size_t numThreads);

#endif // CARTESIAN_HPP

