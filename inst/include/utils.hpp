#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <random>
#include <vector>
#include <tuple>

// Namespace for utility functions
namespace utils {

// Helper function to print tuples
template<class... Args>
void printTuple(const std::tuple<Args...>& t) {
    std::apply([](auto&&... args){((std::cout << args << ", "), ...);}, t);
    std::cout << std::endl;
}

// Function to flatten a vector of vectors of equal lengths
std::vector<double> flatten(const std::vector<std::vector<double>>& vecOfVecs);

// Function to generate random integers
std::vector<int> generateRandomIntegers(int n, int lowerBound, int upperBound);

} // namespace utils

#endif // UTILS_HPP
