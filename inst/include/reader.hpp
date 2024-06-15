#ifndef READER_HPP
#define READER_HPP

#include <nlohmann/json.hpp>
#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>

// Function to extract an unordered map (dictionary) from the input JSON
// @param input JSON object to parse
// @return Unordered map with string keys and values of type T
template<typename T>
std::unordered_map<std::string, T> importList(const nlohmann::json& input);

// Function to extract a list of vectors from the input JSON
// @param input JSON object to parse
// @return Unordered map with string keys and vector of type T as values
template<typename T>
std::unordered_map<std::string, std::vector<T>> importListOfVectors(const nlohmann::json& input);

// Function to print the instance
// @param winners Map of winners with vector of integers
// @param bids List of vectors representing bids
// @param lanes List of vectors representing lanes
// @param contractRates Map of contract rates with vector of doubles
// @param spotRates Vector of spot rates
// @param nSpotCarriers Number of spot carriers
// @param nLanes Number of lanes
void printInstance(const std::unordered_map<std::string, std::vector<int>>& winners,
                   const std::vector<std::vector<int>>& bids,
                   const std::vector<std::vector<int>>& lanes,
                   const std::unordered_map<std::string, std::vector<double>>& contractRates,
                   const std::vector<double>& spotRates,
                   const int& nSpotCarriers,
                   const int& nLanes);

// Function to extract R list from the input JSON
template<typename T>
std::unordered_map<std::string, T> importList(const nlohmann::json& input) {
  std::unordered_map<std::string, T> output;

  // Extract and iterate over the "winner" object
  for (auto it = input.begin(); it != input.end(); ++it) {
    // Extract the key (as string) and value (as array) pair
    std::string key = it.key();
    // Ensure the value is an array with exactly one integer element
    if (!it.value().is_array() || it.value().size() != 1 || !it.value()[0].is_number_integer()) {
      throw std::invalid_argument("Expected an array with a single integer for key: " + key);
    }
    // Attempt to extract the value as type T
    try {
      T value = it.value()[0].get<T>(); // Explicit conversion using get<T>()
      output[key] = value;
    } catch (nlohmann::json::type_error& e) {
      throw std::runtime_error("Type mismatch when converting JSON value to type T: " + std::string(e.what()));
    }
  }

  return output;
}

// Function to extract a list of vectors from the input JSON
template<typename T>
std::unordered_map<std::string, std::vector<T>> importListOfVectors(const nlohmann::json& input) {
  std::unordered_map<std::string, std::vector<T>> output;

  // Check that input is an object
  if (!input.is_object()) {
    throw std::invalid_argument("JSON input is not an object");
  }

  // Iterate over each key-value pair in the JSON object
  for (auto it = input.begin(); it != input.end(); ++it) {
    std::string key = it.key();

    // Check if the value is an array
    if (!it.value().is_array()) {
      throw std::invalid_argument("Expected an array for key: " + key);
    }

    // Attempt to extract the value as a vector of type T
    try {
      std::vector<T> values = it.value().get<std::vector<T>>();
      output[key] = values;
    } catch (nlohmann::json::type_error& e) {
      throw std::runtime_error("Type mismatch when converting JSON array elements to type T for key '" + key + "': " + std::string(e.what()));
    }
  }

  return output;
}

#endif // READER_HPP
