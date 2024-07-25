#include <iostream>
#include <fstream>

#ifdef READER_HPP 
  #include "reader.hpp"
#else
  #include "../inst/include/reader.hpp"
#endif

using namespace std;

// Function to print the instance
void printInstance(const std::unordered_map<std::string, std::vector<int>>& winners,
                   const std::vector<std::vector<int>>& bids,
                   const std::vector<std::vector<int>>& lanes,
                   const std::unordered_map<std::string, std::vector<double>>& contractRates,
                   const std::vector<double>& spotRates,
                   const int& nSpotCarriers,
                   const int& nLanes) {

  // Print winner key-value pairs
  for (const auto& winner : winners) {
    std::string winnerKey = winner.first;
    std::cout << "Carrier index " << winnerKey << " ";

    int k = 0;
    for (size_t bidIndex : winner.second) {
      std::cout << "won bid index " << bidIndex << "!\n";
      const auto& bid = bids[bidIndex - 1];

      for (size_t laneIndex : bid) {
        std::cout << "Lane index " << laneIndex << ": " << std::endl;
        const auto& lane = lanes[laneIndex - 1];
        std::cout << "From origin "     << lane[0] 
                  << " to destination " << lane[1] 
                  << " at " << contractRates.at(winnerKey).at(k++) << " USD per km." << std::endl;
      }
    }

    std::cout << std::endl;
  }

  // Print spot carrier data
  for (int carrierIndex = 0; carrierIndex < nSpotCarriers; carrierIndex++) {
    std::cout << "Spot carrier index " << carrierIndex + 1 << " ";

    for (int laneIndex = 0; laneIndex < nLanes; laneIndex++) {

      std::cout << "Lane index " << laneIndex + 1 << ": " << std::endl;
      const auto& lane = lanes[laneIndex];
      std::cout << "From origin "     << lane[0] 
                << " to destination " << lane[1] 
                << " at " << spotRates[carrierIndex * nLanes + laneIndex] << " USD per km." << std::endl;
    }
  }
}


// int main(int argc, char** argv) {
//   if (argc < 2) {
//     std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
//     return(EXIT_FAILURE);
//   }
// 
//   std::ifstream file(argv[1]);
//   nlohmann::json input;
//   file >> input;
// 
//   // Extract the data
//   const int nCO = input["nCO"][0];
//   const int nL  = input["nL"][0];
//   const std::vector<std::vector<int>> bids = input["B"];
//   const std::vector<std::vector<int>> lanes = input["L"];
//   const std::vector<std::vector<int>> Cb = input["Cb"];
//   const std::vector<std::vector<int>> Co = input["Co"];
// 
//   const std::vector<std::vector<double>> spotRates = input["CTo"];
// 
//   std::unordered_map<std::string, std::vector<int>> winners;
//   std::unordered_map<std::string, std::vector<double>> CTb;
//   // std::unordered_map<std::string, std::vector<int>> Cb;
//   std::unordered_map<std::string, int> carrierIdx;
// 
//   winners = importListOfVectors<int>(input["winner"]);
//   CTb = importListOfVectors<double>(input["CTb_list"]);
//   carrierIdx = importList<int>(input["carrierIdx"]);
// 
//   // Print the strategy of the instance
//   printInstance(winners, bids, lanes, CTb, spotRates[0], nCO, nL);
// 
//   return 0;
// }
