#include "gurobi_c++.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using namespace std;

// Function to flatten a vector of vectors of equal lengths
std::vector<double> flatten(const std::vector<std::vector<double>>& vecOfVecs) {
  std::vector<double> flattened;
  for (const auto& innerVec : vecOfVecs) {
    flattened.insert(flattened.end(), innerVec.begin(), innerVec.end());
  }
  return flattened;
}

// Function to extract R list from the input JSON
template<typename T>
std::unordered_map<std::string, std::vector<T>> importList(const nlohmann::json& input) {
  std::unordered_map<std::string, std::vector<T>> output;

  // Extract and iterate over the "winner" object
  for (auto it = input.begin(); it != input.end(); ++it) {
    // Extract the key (as string) and value (as array) pair
    std::string key = it.key();
    std::vector<T> values = it.value().get<std::vector<T>>();

    // Store the key-value pair into the unordered_map
    output[key] = values;
  }

  return output;
}


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
    for (int bidIndex : winner.second) {
      std::cout << "won bid index " << bidIndex << "!\n";
      const auto& bid = bids[bidIndex - 1];

      for (int laneIndex : bid) {
        std::cout << "Lane index " << laneIndex << ": " << std::endl;
        const auto& lane = lanes[laneIndex - 1];
        std::cout << "From origin " << lane[0] << " to destination " << lane[1] << " at " << contractRates.at(winnerKey)[k++] << " USD per km." << std::endl;
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
      std::cout << "From origin " << lane[0] << " to destination " << lane[1] << " at " << spotRates[carrierIndex * nLanes + laneIndex] << " USD per km." << std::endl;
    }
  }
}

// Function to create decision variables for transportation
GRBVar* createTransportVars(
    GRBModel& model, 
    const int& n, 
    const std::unordered_map<std::string, std::vector<int>>& winners,
    const std::vector<std::vector<int>>& bids,
    const std::vector<std::vector<int>>& lanes,
    const std::unordered_map<std::string, std::vector<double>>& contractRates,
    const std::vector<double>& spotRates,
    const int& nSpotCarriers,
    const int& nLanes) {

  GRBVar* transport = new GRBVar[n];

  int idx = 0;
  // Strategic carrier variables
  for (const auto& winner : winners) {
    std::string winnerKey = winner.first;

    int k = 0;
    for (int bidIndex : winner.second) {
      const auto& bid = bids[bidIndex - 1];

      for (int laneIndex : bid) {
        const auto& lane = lanes[laneIndex - 1];

        ostringstream vname;
        vname << "StrategicCarrier" << winnerKey << ".Bid" << bidIndex << ".Lane" << laneIndex << ".from_" << lane[0] << "_to_" << lane[1];

        transport[idx++] = model.addVar(0.0, GRB_INFINITY, contractRates.at(winnerKey)[k++], GRB_CONTINUOUS, vname.str());
      }
    }
  }

  // Spot carrier variables
  for (int carrierIndex = 0; carrierIndex < nSpotCarriers; carrierIndex++) {

    for (int laneIndex = 0; laneIndex < nLanes; laneIndex++) {
      const auto& lane = lanes[laneIndex];

      ostringstream vname;
      vname << "SpotCarrier" << carrierIndex + 1 << ".Lane" << laneIndex + 1 << ".from_" << lane[0] << "_to_" << lane[1];

      transport[idx++] = model.addVar(0.0, GRB_INFINITY, spotRates[carrierIndex * nLanes + laneIndex], GRB_CONTINUOUS, vname.str());
    }
  }

  return transport;
}

// Volume constraint
void addVolumeConstraint(
    GRBModel& model, 
    GRBVar* transport, 
    const int& n, 
    const double& At) {

  GRBLinExpr expr;
  for (int i = 0; i < n; i++) {
      expr += transport[i];
  }

  model.addConstr(expr == At, "VolumeConstraint");
}

void addCapacityConstraints(
    GRBModel& model, 
    GRBVar* transport,
    const int& stage,
    const std::unordered_map<std::string, std::vector<int>>& winners,
    const std::vector<std::vector<int>>& bids,
    const int& nSpotCarriers,
    const int& nLanes,
    const std::unordered_map<std::string, std::vector<int>>& capacities) {

  int k = 0;

  // Strategic carrier capacity constraints
  for (auto& winner : winners) {
    std::string winnerKey = winner.first;
    const int capacity = capacities.at(winnerKey).at(stage);

    GRBLinExpr expr = 0;
    for (int bidIndex : winner.second) {
      const auto& bid = bids[bidIndex - 1];

      for (int i = 0; i < bid.size(); i++) {
        expr += transport[k++];
      }
    }
    model.addConstr(expr <= capacity, "CapacityConstraint");
  }

  // Spot carrier capacity constraints
  for (int carrierIndex = 0; carrierIndex < nSpotCarriers; carrierIndex++) {
    const int capacity = 20;

    GRBLinExpr expr = 0;
    for (int laneIndex = 0; laneIndex < nLanes; laneIndex++) {
      expr += transport[k++];
    }
    model.addConstr(expr <= capacity, "CapacityConstraint");
  }
}


int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return(EXIT_FAILURE);
  }

  std::ifstream file(argv[1]);
  nlohmann::json input;
  file >> input;

  // Extract the data
  // const int tau = input["tau"][0];
  const int nL_ = input["nL_"][0];
  const int nCO = input["nCO"][0];
  const int nL  = input["nL"][0];
  const std::vector<std::vector<int>> bids = input["B"];
  const std::vector<std::vector<int>> lanes = input["L"];

  std::vector<std::vector<double>> spotRates = input["CTo"];

  std::unordered_map<std::string, std::vector<int>> winners;
  std::unordered_map<std::string, std::vector<double>> CTb;
  std::unordered_map<std::string, std::vector<int>> Cb;

  const int n = nL_ + nCO * nL;
  // Extract and iterate over the "winner" object
  // for (auto it = input["winner"].begin(); it != input["winner"].end(); ++it) {
  //   // Extract the key (as string) and value (as array) pair
  //   std::string key = it.key();
  //   std::vector<int> values = it.value().get<std::vector<int>>();

  //   // Store the key-value pair into the unordered_map
  //   winners[key] = values;
  // }
  winners = importList<int>(input["winner"]);
  CTb = importList<double>(input["CTb_list"]);
  Cb  = importList<int>(input["Cb_list"]);

  // Print the strategy of the instance
  // printInstance(winners, bids, lanes, CTb, spotRates[0], nCO, nL);

  // Create the Gurobi environment
  GRBEnv env;
  GRBVar* transport = nullptr;

  try {
    GRBModel model = GRBModel(env);
    model.set(GRB_StringAttr_ModelName, "Drayage");

    transport = createTransportVars(model, n, winners, bids, lanes, CTb, spotRates[0], nCO, nL);

    addVolumeConstraint(model, transport, n, 40);
    addCapacityConstraints(model, transport, 0, winners, bids, nCO, nL, Cb);

    // Use barrier to solve root relaxation
    model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Solve
    model.optimize();
    int status = model.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL) {
      std::cout << "Optimal solution found!" << std::endl;
      // Print transport names and values
      std::cout << "Transport volumes:" << std::endl;
      for (int i = 0; i < n; ++i) {
        try {
          std::cout << transport[i].get(GRB_StringAttr_VarName) << " = " << transport[i].get(GRB_DoubleAttr_X) << std::endl;
        } catch (GRBException& e) {
          std::cerr << "Error code = " << e.getErrorCode() << std::endl;
          std::cerr << e.getMessage() << std::endl;
        }
      }

    } else if (status == GRB_INFEASIBLE) {
      std::cout << "Model is infeasible" << std::endl;
    } else if (status == GRB_UNBOUNDED) {
      std::cout << "Model is unbounded" << std::endl;
    } else {
      std::cout << "Optimization ended with status " << status << std::endl;
    }

  } catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }


  delete[] transport;

  return 0;
}
