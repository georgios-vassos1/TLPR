#include "gurobi_c++.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using namespace std;
using json  = nlohmann::json;

// Function to print the instance
void printInstance(const std::unordered_map<std::string, std::vector<int>>& winners,
                   const std::vector<std::vector<int>>& bids,
                   const std::vector<std::vector<int>>& lanes) {
  // Print winner key-value pairs
  for (const auto& winner : winners) {
    std::cout << "Carrier index " << winner.first << " ";

    for (int bidIndex : winner.second) {
      std::cout << "won bid index " << bidIndex << "!\n";
      const auto& bid = bids[bidIndex - 1];

      for (int laneIndex : bid) {
        std::cout << "Lane index " << laneIndex << ": " << std::endl;
        const auto& lane = lanes[laneIndex - 1];
        std::cout << "From origin " << lane[0] << " to destination " << lane[1] << std::endl;
      }
    }

    std::cout << std::endl;
  }
}

// Function to create decision variables for transportation
GRBVar* createTransportVars(
    GRBModel& model, 
    const int& n, 
    const std::unordered_map<std::string, std::vector<int>>& winners,
    const std::vector<std::vector<int>>& bids,
    const std::vector<std::vector<int>>& lanes,
    const int& nSpotCarriers,
    const int& nLanes) {

  GRBVar* transport = new GRBVar[n];

  int idx = 0;
  // Strategic carrier variables
  for (const auto& winner : winners) {

    for (int bidIndex : winner.second) {
      const auto& bid = bids[bidIndex - 1];

      for (int laneIndex : bid) {
        const auto& lane = lanes[laneIndex - 1];

        ostringstream vname;
        vname << "Carrier" << winner.first << ".Bid" << bidIndex << ".Lane" << laneIndex << ".from_" << lane[0] << "_to_" << lane[1];

        transport[idx++] = model.addVar(0.0, GRB_INFINITY, 0.1, GRB_CONTINUOUS, vname.str());
      }
    }
  }

  // Spot carrier variables
  for (int carrierIndex = 0; carrierIndex < nSpotCarriers; carrierIndex++) {

    for (int laneIndex = 0; laneIndex < nLanes; laneIndex++) {
      const auto& lane = lanes[laneIndex];

      ostringstream vname;
      vname << "Carrier" << carrierIndex << ".Lane" << laneIndex << ".from_" << lane[0] << "_to_" << lane[1];

      transport[idx++] = model.addVar(0.0, GRB_INFINITY, (double)rand() / RAND_MAX, GRB_CONTINUOUS, vname.str());
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


int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return(EXIT_FAILURE);
  }

  std::ifstream file(argv[1]);
  json input;
  file >> input;

  // Extract the data
  const int nL_ = input["nL_"][0];
  const int nCO = input["nCO"][0];
  const int nL  = input["nL"][0];
  const std::vector<double> costs = input["CTb"];
  const std::vector<std::vector<int>> bids = input["B"];
  const std::vector<std::vector<int>> lanes = input["L"];
  std::unordered_map<std::string, std::vector<int>> winners;

  const int n = nL_ + nCO * nL;
  // Extract and iterate over the "winner" object
  for (auto it = input["winner"].begin(); it != input["winner"].end(); ++it) {
      // Extract the key (as string) and value (as array) pair
      std::string key = it.key();
      std::vector<int> values = it.value().get<std::vector<int>>();

      // Store the key-value pair into the unordered_map
      winners[key] = values;
  }

  // Print the strategy of the instance
  // printInstance(winners, bids, lanes);

  // Create the Gurobi environment
  GRBEnv env;
  GRBVar* transport = nullptr;

  try {
    GRBModel model = GRBModel(env);
    model.set(GRB_StringAttr_ModelName, "Drayage");

    transport = createTransportVars(model, n, winners, bids, lanes, nCO, nL);

    addVolumeConstraint(model, transport, n, 100);

    // Use barrier to solve root relaxation
    model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);

    // Solve
    model.optimize();

    // Print transport names and values
    for (int i = 0; i < n; ++i) {
      try {
        std::cout << transport[i].get(GRB_StringAttr_VarName) << " = " << transport[i].get(GRB_DoubleAttr_X) << std::endl;
      } catch (GRBException& e) {
        std::cerr << "Error code = " << e.getErrorCode() << std::endl;
        std::cerr << e.getMessage() << std::endl;
      }
    }

    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

  } catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }


  delete[] transport;

  return 0;
}
