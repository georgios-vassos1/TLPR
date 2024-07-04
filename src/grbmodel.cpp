#include "gurobi_c++.h"
#include <fstream>
#include <sstream>
#include <thread>
#include <Eigen/Dense>
#include <Rcpp.h>

#include "reader.hpp"
#include "utils.hpp"

using namespace std;


// Function to print the constraints in a gurobi model
void printConstraints(GRBModel& model) {

  GRBConstr* constraints = model.getConstrs();

  int numConstraints = model.get(GRB_IntAttr_NumConstrs);
  for (size_t i = 0; i < numConstraints; ++i) {
    std::cout << constraints[i].get(GRB_StringAttr_ConstrName) << " = " << constraints[i].get(GRB_DoubleAttr_RHS) << std::endl;
  }

}

// Function to print the optimal transport volumes from gurobi model optimization
void printOptimalTransportVolumes(GRBVar* transport, const int& n) {

  std::cout << "Transport volumes:" << std::endl;
  for (size_t i = 0; i < n; ++i) {
    try {
      std::cout << transport[i].get(GRB_StringAttr_VarName) << " = " << transport[i].get(GRB_DoubleAttr_X) << std::endl;
    } catch (GRBException& e) {
      std::cerr << "Error code = " << e.getErrorCode() << std::endl;
      std::cerr << e.getMessage() << std::endl;
    }
  }

}

// Function to create decision variables for transportation
GRBVar* createTransportVars(
    GRBModel& model, 
    const int& n, 
    const std::vector<std::string>& winnerKeys, // This is only used to maintain the original order of the carriers
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
  for (const auto& winnerKey : winnerKeys) {
    // std::string winnerKey = winner.first;

    int k = 0;
    // for (size_t bidIndex : winner.second) {
    for (size_t bidIndex : winners.at(winnerKey)) {
      const auto& bid = bids[bidIndex - 1];

      for (size_t laneIndex : bid) {
        const auto& lane = lanes[laneIndex - 1];

        ostringstream vname;
        vname << "StrategicCarrier" << winnerKey 
              << ".Bid"   << bidIndex 
              << ".Lane"  << laneIndex 
              << ".from_" << lane[0] 
              << "_to_"   << lane[1];

        transport[idx++] = model.addVar(0.0, GRB_INFINITY, contractRates.at(winnerKey).at(k++), GRB_CONTINUOUS, vname.str());
      }
    }
  }

  // Spot carrier variables
  for (size_t carrierIndex = 0; carrierIndex < nSpotCarriers; carrierIndex++) {

    for (size_t laneIndex = 0; laneIndex < nLanes; laneIndex++) {
      const auto& lane = lanes[laneIndex];

      ostringstream vname;
      vname << "SpotCarrier" << carrierIndex + 1 
            << ".Lane" << laneIndex + 1 
            << ".from_" << lane[0] 
            << "_to_" << lane[1];

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
  for (size_t i = 0; i < n; i++) {
      expr += transport[i];
  }

  model.addConstr(expr == At, "VolumeConstraint");
}

void addCapacityConstraints(
    GRBModel& model, 
    GRBVar* transport,
    const std::vector<std::string>& winnerKeys,
    const std::unordered_map<std::string, std::vector<int>>& winners,
    const std::vector<std::vector<int>>& bids,
    const std::unordered_map<std::string, int>& carrierIdx,
    const int& nSpotCarriers,
    const int& nLanes,
    const std::vector<int>& strategicCaps,
    const std::vector<int>& spotCaps) {

  int k = 0;

  // Strategic carrier capacity constraints
  for (auto& winnerKey : winnerKeys) {
    // std::string winnerKey = winner.first;
    const int capacity = strategicCaps.at(carrierIdx.at(winnerKey) - 1);

    GRBLinExpr expr = 0;
    for (size_t bidIndex : winners.at(winnerKey)) {
      const auto& bid = bids[bidIndex - 1];

      for (size_t i = 0; i < bid.size(); ++i) {
        expr += transport[k++];
      }
    }

    ostringstream cname;
    cname << "CapacityConstraint_StrategicCarrier" << winnerKey;

    model.addConstr(expr <= capacity, cname.str());
  }

  // Spot carrier capacity constraints
  for (size_t spotCarrierIdx = 0; spotCarrierIdx < nSpotCarriers; ++spotCarrierIdx) {
    const int capacity = spotCaps.at(spotCarrierIdx);

    GRBLinExpr expr = 0;
    for (size_t laneIndex = 0; laneIndex < nLanes; ++laneIndex) {
      expr += transport[k++];
    }

    ostringstream cname;
    cname << "CapacityConstraint_SpotCarrier" << spotCarrierIdx + 1;

    model.addConstr(expr <= capacity, cname.str());
  }
}

void addStorageLimitConstraints(
    GRBModel& model, 
    GRBVar* transport,
    const std::vector<std::string>& winnerKeys,
    const std::unordered_map<std::string, std::vector<int>>& winners,
    const std::vector<std::vector<int>>& bids,
    const std::vector<std::vector<int>>& lanes,
    const int& nSpotCarriers,
    const int& nLanes,
    const std::vector<int>& nWarehouses, // (nOrigins, nDestinations)
    const std::vector<int>& limits) {

  int position = 0;
  for (size_t m = 0; m < nWarehouses.size(); ++m) {

    for (size_t i = 0; i < nWarehouses[m]; ++i) {
      GRBLinExpr expr = 0;
      int k = 0;

      // Strategic carrier capacity constraints
      for (const auto& winnerKey : winnerKeys) {

        for (size_t bidIndex : winners.at(winnerKey)) {
          const auto& bid = bids[bidIndex - 1];

          for (size_t laneIndex : bid) {
            const auto& lane = lanes[laneIndex - 1];

            if (lane[m] == i + 1) {
              expr += transport[k];
            }
            k++;
          }
        }
      }

      // Spot carrier capacity constraints
      for (size_t carrierIndex = 0; carrierIndex < nSpotCarriers; ++carrierIndex) {

        for (size_t laneIndex = 0; laneIndex < nLanes; ++laneIndex) {
          const auto& lane = lanes[laneIndex];

          if (lane[m] == i + 1) {
            expr += transport[k];
          }
          k++;
        }
      }

      ostringstream cname;
      cname << "StorageLimitConstraint" << m << "_" << i;

      model.addConstr(expr <= limits[position + i], cname.str());
    }
    position += nWarehouses[m];
  }
}

// Optimize the model and reurn objective and decision in R list
// [[Rcpp::export]]
Rcpp::List optimizeModelFromJSON(std::string jsonFile) {

  // Read the JSON file
  std::ifstream file(jsonFile);
  nlohmann::json input;
  file >> input;

  // Extract the data
  const int nOrigins = input["nI"][0];
  const int nDestinations = input["nJ"][0];
  const int nL_ = input["nL_"][0];
  const int nCO = input["nCO"][0];
  const int nL  = input["nL"][0];
  const std::vector<std::vector<int>> bids = input["B"];
  const std::vector<std::vector<int>> lanes = input["L"];
  const std::vector<std::vector<int>> Cb = input["Cb"];
  const std::vector<std::vector<int>> Co = input["Co"];
  const std::vector<std::vector<double>> spotRates = input["CTo"];

  std::unordered_map<std::string, std::vector<int>> winners;
  std::unordered_map<std::string, std::vector<double>> CTb;
  std::unordered_map<std::string, int> carrierIdx;

  const std::vector<std::string> winnerKeys = input["winnerKey"];

  winners = importListOfVectors<int>(input["winner"]);
  CTb = importListOfVectors<double>(input["CTb_list"]);
  carrierIdx = importList<int>(input["carrierIdx"]);

  const int n = nL_ + nCO * nL;
  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  const std::vector<int> limits = utils::generateRandomIntegers(nOrigins + nDestinations, 20, 40);
  // Print limits 
  // std::cout << "Storage limits: ";
  // std::copy(limits.begin(), limits.end(), std::ostream_iterator<int>(std::cout, " "));
  // std::cout << std::endl;

  // Create the Gurobi environment
  GRBEnv env(true);
  // Set the OutputFlag to 0 to turn off logging
  env.set(GRB_IntParam_OutputFlag, 0);

  env.start();

  GRBVar* transport = nullptr;

  // Preallocate the result list
  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["objval"] = 0.0,  // Placeholder for objective value
    Rcpp::_["x"] = Rcpp::NumericVector(n),  // Placeholder for decision variables 'x'
    Rcpp::_["limits"] = limits  // Placeholder for decision variables 'x'
  );

  // Optimize the model
  try {

    GRBModel model = GRBModel(env);
    model.set(GRB_StringAttr_ModelName, "Drayage");

    transport = createTransportVars(model, n, winnerKeys, winners, bids, lanes, CTb, spotRates[0], nCO, nL);

    addVolumeConstraint(model, transport, n, 40);
    addCapacityConstraints(model, transport, winnerKeys,  winners, bids, carrierIdx, nCO, nL, Cb[0], Co[0]);
    addStorageLimitConstraints(model, transport, winnerKeys, winners, bids, lanes, nCO, nL, nWarehouses, limits);

    // Use barrier to solve root relaxation
    model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Solve
    model.optimize();

    int status = model.get(GRB_IntAttr_Status);

    if (status == GRB_OPTIMAL) {
      // std::cout << "Optimal solution found!" << std::endl;

      result["objval"] = model.get(GRB_DoubleAttr_ObjVal);

      Rcpp::NumericVector x(n);
      for (size_t i = 0; i < n; ++i) {
        x[i] = transport[i].get(GRB_DoubleAttr_X);
      }
      result["x"] = x;

    } else if (status == GRB_INFEASIBLE) {
      std::cout << "Model is infeasible" << std::endl;
    } else if (status == GRB_UNBOUNDED) {
      std::cout << "Model is unbounded" << std::endl;
    } else {
      std::cout << "Optimization ended with status " << status << std::endl;
    }

    delete [] transport;

    return result;

  } catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }

  return Rcpp::List::create();
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
  const int nOrigins = input["nI"][0];
  const int nDestinations = input["nJ"][0];
  // const int nL_ = input["nL_"][0];
  // const int nCO = input["nCO"][0];
  // const int nL  = input["nL"][0];
  const std::vector<std::vector<int>> bids = input["B"];
  const std::vector<std::vector<int>> lanes = input["L"];
  const std::vector<std::vector<int>> Cb = input["Cb"];
  const std::vector<std::vector<int>> Co = input["Co"];

  const std::vector<std::vector<double>> spotRates = input["CTo"];

  std::unordered_map<std::string, std::vector<int>> winners;
  std::unordered_map<std::string, std::vector<double>> CTb;
  // std::unordered_map<std::string, std::vector<int>> Cb;
  std::unordered_map<std::string, int> carrierIdx;

  // const int n = nL_ + nCO * nL;
  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  // Extract and iterate over the "winner" object
  // for (auto it = input["winner"].begin(); it != input["winner"].end(); ++it) {
  //   // Extract the key (as string) and value (as array) pair
  //   std::string key = it.key();
  //   std::vector<int> values = it.value().get<std::vector<int>>();

  //   // Store the key-value pair into the unordered_map
  //   winners[key] = values;
  // }
  const std::vector<std::string> winnerKeys = input["winnerKey"];

  winners = importListOfVectors<int>(input["winner"]);
  CTb = importListOfVectors<double>(input["CTb_list"]);
  carrierIdx = importList<int>(input["carrierIdx"]);

  // Print carrier index
  // for (const auto& carrier : carrierIdx) {
  //   std::cout << "Carrier " << carrier.first << " has index " << carrier.second << std::endl;
  // }

  const std::vector<int> limits = utils::generateRandomIntegers(nOrigins + nDestinations, 20, 40);

  // Print the strategy of the instance
  // printInstance(winners, bids, lanes, CTb, spotRates[0], nCO, nL);

  // Create the Gurobi environment
  // GRBEnv env;
  // GRBVar* transport = nullptr;

  // try {
  //   GRBModel model = GRBModel(env);
  //   model.set(GRB_StringAttr_ModelName, "Drayage");

  //   transport = createTransportVars(model, n, winnerKeys, winners, bids, lanes, CTb, spotRates[0], nCO, nL);

  //   addVolumeConstraint(model, transport, n, 40);
  //   addCapacityConstraints(model, transport, winnerKeys, winners, bids, carrierIdx, nCO, nL, Cb[0], Co[0]);
  //   addStorageLimitConstraints(model, transport, winnerKeys, winners, bids, lanes, nCO, nL, nWarehouses, limits);

  //   // Use barrier to solve root relaxation
  //   model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
  //   model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

  //   // Solve
  //   model.optimize();

  //   int status = model.get(GRB_IntAttr_Status);

  //   if (status == GRB_OPTIMAL) {
  //     std::cout << "Optimal solution found!" << std::endl;

  //     // Print transport names and values
  //     printOptimalTransportVolumes(transport, n);

  //     // Print model constraint names and right-hand sides
  //     printConstraints(model);

  //   } else if (status == GRB_INFEASIBLE) {
  //     std::cout << "Model is infeasible" << std::endl;
  //   } else if (status == GRB_UNBOUNDED) {
  //     std::cout << "Model is unbounded" << std::endl;
  //   } else {
  //     std::cout << "Optimization ended with status " << status << std::endl;
  //   }

  // } catch (GRBException e) {
  //   cout << "Error code = " << e.getErrorCode() << endl;
  //   cout << e.getMessage() << endl;
  // } catch (...) {
  //   cout << "Exception during optimization" << endl;
  // }

  // delete[] transport;

  return 0;
}
