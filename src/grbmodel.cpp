#include "gurobi_c++.h"
#include <fstream>
#include <sstream>
#include <thread>
#include <Rcpp.h>

#ifdef READER_HPP 
  #include "reader.hpp"
#else
  #include "../inst/include/reader.hpp"
#endif
#ifdef UTILS_HPP 
  #include "utils.hpp"
#else
  #include "../inst/include/utils.hpp"
#endif
#ifdef CARTESIAN_HPP
  #include "cartesian.hpp"
#else
  #include "../inst/include/cartesian.hpp"
#endif

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

// Function to stack state index vectors
std::vector<std::vector<int>> stackStateIdxVectors(
    const int& nInventoryLevels, 
    const int& nOrigins, 
    const int& nDestinations) {

    std::vector<std::vector<int>> stateIdx;
    stateIdx.reserve(nOrigins + nDestinations);

    const auto idx = utils::createIndexVector(nInventoryLevels);
    const auto jdx = utils::createIndexVector(2 * nInventoryLevels - 1);

    utils::appendIndexVectors(stateIdx, idx, nOrigins);
    utils::appendIndexVectors(stateIdx, jdx, nDestinations);

    return stateIdx;
}

// Function to stack uncertainty index vectors
std::vector<std::vector<int>> stackUncertaintyIdxVectors(
    const int& nFlowLevels, 
    const int& nOrigins, 
    const int& nDestinations, 
    const int& nSpotRates) {

    std::vector<std::vector<int>> uncertaintyIdx;
    uncertaintyIdx.reserve(nOrigins + nDestinations + 1);

    const auto idx = utils::createIndexVector(nFlowLevels);
    const auto kdx = utils::createIndexVector(nSpotRates);

    utils::appendIndexVectors(uncertaintyIdx, idx, nOrigins + nDestinations);
    uncertaintyIdx.push_back(kdx);

    return uncertaintyIdx;
}

// Function to compute the environment
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
std::vector<std::vector<int>> computeEnvironmentSTL(
    const std::string jsonFile, 
    const std::vector<double>& stateSupport,
    const std::vector<double>& uncertaintySupport) {

  std::ifstream file(jsonFile);
  if (!file) {
    Rcpp::stop("Unable to open file.");
  }

  nlohmann::json input;
  file >> input;

  const int nInventoryLevels = input["nSI"][0];
  const int nFlowLevels = input["nQ"][0];
  const int nSpotRates = input["nW"][0];
  const int nOrigins = input["nI"][0];
  const int nDestinations = input["nJ"][0];
  const int nStrategicCarriers = input["nCS"][0];
  const int nSpotCarriers = input["nCO"][0];
  const int storageLimit = input["R"][0];
  const std::vector<std::vector<int>> Cb = input["Cb"];
  const std::vector<std::vector<int>> Co = input["Co"];
  const int nSources = input["nL_"][0];
  const int nLanes  = input["nL"][0];
  const std::vector<std::vector<int>> bids = input["B"];
  const std::vector<std::vector<int>> lanes = input["L"];
  const std::vector<std::vector<double>> spotRates = input["CTo"];
  const std::vector<std::string> winnerKeys = input["winnerKey"];

  std::unordered_map<std::string, std::vector<int>> winners;
  std::unordered_map<std::string, std::vector<double>> CTb;
  std::unordered_map<std::string, int> carrierIdx;

  winners = importListOfVectors<int>(input["winner"]);
  CTb = importListOfVectors<double>(input["CTb_list"]);
  carrierIdx = importList<int>(input["carrierIdx"]);

  const int n = nSources + nSpotCarriers * nLanes;
  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  const std::vector<int> limits = utils::generateRandomIntegers(nOrigins + nDestinations, 10, 20);

  std::vector<std::vector<int>> stateSupports = stackStateIdxVectors(nInventoryLevels, nOrigins, nDestinations);
  std::vector<std::vector<int>> stateIdx = CartesianProductIntSTL(stateSupports);
  const std::vector<double> extendedStateSupport = utils::mirrorAndNegateVector(stateSupport);

  std::vector<std::vector<int>> uncertaintySupports = stackUncertaintyIdxVectors(nFlowLevels, nOrigins, nDestinations, nSpotRates);
  std::vector<std::vector<int>> uncertaintyIdx = CartesianProductIntSTL(uncertaintySupports);

  std::vector<double> rhs(nStrategicCarriers + nSpotCarriers + nDestinations + nOrigins + 1, 0);

  // Create the Gurobi environment
  GRBEnv env(true);
  // Set the OutputFlag to 0 to turn off logging
  env.set(GRB_IntParam_OutputFlag, 0);

  env.start();

  GRBVar* transport = nullptr;

  // Optimize the model
  try {

    GRBModel model = GRBModel(env);
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_StringAttr_ModelName, "Drayage");

    transport = createTransportVars(model, n, winnerKeys, winners, bids, lanes, CTb, spotRates[0], nSpotCarriers, nLanes);

    addCapacityConstraints(model, transport, winnerKeys,  winners, bids, carrierIdx, nSpotCarriers, nLanes, Cb[0], Co[0]);
    addStorageLimitConstraints(model, transport, winnerKeys, winners, bids, lanes, nSpotCarriers, nLanes, nWarehouses, limits);
    addVolumeConstraint(model, transport, n, 10);

    // Use barrier to solve root relaxation
    model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    model.update();

    // int numConstrs = model.get(GRB_IntAttr_NumConstrs);
    GRBConstr* constraints = model.getConstrs(); // Get all constraints in the model
    if (!constraints) {
      throw std::runtime_error("Failed to get constraints from the model.");
    }

    // Compute environment
    size_t t = 0;

    for (size_t idx = 0; idx < nStrategicCarriers; ++idx) {
      rhs[idx] = Cb[t][idx];
    }

    for (size_t idx = 0; idx < nSpotCarriers; ++idx) {
      rhs[nStrategicCarriers + idx] = Co[t][idx];
    }

    rhs[nStrategicCarriers + nSpotCarriers + nDestinations + nOrigins] = 10;

    for (const auto& sdx : stateIdx) {
      // Uncertainty index
      const auto& udx = uncertaintyIdx.at(0);
      // Update right hand side of constraints

      // Disable automatic model update
      model.set(GRB_IntParam_UpdateMode, 0);

      // Set new RHS values
      for (size_t idx = 0; idx < nDestinations; ++idx) {
        size_t p = nStrategicCarriers + nSpotCarriers + idx;
        rhs[p] = storageLimit - std::max(extendedStateSupport.at(sdx.at(nOrigins + idx)), 0.0);
        constraints[p].set(GRB_DoubleAttr_RHS, rhs[p]);
      }

      for (size_t idx = 0; idx < nOrigins; ++idx) {
        size_t p = nStrategicCarriers + nSpotCarriers + nDestinations + idx;
        rhs[p] = stateSupport.at(sdx.at(idx)) + uncertaintySupport.at(udx.at(idx));
        constraints[p].set(GRB_DoubleAttr_RHS, rhs[p]);
      }

      // Re-enable automatic updates and manually update the model
      model.set(GRB_IntParam_UpdateMode, 1);
      model.update();
      // printConstraints(model);
    }
  
    // Solve
    model.optimize();

    int status = model.get(GRB_IntAttr_Status);

    if (status == GRB_OPTIMAL) {
      // std::cout << "Optimal solution found!" << std::endl;

      // result["objval"] = model.get(GRB_DoubleAttr_ObjVal);

      // Rcpp::NumericVector x(n);
      // for (size_t i = 0; i < n; ++i) {
      //   x[i] = transport[i].get(GRB_DoubleAttr_X);
      // }
      // result["x"] = x;

    } else if (status == GRB_INFEASIBLE) {
      std::cout << "Model is infeasible" << std::endl;
    } else if (status == GRB_UNBOUNDED) {
      std::cout << "Model is unbounded" << std::endl;
    } else {
      std::cout << "Optimization ended with status " << status << std::endl;
    }

    delete [] transport;

    return std::vector<std::vector<int>>();

  } catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }

  return std::vector<std::vector<int>>();
}

// Optimize the model and reurn objective and decision in R list
//' @useDynLib TLPR
//' @export
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

    addCapacityConstraints(model, transport, winnerKeys,  winners, bids, carrierIdx, nCO, nL, Cb[0], Co[0]);
    addStorageLimitConstraints(model, transport, winnerKeys, winners, bids, lanes, nCO, nL, nWarehouses, limits);
    addVolumeConstraint(model, transport, n, 10);

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
  // std::unordered_map<std::string, std::vector<int>> Cb;
  std::unordered_map<std::string, int> carrierIdx;

  const int n = nL_ + nCO * nL;
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
  GRBEnv env;
  GRBVar* transport = nullptr;

  try {
    GRBModel model = GRBModel(env);
    model.set(GRB_StringAttr_ModelName, "Drayage");

    transport = createTransportVars(model, n, winnerKeys, winners, bids, lanes, CTb, spotRates[0], nCO, nL);

    addCapacityConstraints(model, transport, winnerKeys, winners, bids, carrierIdx, nCO, nL, Cb[0], Co[0]);
    addStorageLimitConstraints(model, transport, winnerKeys, winners, bids, lanes, nCO, nL, nWarehouses, limits);
    addVolumeConstraint(model, transport, n, 10);

    // Use barrier to solve root relaxation
    model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    model.update();

    const int nInventoryLevels = input["nSI"][0];
    std::vector<double> stateSupport(nInventoryLevels);
    for (size_t i = 0; i < nInventoryLevels; ++i) {
      stateSupport[i] = i;
    }
    const std::vector<double> uncertaintySupport = input["Q"]["vals"].get<std::vector<double>>();

    // std::vector<std::vector<int>> stateSupports = computeEnvironmentSTL(argv[1], stateSupport, uncertaintySupport);

    // // Solve
    // model.optimize();

    // int status = model.get(GRB_IntAttr_Status);

    // if (status == GRB_OPTIMAL) {
    //   std::cout << "Optimal solution found!" << std::endl;

    //   // Print transport names and values
    //   // printOptimalTransportVolumes(transport, n);

    //   // Print model constraint names and right-hand sides
    //   printConstraints(model);

    // } else if (status == GRB_INFEASIBLE) {
    //   std::cout << "Model is infeasible" << std::endl;
    // } else if (status == GRB_UNBOUNDED) {
    //   std::cout << "Model is unbounded" << std::endl;
    // } else {
    //   std::cout << "Optimization ended with status " << status << std::endl;
    // }

  } catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }


  delete[] transport;

  return 0;
}
