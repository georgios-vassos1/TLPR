#include "gurobi_c++.h"
#include <fstream>
#include <sstream>
#include <thread>
#include <Rcpp.h>
#include <chrono>
#include <omp.h>

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
  try {
    // Get the number of constraints
    int numConstraints = model.get(GRB_IntAttr_NumConstrs);

    // Get the constraints
    GRBConstr* constraints = model.getConstrs();

    // Print each constraint's name and RHS value
    std::cout << "Constraints: " << std::endl;
    for (int i = 0; i < numConstraints; ++i) {
      std::cout << constraints[i].get(GRB_StringAttr_ConstrName) << " = " 
                << constraints[i].get(GRB_DoubleAttr_RHS) << std::endl;
    }

    // Clean up dynamically allocated array
    delete[] constraints;

  } catch (GRBException& e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch (...) {
    std::cout << "Exception during optimization" << std::endl;
  }
}

void printObjectiveVector(GRBModel& model) {
  try {
    // Get the number of decision variables
    int numVars = model.get(GRB_IntAttr_NumVars);

    std::cout << "Objective Vector: " << std::endl;
    for (int i = 0; i < numVars; ++i) {
      GRBVar v = model.getVar(i);
      double objCoeff = v.get(GRB_DoubleAttr_Obj);

      std::cout << v.get(GRB_StringAttr_VarName) << ": " << objCoeff << std::endl;
    }

  } catch (GRBException& e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch (...) {
    std::cout << "Exception during optimization" << std::endl;
  }
}

// Function to print the optimal transport volumes from gurobi model optimization
void printOptimalTransportVolumes(GRBVar* transport, const int& n) {

  std::cout << "Transport volumes:" << std::endl;
  for (int i = 0; i < n; ++i) {
    try {
      std::cout << transport[i].get(GRB_StringAttr_VarName) << " = " << transport[i].get(GRB_DoubleAttr_X) << std::endl;
    } catch (GRBException& e) {
      std::cerr << "Error code = " << e.getErrorCode() << std::endl;
      std::cerr << e.getMessage() << std::endl;
    }
  }

}

void updateSpotRates(
    GRBModel& model,
    GRBVar* transport,
    const std::vector<double>& spotRates,
    const int& nStrategicSources,
    const int& nSpotCarriers,
    const int& nLanes) {

  int idx = nStrategicSources;
  // Spot carrier variables
  for (int carrierIndex = 0; carrierIndex < nSpotCarriers; ++carrierIndex) {
    for (int laneIndex = 0; laneIndex < nLanes; ++laneIndex) {
      transport[idx++].set(GRB_DoubleAttr_Obj, spotRates[carrierIndex * nLanes + laneIndex]);
    }
  }

  model.update();
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
  for (int carrierIndex = 0; carrierIndex < nSpotCarriers; carrierIndex++) {

    for (int laneIndex = 0; laneIndex < nLanes; laneIndex++) {
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
  for (int i = 0; i < n; i++) {
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
  for (int spotCarrierIdx = 0; spotCarrierIdx < nSpotCarriers; ++spotCarrierIdx) {
    const int capacity = spotCaps.at(spotCarrierIdx);

    GRBLinExpr expr = 0;
    for (int laneIndex = 0; laneIndex < nLanes; ++laneIndex) {
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

    for (int i = 0; i < nWarehouses[m]; ++i) {
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
      for (int carrierIndex = 0; carrierIndex < nSpotCarriers; ++carrierIndex) {

        for (int laneIndex = 0; laneIndex < nLanes; ++laneIndex) {
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
  // utils::appendIndexVectors(uncertaintyIdx, kdx, nSpotCarriers);
  uncertaintyIdx.push_back(kdx);

  return uncertaintyIdx;
}

// Function to update the state index
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
std::vector<std::vector<int>> updateStateIdx(
    const std::vector<int>& stateIdx, 
    const std::vector<int>& inflowIdx, 
    const std::vector<std::vector<int>>& outflowIndices,
    const std::vector<double>& stateSupport, 
    const std::vector<double>& extendedStateSupport, 
    const std::vector<double>& flowSupport,
    const std::vector<double>& xI,
    const std::vector<double>& xJ,
    const double& storageLimit,
    const std::vector<int>& stateKey,
    const std::vector<int>& flowKey,
    const int& nOrigins,
    const int& nDestinations) {

  std::vector<std::vector<int>> nextStateIdx(outflowIndices.size(), std::vector<int>(2, 0));
  std::vector<double> nextState(nOrigins + nDestinations, 0.0);
  std::vector<double> fdx(nOrigins + nDestinations, 0.0);

  // Calculate initial state for origins
  for (int i = 0; i < nOrigins; ++i) {
    fdx[i] = inflowIdx[i];
    nextState[i] = std::max(std::min(stateSupport[stateIdx[i]] + flowSupport[fdx[i]], storageLimit) - xI[i], 0.0);
  }

  // Iterate over flow indices
  for (size_t k = 0; k < outflowIndices.size(); ++k) {
    std::vector<int> outflowIdx = outflowIndices[k];
    // Update state for destinations based on current flow index
    for (int j = 0; j < nDestinations; ++j) {
      fdx[nOrigins + j] = outflowIdx[j];
      nextState[nOrigins + j] = storageLimit + std::min(std::max(extendedStateSupport[stateIdx[nOrigins + j]] - flowSupport[fdx[nOrigins + j]], -storageLimit) + xJ[j], storageLimit);
    }
    // Calculate the index using inner product
    nextStateIdx[k][0] = std::inner_product(nextState.begin(), nextState.end(), stateKey.begin(), 0);
    // Compute and return scenario index (easy)
    nextStateIdx[k][1] = std::inner_product(fdx.begin(), fdx.end(), flowKey.begin(), 0);
  }

  return nextStateIdx;
}


// Function to compute the environment
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd computeEnvironmentCx(
    const std::string jsonFile, 
    const std::vector<double>& stateSupport,
    const std::vector<double>& flowSupport,
    int numThreads = 8) {

  // Read the JSON file
  std::ifstream file(jsonFile);
  if (!file) {
    Rcpp::stop("Unable to open file.");
  }

  nlohmann::json input;
  file >> input;

  /*
   *  Extract the data from the JSON object
   */
  // Primitive types
  const int nInventoryLevels   = input["R"][0].get<int>() + 1;
  const int nFlowLevels        = input["nQ"][0].get<int>();
  const int nSpotRates         = input["nW"][0].get<int>();
  const int nOrigins           = input["nI"][0].get<int>();
  const int nDestinations      = input["nJ"][0].get<int>();
  const int nStrategicCarriers = input["nCS"][0].get<int>();
  const int nSpotCarriers      = input["nCO"][0].get<int>();
  const int nSources           = input["nL_"][0].get<int>();
  const int nLanes             = input["nL"][0].get<int>();
  const int nActions           = input["R"][0].get<int>() + 1;
  const int storageLimit       = input["R"][0].get<int>();
  const int tau                = input["tau"][0].get<int>();

  // Compound types (data container classes)
  const std::vector<std::vector<int>>& fromOrigin    = input["from_i"];
  const std::vector<std::vector<int>>& toDestination = input["to_j"];
  const std::vector<std::vector<int>>& Cb    = input["Cb"];
  const std::vector<std::vector<int>>& Co    = input["Co"];
  const std::vector<std::vector<int>>& bids  = input["B"];
  const std::vector<std::vector<int>>& lanes = input["L"];
  const std::vector<std::string>& winnerKeys = input["winnerKey"];
  const std::vector<std::vector<double>>& spotRates = input["CTo"];
  const std::vector<double>& spotRateSupport = input["W"]["vals"].get<std::vector<double>>();

  // State and uncertainty cipher keys
  const std::vector<int>& stateKeys = input["stateKeys"];
  const std::vector<int>& flowKeys  = input["flowKeys"];

  // Converting R to C++ types
  std::unordered_map<std::string, std::vector<int>>    winners;
  std::unordered_map<std::string, std::vector<double>> CTb;
  std::unordered_map<std::string, int> carrierIdx;

  winners    = importListOfVectors<int>(input["winner"]);
  CTb        = importListOfVectors<double>(input["CTb_list"]);
  carrierIdx = importList<int>(input["carrierIdx"]);

  // Auxiliary variables
  const int n = nSources + nSpotCarriers * nLanes;

  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  const std::vector<int> limits      = utils::generateRandomIntegers(nOrigins + nDestinations, 10, 20);

  // Generate extended state support data (positive and negative states)
  const std::vector<double> extendedStateSupport = utils::mirrorAndNegateVector(stateSupport);

  // Generate state support data (only positive states)
  std::vector<std::vector<int>> stateSupportStack = stackStateIdxVectors(nInventoryLevels, nOrigins, nDestinations);
  // Generate state index data
  std::vector<std::vector<int>> stateIndices = CartesianProductIntSTL(stateSupportStack);

  // Generate flow index data
  std::vector<int> flowIdxSingle = utils::createIndexVector(nFlowLevels);

  // Inflow indices
  std::vector<std::vector<int>> inflowIdxStack;
  utils::appendIndexVectors(inflowIdxStack, flowIdxSingle, nOrigins);
  std::vector<std::vector<int>> inflowIndices = CartesianProductIntSTL(inflowIdxStack);

  // Outflow indices
  std::vector<std::vector<int>> outflowIdxStack;
  utils::appendIndexVectors(outflowIdxStack, flowIdxSingle, nDestinations);
  std::vector<std::vector<int>> outflowIndices = CartesianProductIntSTL(outflowIdxStack);

  // Generate spot rate index data
  std::vector<int> spotRateIdx = utils::createIndexVector(nSpotRates);

  // Spot rate indices
  std::vector<std::vector<int>> spotRateIdxStack;
  utils::appendIndexVectors(spotRateIdxStack, spotRateIdx, nSpotCarriers);
  std::vector<std::vector<int>> spotRateIndices = CartesianProductIntSTL(spotRateIdxStack);

  // TODO: Get state and scenario keys (maybe C++ utility)

  const int nSdx = stateIndices.size();
  const int nAdx = nActions;
  const int nQdx = inflowIndices.size();
  const int nDdx = outflowIndices.size();
  const int nWdx = spotRateIndices.size();
  const int nScen = nQdx * nDdx * nWdx;

  // Start timing
  auto start = std::chrono::high_resolution_clock::now();

  // Initialize the Eigen matrix with n rows and 5 columns
  Eigen::MatrixXd transit(tau * nSdx * nAdx * nScen, 6);

  // std::cout << "Number of states: " << nSdx << std::endl;
  // std::cout << "Number of actions: " << nAdx << std::endl;
  // std::cout << "Number of inflow indices: " << nQdx << std::endl;
  // std::cout << "Number of outflow indices: " << nDdx << std::endl;
  // std::cout << "Number of spot rates: " << nWdx << std::endl;
  // std::cout << "Number of rows: " << nSdx * nAdx * nScen << std::endl;

  // Set the number of threads
  omp_set_num_threads(numThreads);

  for (int t = 0; t < tau; ++t) {

    // TODO: fix parallelization
    // Enclose the thread-specific code within a block
    #pragma omp parallel
    {
      // Create the Gurobi environment
      GRBEnv env(true);
      // Set the OutputFlag to 0 to turn off logging
      env.set(GRB_IntParam_OutputFlag, 0);
      env.start();

      // Create the Gurobi model
      GRBModel model = GRBModel(env);
      // Set the OutputFlag to 0 to turn off logging
      model.set(GRB_IntParam_OutputFlag, 0);
      model.set(GRB_StringAttr_ModelName, "Drayage");

      // Create the decision variables
      GRBVar* transport = nullptr;
      // Create the constraints
      GRBConstr* constraints = nullptr;
      // Initialize vector for the right hand side (RHS) of constraints
      std::vector<double> rhs(nStrategicCarriers + nSpotCarriers + nDestinations + nOrigins + 1, 0);

      try {

        // Instantiate dummy Gurobi model once to get the decision variables and contraints in place
        transport = createTransportVars(model, n, winnerKeys, winners, bids, lanes, CTb, spotRates[t], nSpotCarriers, nLanes);

        // Do not change the order of calling the constraints
        addCapacityConstraints(model, transport, winnerKeys,  winners, bids, carrierIdx, nSpotCarriers, nLanes, Cb[t], Co[t]);
        addStorageLimitConstraints(model, transport, winnerKeys, winners, bids, lanes, nSpotCarriers, nLanes, nWarehouses, limits);
        addVolumeConstraint(model, transport, n, 10);

        // Use barrier to solve root relaxation
        model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

        model.update();

        // int numConstrs = model.get(GRB_IntAttr_NumConstrs);
        constraints = model.getConstrs(); // Get all constraints in the model
        if (!constraints) {
          throw std::runtime_error("Failed to get constraints from the model.");
        }

        // Compute environment
        for (int idx = 0; idx < nStrategicCarriers; ++idx) {
          rhs[idx] = Cb[t][idx];
        }

        for (int idx = 0; idx < nSpotCarriers; ++idx) {
          rhs[nStrategicCarriers + idx] = Co[t][idx];
        }

        // Compute the environment
        #pragma omp for
        for (int i = 0; i < nSdx; ++i) {
          const auto& stateIdx = stateIndices[i];

          // Update right hand side of constraints

          // Disable automatic model update
          model.set(GRB_IntParam_UpdateMode, 0);

          // Set new RHS values
          for (int idx = 0; idx < nDestinations; ++idx) {
            int p = nStrategicCarriers + nSpotCarriers + nOrigins + idx;

            rhs[p] = storageLimit - std::max(extendedStateSupport[stateIdx[nOrigins + idx]], 0.0);
            constraints[p].set(GRB_DoubleAttr_RHS, rhs[p]);
          }

          for (int j = 0; j < nAdx; ++j) {
            // Disable automatic model update
            model.set(GRB_IntParam_UpdateMode, 0);

            int p = nStrategicCarriers + nSpotCarriers + nOrigins + nDestinations;

            constraints[p].set(GRB_DoubleAttr_RHS, j);

            // Uncertainty index
            for (int k1 = 0; k1 < nQdx; ++k1) {
              const auto& inflowIdx = inflowIndices[k1];

              // Disable automatic model update
              model.set(GRB_IntParam_UpdateMode, 0);

              for (int k1dx = 0; k1dx < nOrigins; ++k1dx) {
                int p = nStrategicCarriers + nSpotCarriers + k1dx;

                rhs[p] = stateSupport[stateIdx[k1dx]] + flowSupport[inflowIdx[k1dx]];
                constraints[p].set(GRB_DoubleAttr_RHS, rhs[p]);
              }

              // Re-enable automatic updates and manually update the model
              model.set(GRB_IntParam_UpdateMode, 1);
              model.update();
              // printConstraints(model);

              for (int k3 = 0; k3 < nWdx; ++k3) {
                const auto& spotRateIdx = spotRateIndices[k3];

                std::vector<double> spotRatesTmp(nSpotCarriers, 0.0);
                for (int k3dx = 0; k3dx < nSpotCarriers; ++k3dx) {
                    spotRatesTmp[k3dx] = spotRateSupport[spotRateIdx[k3dx]];
                }
                updateSpotRates(model, transport, spotRatesTmp, nStrategicCarriers, nSpotCarriers, nLanes);

                // Solve
                model.optimize();

                int status = model.get(GRB_IntAttr_Status);

                if (status == GRB_OPTIMAL) {
                  // std::cout << "Optimal solution found!" << std::endl;

                  // Get the optimal solution
                  double objval = model.get(GRB_DoubleAttr_ObjVal);

                  // Get the optimal allocation vector
                  std::vector<double> x(n);
                  for (int k1dx = 0; k1dx < n; ++k1dx) {
                    x[k1dx] = transport[k1dx].get(GRB_DoubleAttr_X);
                  }

                  // Compute the total outflow from each origin and the total inflow to each destination
                  std::vector<double> xI(nOrigins, 0.0);
                  for (int k1dx = 0; k1dx < nOrigins; ++k1dx) {
                    for (const auto& ldx : fromOrigin[k1dx]) {
                      xI[k1dx] += x[ldx - 1];
                    }
                  }

                  std::vector<double> xJ(nDestinations, 0.0);
                  for (int k1dx = 0; k1dx < nDestinations; ++k1dx) {
                    for (const auto& ldx : toDestination[k1dx]) {
                      xJ[k1dx] += x[ldx - 1];
                    }
                  }

                  // Compute and store the next state
                  // Get next state code index
                  std::vector<std::vector<int>> nextStateIndices(nDdx, std::vector<int>(2, 0));

                  // Iterate over flow indices
                  for (int k2 = 0; k2 < nDdx; ++k2) {
                    std::vector<int> outflowIdx = outflowIndices[k2];
                    // Update state for destinations based on current flow index
                    for (int k2dx = 0; k2dx < nDestinations; ++k2dx) {
                    }

                    // Store (newStateIdx, objective, stateIdx, actionIdx, scenarioIdx) in transit matrix
                    int kdx = (k3 * nQdx + k2) * nDdx + k1;
                    int p = ((t * nSdx + i) * nAdx + j) * nQdx * nDdx * nWdx + kdx;

                    transit(p, 0) = 0;
                    transit(p, 1) = objval;
                    transit(p, 2) = i + 1;
                    transit(p, 3) = j + 1;
                    transit(p, 4) = kdx + 1;
                  }

                }
              }

            }

          }
        }

        delete [] transport;
        delete [] constraints;

      } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Exception during optimization" << endl;
      }
    }
  }

  // End timing
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;

  std::cout << "Time taken: " << elapsed.count() << " seconds" << std::endl;

  return transit;
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
      for (int i = 0; i < n; ++i) {
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
  // const int nFlowLevels = input["nQ"][0];
  const std::vector<double> flowSupport  = input["Q"]["vals"].get<std::vector<double>>();
  const int nOrigins = input["nI"][0];
  const int nDestinations = input["nJ"][0];
  const int nStrategicSources = input["nL_"][0];
  const int nSpotCarriers = input["nCO"][0];
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

  const int n = nStrategicSources + nSpotCarriers * nL;
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

  // std::vector<int> stateIdx = {4, 5, 1, 14, 16};
  // std::vector<int> flowIdx = {0, 0, 0};

  // std::vector<int> flowIdxSingle = utils::createIndexVector(nFlowLevels);
  // std::vector<std::vector<int>> outflowIdxStack;
  // utils::appendIndexVectors(outflowIdxStack, flowIdxSingle, nDestinations);
  // std::vector<std::vector<int>> outflowIndices = CartesianProductIntSTL(outflowIdxStack);

  // std::vector<double> stateSupport = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  // std::vector<double> extendedStateSupport = utils::mirrorAndNegateVector(stateSupport);
  // std::vector<double> xI = {0.0, 0.0, 0.0};
  // std::vector<double> xJ = {0.0, 0.0};
  // double storageLimit = 10.0;
  // std::vector<int> stateKey = {1, 11, 121, 1331, 27951};
  // std::vector<int> flowKey = {1, 3, 9, 27, 81};

  // std::vector<std::vector<int>> nextState = updateStateIdx(stateIdx, flowIdx, outflowIndices, stateSupport, extendedStateSupport, flowSupport, xI, xJ, storageLimit, stateKey, flowKey, nOrigins, nDestinations);
  // // Print next state
  // std::cout << "Next state: " << std::endl;
  // for (const auto& ns : nextState) {
  //   std::copy(ns.begin(), ns.end(), std::ostream_iterator<int>(std::cout, " "));
  //   std::cout << std::endl;
  // }

  // Print the strategy of the instance
  // printInstance(winners, bids, lanes, CTb, spotRates[0], nSpotCarriers, nL);

  // Create the Gurobi environment
  GRBEnv env;
  GRBVar* transport = nullptr;

  try {
    GRBModel model = GRBModel(env);
    model.set(GRB_StringAttr_ModelName, "Drayage");

    transport = createTransportVars(model, n, winnerKeys, winners, bids, lanes, CTb, spotRates[0], nSpotCarriers, nL);

    addCapacityConstraints(model, transport, winnerKeys, winners, bids, carrierIdx, nSpotCarriers, nL, Cb[0], Co[0]);
    addStorageLimitConstraints(model, transport, winnerKeys, winners, bids, lanes, nSpotCarriers, nL, nWarehouses, limits);
    addVolumeConstraint(model, transport, n, 10);

    // Use barrier to solve root relaxation
    model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    model.update();

    // printObjectiveVector(model);
    // updateSpotRates(model, transport, spotRates[1], nStrategicSources, nSpotCarriers, nL);
    // printObjectiveVector(model);

    // const int nInventoryLevels = input["nSI"][0];
    // std::vector<double> stateSupport(nInventoryLevels);
    // for (size_t i = 0; i < nInventoryLevels; ++i) {
    //   stateSupport[i] = i;
    // }

    // This is false assumption (carefully check the state support)
    // const std::vector<double> uncertaintySupport = input["Q"]["vals"].get<std::vector<double>>();

    // // Solve
    // model.optimize();

    // int status = model.get(GRB_IntAttr_Status);

    // if (status == GRB_OPTIMAL) {
    //   std::cout << "Optimal solution found!" << std::endl;

    //   // Print transport names and values
    //   // printOptimalTransportVolumes(transport, n);

    //   // Print model constraint names and right-hand sides
    //   // printConstraints(model);

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
