#include <Highs.h>
#include <atomic>
#include <fstream>
#include <sstream>
#include <thread>
#include <Rcpp.h>
#include <chrono>
#include <omp.h>
#include <unistd.h>
#include <fcntl.h>

#ifdef STANDALONE_BUILD
  #include "reader.hpp"
#else
  #include "../inst/include/reader.hpp"
#endif
#ifdef STANDALONE_BUILD
  #include "utils.hpp"
#else
  #include "../inst/include/utils.hpp"
#endif
#ifdef STANDALONE_BUILD
  #include "cartesian.hpp"
#else
  #include "../inst/include/cartesian.hpp"
#endif

// Explicit using-declarations instead of 'using namespace std' to avoid
// pulling in the entire standard namespace into library translation units.
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::ifstream;

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

// Function to print the objective vector in a HiGHS model
void printObjectiveVector(Highs& highs) {
  const HighsLp& lp = highs.getLp();
  int numVars = static_cast<int>(lp.col_cost_.size());
  std::cout << "Number of decision variables: " << numVars << std::endl;
  std::cout << "Objective Vector: " << std::endl;
  for (int i = 0; i < numVars; ++i) {
    std::cout << "x[" << i << "]: " << lp.col_cost_[i] << std::endl;
  }
}

// Function to print the constraints in a HiGHS model
void printConstraints(Highs& highs) {
  const HighsLp& lp = highs.getLp();
  int numConstraints = static_cast<int>(lp.row_lower_.size());
  std::cout << "Constraints: " << std::endl;
  for (int i = 0; i < numConstraints; ++i) {
    std::cout << "row[" << i << "]: lower=" << lp.row_lower_[i]
              << " upper=" << lp.row_upper_[i] << std::endl;
  }
}

// Function to print the optimal transport volumes
void printOptimalTransportVolumes(const Highs& highs, const std::vector<int>& colIdx, const int& n) {
  const auto& sol = highs.getSolution();
  std::cout << "Transport volumes:" << std::endl;
  for (int i = 0; i < n; ++i) {
    std::cout << "x[" << colIdx[i] << "] = " << sol.col_value[colIdx[i]] << std::endl;
  }
}

// Function to update spot rates in the objective vector
void updateSpotRates(
    Highs& highs,
    const std::vector<int>& colIdx,
    const std::vector<double>& spotRates,
    const int& nStrategicSources,
    const int& nSpotCarriers,
    const int& nLanes) {

  int idx = nStrategicSources;
  for (int carrierIndex = 0; carrierIndex < nSpotCarriers; ++carrierIndex) {
    for (int laneIndex = 0; laneIndex < nLanes; ++laneIndex) {
      highs.changeColCost(colIdx[idx++], spotRates[carrierIndex * nLanes + laneIndex]);
    }
  }
}

// Expand a leading '~' to the user's home directory.
// std::ifstream does not perform shell tilde-expansion.
static std::string expand_path(const std::string& path) {
  if (!path.empty() && path[0] == '~') {
    const char* home = std::getenv("HOME");
    if (home) return std::string(home) + path.substr(1);
  }
  return path;
}

// Function to create decision variables for transportation
std::vector<int> createTransportVars(
    Highs& highs,
    const int& n,
    const std::vector<std::string>& winnerKeys,
    const std::unordered_map<std::string, std::vector<int>>& winners,
    const std::vector<std::vector<int>>& bids,
    const std::vector<std::vector<int>>& lanes,
    const std::unordered_map<std::string, std::vector<double>>& contractRates,
    const std::vector<double>& spotRates,
    const int& nSpotCarriers,
    const int& nLanes) {

  std::vector<int> colIdx;
  colIdx.reserve(n);

  // Strategic carrier variables
  for (const auto& winnerKey : winnerKeys) {
    int k = 0;
    for (size_t bidIndex : winners.at(winnerKey)) {
      const auto& bid = bids[bidIndex - 1];
      for (size_t laneIndex : bid) {
        int col = highs.getNumCol();
        highs.addVar(0.0, kHighsInf);
        highs.changeColCost(col, contractRates.at(winnerKey).at(k++));
        colIdx.push_back(col);
      }
    }
  }

  // Spot carrier variables
  for (int carrierIndex = 0; carrierIndex < nSpotCarriers; ++carrierIndex) {
    for (int laneIndex = 0; laneIndex < nLanes; ++laneIndex) {
      int col = highs.getNumCol();
      highs.addVar(0.0, kHighsInf);
      highs.changeColCost(col, spotRates[carrierIndex * nLanes + laneIndex]);
      colIdx.push_back(col);
    }
  }

  return colIdx;
}

void addCapacityConstraints(
    Highs& highs,
    const std::vector<int>& colIdx,
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
  for (const auto& winnerKey : winnerKeys) {
    const int capacity = strategicCaps.at(carrierIdx.at(winnerKey) - 1);
    std::vector<int>    ids;
    std::vector<double> vals;
    for (size_t bidIndex : winners.at(winnerKey)) {
      const auto& bid = bids[bidIndex - 1];
      for (size_t i = 0; i < bid.size(); ++i) {
        ids.push_back(colIdx[k++]);
        vals.push_back(1.0);
      }
    }
    highs.addRow(-kHighsInf, capacity, static_cast<int>(ids.size()), ids.data(), vals.data());
  }

  // Spot carrier capacity constraints
  for (int spotCarrierIdx = 0; spotCarrierIdx < nSpotCarriers; ++spotCarrierIdx) {
    const int capacity = spotCaps.at(spotCarrierIdx);
    std::vector<int>    ids;
    std::vector<double> vals;
    for (int laneIndex = 0; laneIndex < nLanes; ++laneIndex) {
      ids.push_back(colIdx[k++]);
      vals.push_back(1.0);
    }
    highs.addRow(-kHighsInf, capacity, static_cast<int>(ids.size()), ids.data(), vals.data());
  }
}

void addStorageLimitConstraints(
    Highs& highs,
    const std::vector<int>& colIdx,
    const std::vector<std::string>& winnerKeys,
    const std::unordered_map<std::string, std::vector<int>>& winners,
    const std::vector<std::vector<int>>& bids,
    const std::vector<std::vector<int>>& lanes,
    const int& nSpotCarriers,
    const int& nLanes,
    const std::vector<int>& nWarehouses,
    const std::vector<int>& limits) {

  int position = 0;
  for (size_t m = 0; m < nWarehouses.size(); ++m) {
    for (int i = 0; i < nWarehouses[m]; ++i) {
      std::vector<int>    ids;
      std::vector<double> vals;
      int k = 0;

      // Strategic carrier contributions
      for (const auto& winnerKey : winnerKeys) {
        for (size_t bidIndex : winners.at(winnerKey)) {
          const auto& bid = bids[bidIndex - 1];
          for (size_t laneIndex : bid) {
            const auto& lane = lanes[laneIndex - 1];
            if (lane[m] == i + 1) {
              ids.push_back(colIdx[k]);
              vals.push_back(1.0);
            }
            k++;
          }
        }
      }

      // Spot carrier contributions
      for (int carrierIndex = 0; carrierIndex < nSpotCarriers; ++carrierIndex) {
        for (int laneIndex = 0; laneIndex < nLanes; ++laneIndex) {
          const auto& lane = lanes[laneIndex];
          if (lane[m] == i + 1) {
            ids.push_back(colIdx[k]);
            vals.push_back(1.0);
          }
          k++;
        }
      }

      highs.addRow(-kHighsInf, limits[position + i],
                   static_cast<int>(ids.size()), ids.data(), vals.data());
    }
    position += nWarehouses[m];
  }
}

// Volume constraint
void addVolumeConstraint(
    Highs& highs,
    const std::vector<int>& colIdx,
    const int& n,
    const double& At) {

  std::vector<double> vals(n, 1.0);
  highs.addRow(At, At, n, colIdx.data(), vals.data());
}

// Function to compute the environment
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd computeEnvironmentCx(
    const std::string jsonFile,
    const int& t,
    const std::vector<double>& stateSupport,
    const std::vector<double>& flowSupport,
    int numThreads = 8) {

  // Read the JSON file
  std::ifstream file(expand_path(jsonFile));
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
  const int nServices          = input["nL_"][0].get<int>();
  const int nLanes             = input["nL"][0].get<int>();
  const int nActions           = input["R"][0].get<int>() + 1;
  const int storageLimit       = input["R"][0].get<int>();

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
  const int n = nServices + nSpotCarriers * nLanes;

  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  const std::vector<int> limits(nOrigins + nDestinations, storageLimit);

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

  const int nSdx = static_cast<int>(stateIndices.size());
  const int nAdx = nActions;
  const int nQdx = static_cast<int>(inflowIndices.size());
  const int nDdx = static_cast<int>(outflowIndices.size());
  const int nWdx = static_cast<int>(spotRateIndices.size());
  const long long nScen = static_cast<long long>(nQdx) * nDdx * nWdx;
  const long long nTransit = static_cast<long long>(nSdx) * nAdx * nScen;

  // Start timing
  auto start = std::chrono::high_resolution_clock::now();

  // Initialize the Eigen matrix with n rows and 6 columns
  Eigen::MatrixXd transit(nTransit, 6);
  transit.setConstant(std::numeric_limits<double>::quiet_NaN());

  // Build the base model once, serially.
  Highs baseHighs;
  baseHighs.setOptionValue("output_flag", false);
  baseHighs.setOptionValue("threads", 1);

  std::vector<int> baseColIdx = createTransportVars(baseHighs, n, winnerKeys, winners, bids, lanes, CTb, spotRates[t], nSpotCarriers, nLanes);
  try {
    addCapacityConstraints(baseHighs, baseColIdx, winnerKeys, winners, bids, carrierIdx, nSpotCarriers, nLanes, Cb[t], Co[t]);
    addStorageLimitConstraints(baseHighs, baseColIdx, winnerKeys, winners, bids, lanes, nSpotCarriers, nLanes, nWarehouses, limits);
    addVolumeConstraint(baseHighs, baseColIdx, n, 0); // placeholder; overridden per action in loop
  } catch (const std::exception& e) {
    Rcpp::stop("Error building base model: %s", e.what());
  }

  // Extract LP for OMP cloning
  const HighsLp baseLp = baseHighs.getLp();

  // Capacity RHS values are constant across all states and actions within t.
  std::vector<double> baseRhs(nStrategicCarriers + nSpotCarriers + nDestinations + nOrigins + 1, 0.0);
  for (int idx = 0; idx < nStrategicCarriers; ++idx)
    baseRhs[idx] = Cb[t][idx];
  for (int idx = 0; idx < nSpotCarriers; ++idx)
    baseRhs[nStrategicCarriers + idx] = Co[t][idx];

  // Error reporting across OMP threads
  std::atomic<bool> anyThreadFailed{false};
  std::string       threadErrorMsg;

  #pragma omp parallel num_threads(numThreads)
  {
    Highs threadHighs;
    threadHighs.setOptionValue("output_flag", false);
    threadHighs.setOptionValue("threads", 1);
    threadHighs.passModel(baseLp);
    const std::vector<int>& colIdx = baseColIdx;

    // Thread-local working storage
    std::vector<double> rhs           = baseRhs;
    std::vector<int>    fdx           (nOrigins + nDestinations + nSpotCarriers, 0);
    std::vector<double> nextState     (nOrigins + nDestinations, 0.0);
    std::vector<double> spotRatesTmp  (nSpotCarriers * nLanes, 0.0);
    std::vector<double> x             (n, 0.0);
    std::vector<double> xI            (nOrigins, 0.0);
    std::vector<double> xJ            (nDestinations, 0.0);
    int                 nextStateIdx = 0;

    const int volumeRow = nStrategicCarriers + nSpotCarriers + nOrigins + nDestinations;

    try {
      #pragma omp for schedule(dynamic)
      for (int i = 0; i < nSdx; ++i) {
        const auto& stateIdx = stateIndices[i];

        for (int idx = 0; idx < nDestinations; ++idx) {
          int p = nStrategicCarriers + nSpotCarriers + nOrigins + idx;
          rhs[p] = storageLimit - extendedStateSupport[stateIdx[nOrigins + idx]];
          threadHighs.changeRowBounds(p, -kHighsInf, rhs[p]);
        }

        for (int j = 0; j < nAdx; ++j) {
          threadHighs.changeRowBounds(volumeRow, static_cast<double>(j), static_cast<double>(j));

          for (int k1 = 0; k1 < nQdx; ++k1) {
            const auto& inflowIdx = inflowIndices[k1];

            for (int k1dx = 0; k1dx < nOrigins; ++k1dx) {
              int p = nStrategicCarriers + nSpotCarriers + k1dx;
              rhs[p] = stateSupport[stateIdx[k1dx]] + flowSupport[inflowIdx[k1dx]];
              threadHighs.changeRowBounds(p, -kHighsInf, rhs[p]);
            }

            std::fill(fdx.begin(), fdx.end(), 0);

            for (int k3 = 0; k3 < nWdx; ++k3) {
              const auto& spotRateIdx = spotRateIndices[k3];

              std::fill(spotRatesTmp.begin(), spotRatesTmp.end(), 0.0);
              for (int k3dx = 0; k3dx < nSpotCarriers; ++k3dx) {
                fdx[nOrigins + nDestinations + k3dx] = spotRateIdx[k3dx];
                for (int ldx = 0; ldx < nLanes; ++ldx)
                  spotRatesTmp[k3dx * nLanes + ldx] = spotRateSupport[spotRateIdx[k3dx]];
              }
              updateSpotRates(threadHighs, colIdx, spotRatesTmp, nServices, nSpotCarriers, nLanes);

              threadHighs.run();

              if (threadHighs.getModelStatus() == HighsModelStatus::kOptimal) {
                double objval = threadHighs.getObjectiveValue();
                const auto& sol = threadHighs.getSolution();

                for (int k1dx = 0; k1dx < n; ++k1dx)
                  x[k1dx] = sol.col_value[colIdx[k1dx]];

                std::fill(xI.begin(), xI.end(), 0.0);
                for (int k1dx = 0; k1dx < nOrigins; ++k1dx)
                  for (const auto& ldx : fromOrigin[k1dx])
                    xI[k1dx] += x[ldx - 1];

                std::fill(xJ.begin(), xJ.end(), 0.0);
                for (int k1dx = 0; k1dx < nDestinations; ++k1dx)
                  for (const auto& ldx : toDestination[k1dx])
                    xJ[k1dx] += x[ldx - 1];

                for (int k1dx = 0; k1dx < nOrigins; ++k1dx) {
                  fdx[k1dx] = inflowIdx[k1dx];
                  nextState[k1dx] = std::max<double>(
                    std::min<double>(stateSupport[stateIdx[k1dx]] + flowSupport[inflowIdx[k1dx]] - xI[k1dx], storageLimit),
                    0.0);
                }

                for (int k2 = 0; k2 < nDdx; ++k2) {
                  const std::vector<int>& outflowIdx = outflowIndices[k2];
                  for (int k2dx = 0; k2dx < nDestinations; ++k2dx) {
                    fdx[nOrigins + k2dx] = outflowIdx[k2dx];
                    nextState[nOrigins + k2dx] = storageLimit + std::min<double>(
                      std::max<double>(extendedStateSupport[stateIdx[nOrigins + k2dx]] - flowSupport[outflowIdx[k2dx]] + xJ[k2dx], -storageLimit),
                      storageLimit);
                  }
                  nextStateIdx = std::inner_product(nextState.begin(), nextState.end(), stateKeys.begin(), 0);

                  int kdx = std::inner_product(fdx.begin(), fdx.end(), flowKeys.begin(), 0);
                  long long p = (static_cast<long long>(i) * nAdx + j) * nQdx * nDdx * nWdx + kdx;

                  // Transit matrix columns:
                  //   [0] next_state_idx (1-based)   [1] transport_cost
                  //   [2] state_idx      (1-based)   [3] action_idx (1-based)
                  //   [4] scenario_idx   (1-based)   [5] period (1-based, = t+1)
                  transit(p, 0) = nextStateIdx + 1;
                  transit(p, 1) = objval;
                  transit(p, 2) = i + 1;
                  transit(p, 3) = j + 1;
                  transit(p, 4) = kdx + 1;
                  transit(p, 5) = t + 1;
                }
              }
            }
          }
        }
      }
    } catch (const std::exception& e) {
      #pragma omp critical
      {
        if (!anyThreadFailed.load()) {
          threadErrorMsg = std::string("Error during optimization: ") + e.what();
          anyThreadFailed.store(true);
        }
      }
    } catch (...) {
      #pragma omp critical
      {
        if (!anyThreadFailed.load()) {
          threadErrorMsg = "Unknown exception during optimization";
          anyThreadFailed.store(true);
        }
      }
    }

  } // end omp parallel

  if (anyThreadFailed.load()) {
    Rcpp::stop(threadErrorMsg);
  }

  // End timing
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;

  std::cout << "Time taken: " << elapsed.count() << " seconds" << std::endl;

  return transit;
}

//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Rcpp::List bellmanUpdateCx(
    const std::string  jsonFile,
    const int&         t,
    const std::vector<double>& stateSupport,
    const std::vector<double>& flowSupport,
    const std::vector<double>& scnpb,
    const std::vector<double>& alpha,
    const std::vector<double>& V_next,
    int numThreads = 8) {

  // Read the JSON file
  std::ifstream file(expand_path(jsonFile));
  if (!file) {
    Rcpp::stop("Unable to open file.");
  }

  nlohmann::json input;
  file >> input;

  // Extract the data from the JSON object
  const int nInventoryLevels   = input["R"][0].get<int>() + 1;
  const int nFlowLevels        = input["nQ"][0].get<int>();
  const int nSpotRates         = input["nW"][0].get<int>();
  const int nOrigins           = input["nI"][0].get<int>();
  const int nDestinations      = input["nJ"][0].get<int>();
  const int nStrategicCarriers = input["nCS"][0].get<int>();
  const int nSpotCarriers      = input["nCO"][0].get<int>();
  const int nServices          = input["nL_"][0].get<int>();
  const int nLanes             = input["nL"][0].get<int>();
  const int nActions           = input["R"][0].get<int>() + 1;
  const int storageLimit       = input["R"][0].get<int>();

  const std::vector<std::vector<int>>& fromOrigin    = input["from_i"];
  const std::vector<std::vector<int>>& toDestination = input["to_j"];
  const std::vector<std::vector<int>>& Cb    = input["Cb"];
  const std::vector<std::vector<int>>& Co    = input["Co"];
  const std::vector<std::vector<int>>& bids  = input["B"];
  const std::vector<std::vector<int>>& lanes = input["L"];
  const std::vector<std::string>& winnerKeys = input["winnerKey"];
  const std::vector<std::vector<double>>& spotRates = input["CTo"];
  const std::vector<double>& spotRateSupport = input["W"]["vals"].get<std::vector<double>>();

  const std::vector<int>& stateKeys = input["stateKeys"];
  const std::vector<int>& flowKeys  = input["flowKeys"];

  std::unordered_map<std::string, std::vector<int>>    winners;
  std::unordered_map<std::string, std::vector<double>> CTb;
  std::unordered_map<std::string, int> carrierIdx;

  winners    = importListOfVectors<int>(input["winner"]);
  CTb        = importListOfVectors<double>(input["CTb_list"]);
  carrierIdx = importList<int>(input["carrierIdx"]);

  const int n = nServices + nSpotCarriers * nLanes;
  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  const std::vector<int> limits(nOrigins + nDestinations, storageLimit);

  const std::vector<double> extendedStateSupport = utils::mirrorAndNegateVector(stateSupport);

  std::vector<std::vector<int>> stateSupportStack = stackStateIdxVectors(nInventoryLevels, nOrigins, nDestinations);
  std::vector<std::vector<int>> stateIndices = CartesianProductIntSTL(stateSupportStack);

  std::vector<int> flowIdxSingle = utils::createIndexVector(nFlowLevels);

  std::vector<std::vector<int>> inflowIdxStack;
  utils::appendIndexVectors(inflowIdxStack, flowIdxSingle, nOrigins);
  std::vector<std::vector<int>> inflowIndices = CartesianProductIntSTL(inflowIdxStack);

  std::vector<std::vector<int>> outflowIdxStack;
  utils::appendIndexVectors(outflowIdxStack, flowIdxSingle, nDestinations);
  std::vector<std::vector<int>> outflowIndices = CartesianProductIntSTL(outflowIdxStack);

  std::vector<int> spotRateIdx = utils::createIndexVector(nSpotRates);

  std::vector<std::vector<int>> spotRateIdxStack;
  utils::appendIndexVectors(spotRateIdxStack, spotRateIdx, nSpotCarriers);
  std::vector<std::vector<int>> spotRateIndices = CartesianProductIntSTL(spotRateIdxStack);

  const int nSdx = static_cast<int>(stateIndices.size());
  const int nAdx = nActions;
  const int nQdx = static_cast<int>(inflowIndices.size());
  const int nDdx = static_cast<int>(outflowIndices.size());
  const int nWdx = static_cast<int>(spotRateIndices.size());

  Highs baseHighs;
  baseHighs.setOptionValue("output_flag", false);
  baseHighs.setOptionValue("threads", 1);

  std::vector<int> baseColIdx = createTransportVars(baseHighs, n, winnerKeys, winners, bids, lanes, CTb, spotRates[t], nSpotCarriers, nLanes);
  try {
    addCapacityConstraints(baseHighs, baseColIdx, winnerKeys, winners, bids, carrierIdx, nSpotCarriers, nLanes, Cb[t], Co[t]);
    addStorageLimitConstraints(baseHighs, baseColIdx, winnerKeys, winners, bids, lanes, nSpotCarriers, nLanes, nWarehouses, limits);
    addVolumeConstraint(baseHighs, baseColIdx, n, 0); // placeholder; overridden per action in loop
  } catch (const std::exception& e) {
    Rcpp::stop("Error building base model: %s", e.what());
  }

  const HighsLp baseLp = baseHighs.getLp();

  std::vector<double> baseRhs(nStrategicCarriers + nSpotCarriers + nDestinations + nOrigins + 1, 0.0);
  for (int idx = 0; idx < nStrategicCarriers; ++idx)
    baseRhs[idx] = Cb[t][idx];
  for (int idx = 0; idx < nSpotCarriers; ++idx)
    baseRhs[nStrategicCarriers + idx] = Co[t][idx];

  // Per-(i,j) accumulators
  std::vector<double> sumW   (nSdx * nAdx, 0.0);
  std::vector<double> sumCost(nSdx * nAdx, 0.0);
  std::vector<double> sumV   (nSdx * nAdx, 0.0);
  std::vector<int>    hasFeas(nSdx * nAdx, 0);

  std::atomic<bool> anyThreadFailed{false};
  std::string       threadErrorMsg;

  #pragma omp parallel num_threads(numThreads)
  {
    Highs threadHighs;
    threadHighs.setOptionValue("output_flag", false);
    threadHighs.setOptionValue("threads", 1);
    threadHighs.passModel(baseLp);
    const std::vector<int>& colIdx = baseColIdx;

    std::vector<double> rhs           = baseRhs;
    std::vector<int>    fdx           (nOrigins + nDestinations + nSpotCarriers, 0);
    std::vector<double> nextState     (nOrigins + nDestinations, 0.0);
    std::vector<double> spotRatesTmp  (nSpotCarriers * nLanes, 0.0);
    std::vector<double> x             (n, 0.0);
    std::vector<double> xI            (nOrigins, 0.0);
    std::vector<double> xJ            (nDestinations, 0.0);
    int                 nextStateIdx = 0;

    const int volumeRow = nStrategicCarriers + nSpotCarriers + nOrigins + nDestinations;

    try {
      #pragma omp for schedule(dynamic)
      for (int i = 0; i < nSdx; ++i) {
        const auto& stateIdx = stateIndices[i];

        for (int idx = 0; idx < nDestinations; ++idx) {
          int p = nStrategicCarriers + nSpotCarriers + nOrigins + idx;
          rhs[p] = storageLimit - extendedStateSupport[stateIdx[nOrigins + idx]];
          threadHighs.changeRowBounds(p, -kHighsInf, rhs[p]);
        }

        for (int j = 0; j < nAdx; ++j) {
          threadHighs.changeRowBounds(volumeRow, static_cast<double>(j), static_cast<double>(j));

          for (int k1 = 0; k1 < nQdx; ++k1) {
            const auto& inflowIdx = inflowIndices[k1];

            for (int k1dx = 0; k1dx < nOrigins; ++k1dx) {
              int p = nStrategicCarriers + nSpotCarriers + k1dx;
              rhs[p] = stateSupport[stateIdx[k1dx]] + flowSupport[inflowIdx[k1dx]];
              threadHighs.changeRowBounds(p, -kHighsInf, rhs[p]);
            }

            std::fill(fdx.begin(), fdx.end(), 0);

            for (int k3 = 0; k3 < nWdx; ++k3) {
              const auto& spotRateIdx = spotRateIndices[k3];

              std::fill(spotRatesTmp.begin(), spotRatesTmp.end(), 0.0);
              for (int k3dx = 0; k3dx < nSpotCarriers; ++k3dx) {
                fdx[nOrigins + nDestinations + k3dx] = spotRateIdx[k3dx];
                for (int ldx = 0; ldx < nLanes; ++ldx)
                  spotRatesTmp[k3dx * nLanes + ldx] = spotRateSupport[spotRateIdx[k3dx]];
              }
              updateSpotRates(threadHighs, colIdx, spotRatesTmp, nServices, nSpotCarriers, nLanes);

              threadHighs.run();

              if (threadHighs.getModelStatus() == HighsModelStatus::kOptimal) {
                double objval = threadHighs.getObjectiveValue();
                const auto& sol = threadHighs.getSolution();

                for (int k1dx = 0; k1dx < n; ++k1dx)
                  x[k1dx] = sol.col_value[colIdx[k1dx]];

                std::fill(xI.begin(), xI.end(), 0.0);
                for (int k1dx = 0; k1dx < nOrigins; ++k1dx)
                  for (const auto& ldx : fromOrigin[k1dx])
                    xI[k1dx] += x[ldx - 1];

                std::fill(xJ.begin(), xJ.end(), 0.0);
                for (int k1dx = 0; k1dx < nDestinations; ++k1dx)
                  for (const auto& ldx : toDestination[k1dx])
                    xJ[k1dx] += x[ldx - 1];

                for (int k1dx = 0; k1dx < nOrigins; ++k1dx) {
                  fdx[k1dx] = inflowIdx[k1dx];
                  nextState[k1dx] = std::max<double>(
                    std::min<double>(stateSupport[stateIdx[k1dx]] + flowSupport[inflowIdx[k1dx]] - xI[k1dx], storageLimit),
                    0.0);
                }

                for (int k2 = 0; k2 < nDdx; ++k2) {
                  const std::vector<int>& outflowIdx = outflowIndices[k2];
                  for (int k2dx = 0; k2dx < nDestinations; ++k2dx) {
                    fdx[nOrigins + k2dx] = outflowIdx[k2dx];
                    nextState[nOrigins + k2dx] = storageLimit + std::min<double>(
                      std::max<double>(extendedStateSupport[stateIdx[nOrigins + k2dx]] - flowSupport[outflowIdx[k2dx]] + xJ[k2dx], -storageLimit),
                      storageLimit);
                  }
                  nextStateIdx = std::inner_product(nextState.begin(), nextState.end(), stateKeys.begin(), 0);

                  int kdx = std::inner_product(fdx.begin(), fdx.end(), flowKeys.begin(), 0);
                  const int ij = i * nAdx + j;
                  double pb = scnpb[kdx];
                  sumW   [ij] += pb;
                  sumCost[ij] += pb * objval;
                  sumV   [ij] += pb * V_next[nextStateIdx];
                  hasFeas[ij]  = 1;
                }
              }
            }
          }
        }
      }
    } catch (const std::exception& e) {
      #pragma omp critical
      {
        if (!anyThreadFailed.load()) {
          threadErrorMsg = std::string("Error during optimization: ") + e.what();
          anyThreadFailed.store(true);
        }
      }
    } catch (...) {
      #pragma omp critical
      {
        if (!anyThreadFailed.load()) {
          threadErrorMsg = "Unknown exception during optimization";
          anyThreadFailed.store(true);
        }
      }
    }

  } // end omp parallel

  if (anyThreadFailed.load()) {
    Rcpp::stop(threadErrorMsg);
  }

  // Post-OMP Bellman update (serial)
  Rcpp::NumericVector Q_t      (nSdx * nAdx, NA_REAL);
  Rcpp::NumericVector V_t      (nSdx,        NA_REAL);
  Rcpp::NumericVector pi_star_t(nSdx * nAdx, 0.0);
  Rcpp::NumericVector pi_rand_t(nSdx * nAdx, 0.0);

  for (int i = 0; i < nSdx; ++i) {
    const auto& stateIdx = stateIndices[i];

    double holdCost = 0.0;
    for (int k = 0; k < nOrigins; ++k)
      holdCost += stateSupport[stateIdx[k]] * alpha[k];
    for (int k = 0; k < nDestinations; ++k) {
      double sj = extendedStateSupport[stateIdx[nOrigins + k]];
      holdCost += std::max(sj, 0.0) * alpha[nOrigins + k];
      holdCost -= std::min(sj, 0.0) * alpha[nOrigins + nDestinations + k];
    }

    double Qmax = -std::numeric_limits<double>::infinity();
    double Qmin =  std::numeric_limits<double>::infinity();
    for (int j = 0; j < nAdx; ++j) {
      int ij = i * nAdx + j;
      if (!hasFeas[ij]) continue;
      double q = -holdCost - sumCost[ij] + sumV[ij];
      Q_t[ij] = q;
      Qmax = std::max(Qmax, q);
      Qmin = std::min(Qmin, q);
    }

    if (std::isfinite(Qmax)) V_t[i] = Qmax;

    double randSum = 0.0;
    for (int j = 0; j < nAdx; ++j) {
      int ij = i * nAdx + j;
      if (!Rcpp::NumericVector::is_na(Q_t[ij]))
        randSum += Q_t[ij] - Qmin + 1.0;
    }
    for (int j = 0; j < nAdx; ++j) {
      int ij = i * nAdx + j;
      if (Rcpp::NumericVector::is_na(Q_t[ij])) continue;
      pi_star_t[ij] = (Q_t[ij] == Qmax) ? 1.0 : 0.0;
      pi_rand_t[ij] = (randSum > 0.0) ? (Q_t[ij] - Qmin + 1.0) / randSum : 0.0;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("V_t")       = V_t,
    Rcpp::Named("Q_t")       = Q_t,
    Rcpp::Named("pi_star_t") = pi_star_t,
    Rcpp::Named("pi_rand_t") = pi_rand_t
  );
}

// Create a HiGHS model (replaces createGRBmodel)
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
SEXP createHIGHSmodel() {
  Highs* highs = new Highs();
  highs->setOptionValue("output_flag", false);
  return Rcpp::XPtr<Highs>(highs, true);
}

template <typename T>
std::unordered_map<std::string, std::vector<T>> rListToMap(Rcpp::List rlist) {
  std::unordered_map<std::string, std::vector<T>> result;

  Rcpp::CharacterVector names = rlist.names();

  for (int i = 0; i < rlist.size(); i++) {
    std::string name = Rcpp::as<std::string>(names[i]);
    Rcpp::Vector< Rcpp::traits::r_sexptype_traits<T>::rtype > vec = rlist[i];
    result[name] = Rcpp::as<std::vector<T>>(vec);
  }

  return result;
}

// Overload for scalars
template <typename T>
std::unordered_map<std::string, T> rListToMapScalar(Rcpp::List rlist) {
  std::unordered_map<std::string, T> result;

  Rcpp::CharacterVector names = rlist.names();

  for (int i = 0; i < rlist.size(); i++) {
    std::string name = Rcpp::as<std::string>(names[i]);
    T value = Rcpp::as<T>(rlist[i]);
    result[name] = value;
  }

  return result;
}

// Function to expose createTransportVars to R
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
SEXP createTransportVarsCx(
    SEXP model_ptr,
    const int& n,
    const std::vector<std::string>& winnerKeys,
    const Rcpp::List& winners,
    const std::vector<std::vector<int>>& bids,
    const std::vector<std::vector<int>>& lanes,
    const Rcpp::List& contractRates,
    const std::vector<double>& spotRates,
    const int& nSpotCarriers,
    const int& nLanes) {

  Rcpp::XPtr<Highs> highs(model_ptr);
  if (!highs) Rcpp::stop("model_ptr is NULL or expired");

  // Check if model already has variables of size n
  int numVars = static_cast<int>(highs->getLp().col_cost_.size());
  if (numVars == n) {
    Rcpp::Rcout << "Model already has an objective vector of size " << n << std::endl;
    // Return existing column indices (0..n-1)
    std::vector<int>* colIdx = new std::vector<int>(n);
    std::iota(colIdx->begin(), colIdx->end(), 0);
    return Rcpp::XPtr<std::vector<int>>(colIdx, true);
  }

  std::unordered_map<std::string, std::vector<int>>    winnersCx    = rListToMap<int>(winners);
  std::unordered_map<std::string, std::vector<double>> contractRatesCx = rListToMap<double>(contractRates);

  std::vector<int>* colIdx = new std::vector<int>(
      createTransportVars(*highs, n, winnerKeys, winnersCx, bids, lanes, contractRatesCx, spotRates, nSpotCarriers, nLanes));

  return Rcpp::XPtr<std::vector<int>>(colIdx, true);
}

// Function to add capacity constraints
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void addCapacityConstraintsCx(
    SEXP model_ptr,
    SEXP transport_ptr,
    const std::vector<std::string>& winnerKeys,
    const Rcpp::List& winners,
    const std::vector<std::vector<int>>& bids,
    const Rcpp::List& carrierIdx,
    const std::vector<int>& strategicCapacities,
    const std::vector<int>& spotCapacities,
    const int& nSpotSources,
    const int& nLanes) {

  Rcpp::XPtr<Highs>              highs(model_ptr);
  Rcpp::XPtr<std::vector<int>>   colIdx(transport_ptr);
  if (!highs)  Rcpp::stop("model_ptr is NULL or expired");
  if (!colIdx) Rcpp::stop("transport_ptr is NULL or expired");

  std::unordered_map<std::string, std::vector<int>> winnersCx    = rListToMap<int>(winners);
  std::unordered_map<std::string, int>              carrierIdxCx = rListToMapScalar<int>(carrierIdx);

  addCapacityConstraints(
      *highs, *colIdx,
      winnerKeys, winnersCx, bids, carrierIdxCx,
      nSpotSources, nLanes,
      strategicCapacities, spotCapacities);
}

// Function to add storage limit constraints
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void addStorageLimitConstraintsCx(
    SEXP model_ptr,
    SEXP transport_ptr,
    const std::vector<std::string>& winnerKeys,
    const Rcpp::List& winners,
    const std::vector<std::vector<int>>& bids,
    const std::vector<std::vector<int>>& lanes,
    const int& nSpotCarriers,
    const int& nLanes,
    const std::vector<int>& nWarehouses,
    const std::vector<int>& limits) {

  Rcpp::XPtr<Highs>            highs(model_ptr);
  Rcpp::XPtr<std::vector<int>> colIdx(transport_ptr);
  if (!highs)  Rcpp::stop("model_ptr is NULL or expired");
  if (!colIdx) Rcpp::stop("transport_ptr is NULL or expired");

  std::unordered_map<std::string, std::vector<int>> winnersCx = rListToMap<int>(winners);

  addStorageLimitConstraints(
      *highs, *colIdx,
      winnerKeys, winnersCx, bids, lanes,
      nSpotCarriers, nLanes,
      nWarehouses, limits);
}

// Function to add volume constraint
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void addVolumeConstraintCx(SEXP model_ptr, SEXP transport_ptr, const int& n, const double& At) {
  Rcpp::XPtr<Highs>            highs(model_ptr);
  Rcpp::XPtr<std::vector<int>> colIdx(transport_ptr);
  if (!highs)  Rcpp::stop("model_ptr is NULL or expired");
  if (!colIdx) Rcpp::stop("transport_ptr is NULL or expired");

  addVolumeConstraint(*highs, *colIdx, n, At);
}

// Function to expose printing of objective vector to R
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void printObjectiveVectorCx(SEXP model_ptr) {
  Rcpp::XPtr<Highs> highs(model_ptr);
  if (!highs) Rcpp::stop("model_ptr is NULL or expired");
  printObjectiveVector(*highs);
}

// Function to expose printing of constraints to R
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void printConstraintsCx(SEXP model_ptr) {
  Rcpp::XPtr<Highs> highs(model_ptr);
  if (!highs) Rcpp::stop("model_ptr is NULL or expired");
  printConstraints(*highs);
}

// Optimize the model and return objective and decision in R list
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Rcpp::List optimizeModelFromJSON(
    std::string jsonFile,
    const int& t,
    const std::vector<double>& spotRates,
    const std::vector<int>& storage_limits,
    const int& volume) {

  // Read the JSON file
  std::ifstream file(expand_path(jsonFile));
  if (!file) Rcpp::stop("Cannot open JSON file: %s", jsonFile.c_str());
  nlohmann::json input;
  file >> input;

  // Extract the data
  const int nOrigins        = input["nI"][0];
  const int nDestinations   = input["nJ"][0];
  const int nStrategicLanes = input["nL_"][0];
  const int nSpotSources    = input["nCO"][0];
  const int nLanes          = input["nL"][0];
  const std::vector<std::vector<int>> bids                = input["B"];
  const std::vector<std::vector<int>> lanes               = input["L"];
  const std::vector<std::vector<int>> strategicCapacities = input["Cb"];
  const std::vector<std::vector<int>> spotCapacities      = input["Co"];

  std::unordered_map<std::string, std::vector<int>>    winners;
  std::unordered_map<std::string, std::vector<double>> strategicRates;
  std::unordered_map<std::string, int> carrierIdx;

  const std::vector<std::string> winnerKeys = input["winnerKey"];

  winners        = importListOfVectors<int>(input["winner"]);
  strategicRates = importListOfVectors<double>(input["CTb_list"]);
  carrierIdx     = importList<int>(input["carrierIdx"]);

  const int n = nStrategicLanes + nSpotSources * nLanes;
  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  const std::vector<int>& limits = storage_limits;

  // Preallocate the result list
  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["objval"] = 0.0,
    Rcpp::_["x"] = Rcpp::NumericVector(n),
    Rcpp::_["status"] = "INITIAL"
  );

  try {
    Highs highs;
    highs.setOptionValue("output_flag", false);

    std::vector<int> colIdx = createTransportVars(highs, n, winnerKeys, winners, bids, lanes, strategicRates, spotRates, nSpotSources, nLanes);

    addCapacityConstraints(highs, colIdx, winnerKeys, winners, bids, carrierIdx, nSpotSources, nLanes, strategicCapacities[t], spotCapacities[t]);
    addStorageLimitConstraints(highs, colIdx, winnerKeys, winners, bids, lanes, nSpotSources, nLanes, nWarehouses, storage_limits);
    addVolumeConstraint(highs, colIdx, n, volume);

    highs.run();

    if (highs.getModelStatus() == HighsModelStatus::kOptimal) {
      result["status"] = std::string("OPTIMAL");
      result["objval"] = highs.getObjectiveValue();

      const auto& sol = highs.getSolution();
      Rcpp::NumericVector x(n);
      for (int i = 0; i < n; ++i)
        x[i] = sol.col_value[colIdx[i]];
      result["x"] = x;

    } else if (highs.getModelStatus() == HighsModelStatus::kInfeasible) {
      result["status"] = std::string("INFEASIBLE");
      Rcpp::Rcout << "Model is infeasible" << std::endl;
    } else if (highs.getModelStatus() == HighsModelStatus::kUnbounded) {
      result["status"] = std::string("UNBOUNDED");
      Rcpp::Rcout << "Model is unbounded" << std::endl;
    } else {
      result["status"] = std::string("OTHER");
      Rcpp::Rcout << "Optimization ended with non-optimal status" << std::endl;
    }

    return result;

  } catch (const std::exception& e) {
    Rcpp::Rcout << "Error: " << e.what() << std::endl;
  } catch (...) {
    Rcpp::Rcout << "Exception during optimization" << std::endl;
  }

  return result;
}

// Function to update the state index
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
std::vector<int> updateStateIdx(
    const std::vector<int>& stateIdx,
    const std::vector<int>& inflowIdx,
    const std::vector<std::vector<int>>& outflowIndices,
    const std::vector<double>& stateSupport,
    const std::vector<double>& extendedStateSupport,
    const std::vector<double>& flowSupport,
    const std::vector<double>& xI,
    const std::vector<double>& xJ,
    const double& storageLimit,
    const std::vector<int>& stateKeys,
    const int& nOrigins,
    const int& nDestinations) {

  std::vector<int> nextStateIdx(outflowIndices.size(), 0);
  std::vector<double> nextState(nOrigins + nDestinations, 0.0);

  for (int i = 0; i < nOrigins; ++i) {
    nextState[i] = std::max(std::min(stateSupport[stateIdx[i]] + flowSupport[inflowIdx[i]] - xI[i], storageLimit), 0.0);
  }

  for (size_t k = 0; k < outflowIndices.size(); ++k) {
    const std::vector<int>& outflowIdx = outflowIndices[k];
    for (int j = 0; j < nDestinations; ++j) {
      nextState[nOrigins + j] = storageLimit + std::min(std::max(extendedStateSupport[stateIdx[nOrigins + j]] - flowSupport[outflowIdx[j]] + xJ[j], -storageLimit), storageLimit);
    }
    nextStateIdx[k] = std::inner_product(nextState.begin(), nextState.end(), stateKeys.begin(), 0);
  }

  return nextStateIdx;
}

// Redirect C-level stdout to /dev/null; returns the saved fd for restoration.
// Used to suppress the HiGHS startup banner and deprecated-API warning printed
// by the 'highs' R package before option-setting takes effect.
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
int begin_suppress_stdout() {
  fflush(stdout);
  int saved = dup(1);
  if (saved == -1) return -1;
  int devnull = open("/dev/null", O_WRONLY);
  if (devnull == -1) { close(saved); return -1; }
  if (dup2(devnull, 1) == -1) { close(devnull); close(saved); return -1; }
  close(devnull);
  return saved;
}

//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void end_suppress_stdout(int saved_fd) {
  if (saved_fd == -1) return;
  fflush(stdout);
  dup2(saved_fd, 1);
  close(saved_fd);
}

#ifdef STANDALONE_BUILD
int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return(EXIT_FAILURE);
  }

  std::ifstream file(argv[1]);
  nlohmann::json input;
  file >> input;

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
  std::unordered_map<std::string, int> carrierIdx;

  const int n = nStrategicSources + nSpotCarriers * nL;
  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  const std::vector<std::string> winnerKeys = input["winnerKey"];

  winners = importListOfVectors<int>(input["winner"]);
  CTb = importListOfVectors<double>(input["CTb_list"]);
  carrierIdx = importList<int>(input["carrierIdx"]);

  const std::vector<int> limits = utils::generateRandomIntegers(nOrigins + nDestinations, 20, 40);

  try {
    Highs highs;
    std::vector<int> colIdx = createTransportVars(highs, n, winnerKeys, winners, bids, lanes, CTb, spotRates[0], nSpotCarriers, nL);

    addCapacityConstraints(highs, colIdx, winnerKeys, winners, bids, carrierIdx, nSpotCarriers, nL, Cb[0], Co[0]);
    addStorageLimitConstraints(highs, colIdx, winnerKeys, winners, bids, lanes, nSpotCarriers, nL, nWarehouses, limits);
    addVolumeConstraint(highs, colIdx, n, 10);

  } catch (const std::exception& e) {
    cout << "Error: " << e.what() << endl;
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}
#endif // STANDALONE_BUILD
