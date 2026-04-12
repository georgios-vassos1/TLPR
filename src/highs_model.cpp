#include <Highs.h>
#include <atomic>
#include <fstream>
#include <Rcpp.h>
#include <omp.h>
#include <unistd.h>
#include <fcntl.h>

#ifdef STANDALONE_BUILD
  #include "reader.hpp"
  #include "utils.hpp"
  #include "cartesian.hpp"
  #include "hilbert.hpp"
#else
  #include "../inst/include/reader.hpp"
  #include "../inst/include/utils.hpp"
  #include "../inst/include/cartesian.hpp"
  #include "../inst/include/hilbert.hpp"
#endif

// Explicit using-declarations instead of 'using namespace std' to avoid
// pulling in the entire standard namespace into library translation units.
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;

// State traversal order for bellmanUpdateImpl.
// Hilbert order improves warm-start basis reuse by visiting geometrically
// adjacent states consecutively.  Lexicographic is the default (no reordering).
enum class StateOrder { Lexicographic, Hilbert };

// Sort stateIndices in-place by dual Hilbert key.
// Origin sub-state  (dims 0..nI-1) encoded with bitsFor(R).
// Destination sub-state (dims nI..nI+nJ-1) encoded after S+/S- split.
// Ties in origin key are broken by destination key.
static void sortStatesByHilbert(std::vector<std::vector<int>>& stateIndices, int nOrigins,
                                int nDestinations, int R) {
  const int b = hilbert::bitsFor(R);
  std::sort(stateIndices.begin(), stateIndices.end(),
            [&](const std::vector<int>& a, const std::vector<int>& bv) {
              // Origin Hilbert key
              const uint64_t ha = hilbert::encode(a.data(), nOrigins, b);
              const uint64_t hb = hilbert::encode(bv.data(), nOrigins, b);
              if (ha != hb) {
                return ha < hb;
              }
              // Destination Hilbert key via S+/S- split
              std::vector<int> da, db;
              da.reserve(2 * nDestinations);
              db.reserve(2 * nDestinations);
              for (int j = 0; j < nDestinations; ++j) {
                auto [spa, sma] = utils::splitSignedIdx(a[nOrigins + j], R);
                auto [spb, smb] = utils::splitSignedIdx(bv[nOrigins + j], R);
                da.push_back(spa);
                da.push_back(sma);
                db.push_back(spb);
                db.push_back(smb);
              }
              return hilbert::encode(da, b) < hilbert::encode(db, b);
            });
}

// Function to stack state index vectors
std::vector<std::vector<int>> stackStateIdxVectors(int nInventoryLevels, int nOrigins,
                                                   int nDestinations) {

  std::vector<std::vector<int>> stateIdx;
  stateIdx.reserve(nOrigins + nDestinations);

  const auto idx = utils::createIndexVector(nInventoryLevels);
  const auto jdx = utils::createIndexVector(2 * nInventoryLevels - 1);

  utils::appendIndexVectors(stateIdx, idx, nOrigins);
  utils::appendIndexVectors(stateIdx, jdx, nDestinations);

  return stateIdx;
}

// Function to print the objective vector in a HiGHS model
void printObjectiveVector(const Highs& highs) {
  const HighsLp& lp = highs.getLp();
  int numVars = static_cast<int>(lp.col_cost_.size());
  Rcpp::Rcout << "Number of decision variables: " << numVars << '\n';
  Rcpp::Rcout << "Objective Vector: " << '\n';
  for (int i = 0; i < numVars; ++i) {
    Rcpp::Rcout << "x[" << i << "]: " << lp.col_cost_[i] << '\n';
  }
}

// Function to print the constraints in a HiGHS model
void printConstraints(const Highs& highs) {
  const HighsLp& lp = highs.getLp();
  int numConstraints = static_cast<int>(lp.row_lower_.size());
  Rcpp::Rcout << "Constraints: " << '\n';
  for (int i = 0; i < numConstraints; ++i) {
    Rcpp::Rcout << "row[" << i << "]: lower=" << lp.row_lower_[i] << " upper=" << lp.row_upper_[i]
                << '\n';
  }
}

// Function to print the optimal transport volumes
void printOptimalTransportVolumes(const Highs& highs, const std::vector<int>& colIdx, int n) {
  const auto& sol = highs.getSolution();
  Rcpp::Rcout << "Transport volumes:" << '\n';
  for (int i = 0; i < n; ++i) {
    Rcpp::Rcout << "x[" << colIdx[i] << "] = " << sol.col_value[colIdx[i]] << '\n';
  }
}

// Function to update spot rates in the objective vector
void updateSpotRates(Highs& highs, const std::vector<int>& colIdx,
                     const std::vector<double>& spotRates, int nStrategicSources, int nSpotCarriers,
                     int nLanes) {

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
    if (home != nullptr) {
      return std::string(home) + path.substr(1);
    }
  }
  return path;
}

// ─────────────────────────────────────────────────────────────────────────
//  ProblemData: holds all fields parsed from an instance JSON so that the
//  data can be loaded once and reused across multiple C++ solver calls,
//  eliminating the per-call file I/O and JSON parsing overhead.
// ─────────────────────────────────────────────────────────────────────────
struct ProblemData {
  // Scalar dimensions
  int R{};   // storage limit (= max inventory = nActions - 1)
  int nQ{};  // number of flow levels
  int nW{};  // number of spot-rate support points
  int nI{};  // number of origins
  int nJ{};  // number of destinations
  int nCS{}; // number of strategic carriers
  int nCO{}; // number of spot carriers
  int nL_{}; // number of strategic services
  int nL{};  // total number of lanes

  // Network topology
  std::vector<std::vector<int>> fromOrigin{};    // from_i
  std::vector<std::vector<int>> toDestination{}; // to_j
  std::vector<std::vector<int>> Cb{};            // strategic capacities [t][carrier]
  std::vector<std::vector<int>> Co{};            // spot capacities      [t][carrier]
  std::vector<std::vector<int>> bids{};          // B
  std::vector<std::vector<int>> lanes{};         // L
  std::vector<std::string> winnerKeys{};         // winnerKey
  std::vector<std::vector<double>> spotRates{};  // CTo [t][carrier*lane]
  std::vector<double> spotRateSupport{};         // W.vals

  // Cipher keys for state / scenario indexing
  std::vector<int> stateKeys{};
  std::vector<int> flowKeys{};

  // Carrier data (from R named-list fields)
  std::unordered_map<std::string, std::vector<int>> winners{};
  std::unordered_map<std::string, std::vector<double>> CTb{};
  std::unordered_map<std::string, int> carrierIdx{};
};

// Parse a JSON instance file into a ProblemData value.
// Shared by loadProblemDataCx (XPtr path) and the legacy JSON wrappers.
[[nodiscard]] static ProblemData parseProblemData(const std::string& jsonFile) {
  std::ifstream file(expand_path(jsonFile));
  if (!file) {
    Rcpp::stop("Cannot open JSON file: %s", jsonFile.c_str());
  }

  nlohmann::json input;
  file >> input;

  ProblemData pd;

  pd.R = input["R"][0].get<int>();
  pd.nQ = input["nQ"][0].get<int>();
  pd.nW = input["nW"][0].get<int>();
  pd.nI = input["nI"][0].get<int>();
  pd.nJ = input["nJ"][0].get<int>();
  pd.nCS = input["nCS"][0].get<int>();
  pd.nCO = input["nCO"][0].get<int>();
  pd.nL_ = input["nL_"][0].get<int>();
  pd.nL = input["nL"][0].get<int>();

  pd.fromOrigin = input["from_i"].get<std::vector<std::vector<int>>>();
  pd.toDestination = input["to_j"].get<std::vector<std::vector<int>>>();
  pd.Cb = input["Cb"].get<std::vector<std::vector<int>>>();
  pd.Co = input["Co"].get<std::vector<std::vector<int>>>();
  pd.bids = input["B"].get<std::vector<std::vector<int>>>();
  pd.lanes = input["L"].get<std::vector<std::vector<int>>>();
  pd.winnerKeys = input["winnerKey"].get<std::vector<std::string>>();
  pd.spotRates = input["CTo"].get<std::vector<std::vector<double>>>();
  pd.spotRateSupport = input["W"]["vals"].get<std::vector<double>>();

  pd.stateKeys = input["stateKeys"].get<std::vector<int>>();
  pd.flowKeys = input["flowKeys"].get<std::vector<int>>();

  pd.winners = importListOfVectors<int>(input["winner"]);
  pd.CTb = importListOfVectors<double>(input["CTb_list"]);
  pd.carrierIdx = importList<int>(input["carrierIdx"]);

  return pd;
}

// Load problem data from a JSON file and return an external pointer.
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
SEXP loadProblemDataCx(const std::string& jsonFile) {
  return Rcpp::XPtr<ProblemData>(new ProblemData(parseProblemData(jsonFile)), true);
}

// Function to create decision variables for transportation
[[nodiscard]] std::vector<int>
createTransportVars(Highs& highs, int n, const std::vector<std::string>& winnerKeys,
                    const std::unordered_map<std::string, std::vector<int>>& winners,
                    const std::vector<std::vector<int>>& bids,
                    [[maybe_unused]] const std::vector<std::vector<int>>& lanes,
                    const std::unordered_map<std::string, std::vector<double>>& contractRates,
                    const std::vector<double>& spotRates, int nSpotCarriers, int nLanes) {

  std::vector<int> colIdx;
  colIdx.reserve(n);

  // Strategic carrier variables
  for (const auto& winnerKey : winnerKeys) {
    int k = 0;
    for (int bidIndex : winners.at(winnerKey)) {
      const auto& bid = bids[bidIndex - 1];
      for ([[maybe_unused]] int laneIndex : bid) {
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

void addCapacityConstraints(Highs& highs, const std::vector<int>& colIdx,
                            const std::vector<std::string>& winnerKeys,
                            const std::unordered_map<std::string, std::vector<int>>& winners,
                            const std::vector<std::vector<int>>& bids,
                            const std::unordered_map<std::string, int>& carrierIdx,
                            int nSpotCarriers, int nLanes, const std::vector<int>& strategicCaps,
                            const std::vector<int>& spotCaps) {

  int k = 0;

  // Strategic carrier capacity constraints
  for (const auto& winnerKey : winnerKeys) {
    const int capacity = strategicCaps.at(carrierIdx.at(winnerKey) - 1);
    std::vector<int> ids;
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
    std::vector<int> ids;
    std::vector<double> vals;
    for (int laneIndex = 0; laneIndex < nLanes; ++laneIndex) {
      ids.push_back(colIdx[k++]);
      vals.push_back(1.0);
    }
    highs.addRow(-kHighsInf, capacity, static_cast<int>(ids.size()), ids.data(), vals.data());
  }
}

void addStorageLimitConstraints(Highs& highs, const std::vector<int>& colIdx,
                                const std::vector<std::string>& winnerKeys,
                                const std::unordered_map<std::string, std::vector<int>>& winners,
                                const std::vector<std::vector<int>>& bids,
                                const std::vector<std::vector<int>>& lanes, int nSpotCarriers,
                                int nLanes, const std::vector<int>& nWarehouses,
                                const std::vector<int>& limits) {

  int position = 0;
  for (size_t m = 0; m < nWarehouses.size(); ++m) {
    for (int i = 0; i < nWarehouses[m]; ++i) {
      std::vector<int> ids;
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

      highs.addRow(-kHighsInf, limits[position + i], static_cast<int>(ids.size()), ids.data(),
                   vals.data());
    }
    position += nWarehouses[m];
  }
}

// Volume constraint
void addVolumeConstraint(Highs& highs, const std::vector<int>& colIdx, int n, double At) {

  std::vector<double> vals(n, 1.0);
  highs.addRow(At, At, n, colIdx.data(), vals.data());
}

// Canonical implementation — takes a resolved ProblemData reference.
// Called by both the JSON wrapper and the XPtr wrapper below.
[[nodiscard]] static Eigen::MatrixXd computeEnvironmentImpl(const ProblemData& pd, int t,
                                                            const std::vector<double>& stateSupport,
                                                            const std::vector<double>& flowSupport,
                                                            int numThreads) {

  const int nInventoryLevels = pd.R + 1;
  const int nFlowLevels = pd.nQ;
  const int nSpotRates = pd.nW;
  const int nOrigins = pd.nI;
  const int nDestinations = pd.nJ;
  const int nStrategicCarriers = pd.nCS;
  const int nSpotCarriers = pd.nCO;
  const int nServices = pd.nL_;
  const int nLanes = pd.nL;
  const int nActions = pd.R + 1;
  const int storageLimit = pd.R;

  const auto& fromOrigin = pd.fromOrigin;
  const auto& toDestination = pd.toDestination;
  const auto& Cb = pd.Cb;
  const auto& Co = pd.Co;
  const auto& bids = pd.bids;
  const auto& lanes = pd.lanes;
  const auto& winnerKeys = pd.winnerKeys;
  const auto& spotRates = pd.spotRates;
  const auto& spotRateSupport = pd.spotRateSupport;
  const auto& stateKeys = pd.stateKeys;
  const auto& flowKeys = pd.flowKeys;
  const auto& winners = pd.winners;
  const auto& CTb = pd.CTb;
  const auto& carrierIdx = pd.carrierIdx;

  const int n = nServices + nSpotCarriers * nLanes;

  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  const std::vector<int> limits(nOrigins + nDestinations, storageLimit);

  // Generate extended state support data (positive and negative states)
  const std::vector<double> extendedStateSupport = utils::mirrorAndNegateVector(stateSupport);

  // Generate state support data (only positive states)
  std::vector<std::vector<int>> stateSupportStack =
      stackStateIdxVectors(nInventoryLevels, nOrigins, nDestinations);
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

  // Initialize the Eigen matrix with n rows and 6 columns
  Eigen::MatrixXd transit(nTransit, 6);
  transit.setConstant(std::numeric_limits<double>::quiet_NaN());

  // Build the base model once, serially.
  Highs baseHighs;
  baseHighs.setOptionValue("output_flag", false);
  baseHighs.setOptionValue("threads", 1);

  std::vector<int> baseColIdx = createTransportVars(baseHighs, n, winnerKeys, winners, bids, lanes,
                                                    CTb, spotRates[t], nSpotCarriers, nLanes);
  try {
    addCapacityConstraints(baseHighs, baseColIdx, winnerKeys, winners, bids, carrierIdx,
                           nSpotCarriers, nLanes, Cb[t], Co[t]);
    addStorageLimitConstraints(baseHighs, baseColIdx, winnerKeys, winners, bids, lanes,
                               nSpotCarriers, nLanes, nWarehouses, limits);
    addVolumeConstraint(baseHighs, baseColIdx, n, 0); // placeholder; overridden per action in loop
  } catch (const std::exception& e) {
    Rcpp::stop("Error building base model: %s", e.what());
  }

  // Extract LP for OMP cloning
  const HighsLp baseLp = baseHighs.getLp();

  // Capacity RHS values are constant across all states and actions within t.
  std::vector<double> baseRhs(nStrategicCarriers + nSpotCarriers + nDestinations + nOrigins + 1,
                              0.0);
  for (int idx = 0; idx < nStrategicCarriers; ++idx) {
    baseRhs[idx] = Cb[t][idx];
  }
  for (int idx = 0; idx < nSpotCarriers; ++idx) {
    baseRhs[nStrategicCarriers + idx] = Co[t][idx];
  }

  // Error reporting across OMP threads
  std::atomic<bool> anyThreadFailed{false};
  std::string threadErrorMsg;

#pragma omp parallel num_threads(numThreads)
  {
    Highs threadHighs;
    threadHighs.setOptionValue("output_flag", false);
    threadHighs.setOptionValue("threads", 1);
    threadHighs.passModel(baseLp);
    const std::vector<int>& colIdx = baseColIdx;

    // Thread-local working storage
    std::vector<double> rhs = baseRhs;
    std::vector<int> fdx(nOrigins + nDestinations + nSpotCarriers, 0);
    std::vector<double> nextState(nOrigins + nDestinations, 0.0);
    std::vector<double> spotRatesTmp(nSpotCarriers * nLanes, 0.0);
    std::vector<double> x(n, 0.0);
    std::vector<double> xI(nOrigins, 0.0);
    std::vector<double> xJ(nDestinations, 0.0);
    int nextStateIdx = 0;
    HighsBasis lastBasis;
    bool hasBasis = false;

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
                for (int ldx = 0; ldx < nLanes; ++ldx) {
                  spotRatesTmp[k3dx * nLanes + ldx] = spotRateSupport[spotRateIdx[k3dx]];
                }
              }
              updateSpotRates(threadHighs, colIdx, spotRatesTmp, nServices, nSpotCarriers, nLanes);

              if (hasBasis) {
                threadHighs.setBasis(lastBasis);
              }
              threadHighs.run();

              if (threadHighs.getModelStatus() == HighsModelStatus::kOptimal) {
                lastBasis = threadHighs.getBasis();
                hasBasis = true;
                double objval = threadHighs.getObjectiveValue();
                const auto& sol = threadHighs.getSolution();

                for (int k1dx = 0; k1dx < n; ++k1dx) {
                  x[k1dx] = sol.col_value[colIdx[k1dx]];
                }

                std::fill(xI.begin(), xI.end(), 0.0);
                for (int k1dx = 0; k1dx < nOrigins; ++k1dx) {
                  for (const auto& ldx : fromOrigin[k1dx]) {
                    xI[k1dx] += x[ldx - 1];
                  }
                }

                std::fill(xJ.begin(), xJ.end(), 0.0);
                for (int k1dx = 0; k1dx < nDestinations; ++k1dx) {
                  for (const auto& ldx : toDestination[k1dx]) {
                    xJ[k1dx] += x[ldx - 1];
                  }
                }

                for (int k1dx = 0; k1dx < nOrigins; ++k1dx) {
                  fdx[k1dx] = inflowIdx[k1dx];
                  nextState[k1dx] =
                      std::max<double>(std::min<double>(stateSupport[stateIdx[k1dx]] +
                                                            flowSupport[inflowIdx[k1dx]] - xI[k1dx],
                                                        storageLimit),
                                       0.0);
                }

                for (int k2 = 0; k2 < nDdx; ++k2) {
                  const std::vector<int>& outflowIdx = outflowIndices[k2];
                  for (int k2dx = 0; k2dx < nDestinations; ++k2dx) {
                    fdx[nOrigins + k2dx] = outflowIdx[k2dx];
                    nextState[nOrigins + k2dx] =
                        storageLimit +
                        std::min<double>(
                            std::max<double>(extendedStateSupport[stateIdx[nOrigins + k2dx]] -
                                                 flowSupport[outflowIdx[k2dx]] + xJ[k2dx],
                                             -storageLimit),
                            storageLimit);
                  }
                  nextStateIdx =
                      std::inner_product(nextState.begin(), nextState.end(), stateKeys.begin(), 0);

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

  return transit;
}

// JSON wrapper — parses file then delegates to the canonical impl.
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd computeEnvironmentCx(const std::string& jsonFile, int t,
                                     const std::vector<double>& stateSupport,
                                     const std::vector<double>& flowSupport, int numThreads = 8) {
  return computeEnvironmentImpl(parseProblemData(jsonFile), t, stateSupport, flowSupport,
                                numThreads);
}

// XPtr wrapper — dereferences pointer then delegates to the canonical impl.
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd computeEnvironmentPtr(SEXP problem_ptr, int t,
                                      const std::vector<double>& stateSupport,
                                      const std::vector<double>& flowSupport, int numThreads = 8) {
  Rcpp::XPtr<ProblemData> pd(problem_ptr);
  if (!pd) {
    Rcpp::stop("problem_ptr is NULL or expired");
  }
  return computeEnvironmentImpl(*pd, t, stateSupport, flowSupport, numThreads);
}

// Canonical implementation — takes a resolved ProblemData reference.
// Called by both the JSON wrapper and the XPtr wrapper below.
[[nodiscard]] static Rcpp::List
bellmanUpdateImpl(const ProblemData& pd, int t, const std::vector<double>& stateSupport,
                  const std::vector<double>& flowSupport, const std::vector<double>& scnpb,
                  const std::vector<double>& alpha, const std::vector<double>& V_next,
                  int numThreads, StateOrder order = StateOrder::Lexicographic, int chunkSize = 32,
                  const std::vector<int>& stateSubset = {}) {

  const int nInventoryLevels = pd.R + 1;
  const int nFlowLevels = pd.nQ;
  const int nSpotRates = pd.nW;
  const int nOrigins = pd.nI;
  const int nDestinations = pd.nJ;
  const int nStrategicCarriers = pd.nCS;
  const int nSpotCarriers = pd.nCO;
  const int nServices = pd.nL_;
  const int nLanes = pd.nL;
  const int nActions = pd.R + 1;
  const int storageLimit = pd.R;

  const auto& fromOrigin = pd.fromOrigin;
  const auto& toDestination = pd.toDestination;
  const auto& Cb = pd.Cb;
  const auto& Co = pd.Co;
  const auto& bids = pd.bids;
  const auto& lanes = pd.lanes;
  const auto& winnerKeys = pd.winnerKeys;
  const auto& spotRates = pd.spotRates;
  const auto& spotRateSupport = pd.spotRateSupport;
  const auto& stateKeys = pd.stateKeys;
  const auto& flowKeys = pd.flowKeys;
  const auto& winners = pd.winners;
  const auto& CTb = pd.CTb;
  const auto& carrierIdx = pd.carrierIdx;

  const int n = nServices + nSpotCarriers * nLanes;
  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  const std::vector<int> limits(nOrigins + nDestinations, storageLimit);

  const std::vector<double> extendedStateSupport = utils::mirrorAndNegateVector(stateSupport);

  std::vector<std::vector<int>> stateSupportStack =
      stackStateIdxVectors(nInventoryLevels, nOrigins, nDestinations);
  std::vector<std::vector<int>> stateIndices = CartesianProductIntSTL(stateSupportStack);

  // Full grid size — used for output array allocation so that callers can
  // always index by the canonical lex state index regardless of sampling.
  const int nSdxFull = static_cast<int>(stateIndices.size());

  // Optionally restrict to a caller-supplied subset of states (0-based lex
  // indices).  Enables sampled Approximate Value Iteration: only the LP
  // subproblems for the sampled states are solved; output arrays retain full
  // nSdxFull length with NA/0 for unsampled positions.
  if (!stateSubset.empty()) {
    std::vector<std::vector<int>> filtered;
    filtered.reserve(stateSubset.size());
    for (int idx : stateSubset) {
      if (idx >= 0 && idx < nSdxFull) {
        filtered.push_back(stateIndices[idx]);
      }
    }
    stateIndices = std::move(filtered);
  }

  if (order == StateOrder::Hilbert) {
    sortStatesByHilbert(stateIndices, nOrigins, nDestinations, storageLimit);
  }

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

  // Precompute canonical lexicographic index for each state in traversal order.
  // In Lexicographic mode lexPos[i] == i; in Hilbert mode it gives the correct
  // mixed-radix position so that V_t and accumulators share the same indexing
  // as V_next[nextStateIdx].
  std::vector<int> lexPos(nSdx);
  for (int i = 0; i < nSdx; ++i) {
    lexPos[i] =
        std::inner_product(stateIndices[i].begin(), stateIndices[i].end(), stateKeys.begin(), 0);
  }

  Highs baseHighs;
  baseHighs.setOptionValue("output_flag", false);
  baseHighs.setOptionValue("threads", 1);

  std::vector<int> baseColIdx = createTransportVars(baseHighs, n, winnerKeys, winners, bids, lanes,
                                                    CTb, spotRates[t], nSpotCarriers, nLanes);
  try {
    addCapacityConstraints(baseHighs, baseColIdx, winnerKeys, winners, bids, carrierIdx,
                           nSpotCarriers, nLanes, Cb[t], Co[t]);
    addStorageLimitConstraints(baseHighs, baseColIdx, winnerKeys, winners, bids, lanes,
                               nSpotCarriers, nLanes, nWarehouses, limits);
    addVolumeConstraint(baseHighs, baseColIdx, n, 0); // placeholder; overridden per action in loop
  } catch (const std::exception& e) {
    Rcpp::stop("Error building base model: %s", e.what());
  }

  const HighsLp baseLp = baseHighs.getLp();

  std::vector<double> baseRhs(nStrategicCarriers + nSpotCarriers + nDestinations + nOrigins + 1,
                              0.0);
  for (int idx = 0; idx < nStrategicCarriers; ++idx) {
    baseRhs[idx] = Cb[t][idx];
  }
  for (int idx = 0; idx < nSpotCarriers; ++idx) {
    baseRhs[nStrategicCarriers + idx] = Co[t][idx];
  }

  // Per-(i,j) accumulators — sized to the full grid so that lexPos[i] (which
  // maps to the canonical full-grid position) is always a valid index.
  std::vector<double> sumW(nSdxFull * nAdx, 0.0);
  std::vector<double> sumCost(nSdxFull * nAdx, 0.0);
  std::vector<double> sumV(nSdxFull * nAdx, 0.0);
  std::vector<int> hasFeas(nSdxFull * nAdx, 0);

  std::atomic<bool> anyThreadFailed{false};
  std::string threadErrorMsg;

#pragma omp parallel num_threads(numThreads)
  {
    Highs threadHighs;
    threadHighs.setOptionValue("output_flag", false);
    threadHighs.setOptionValue("threads", 1);
    threadHighs.passModel(baseLp);
    const std::vector<int>& colIdx = baseColIdx;

    std::vector<double> rhs = baseRhs;
    std::vector<int> fdx(nOrigins + nDestinations + nSpotCarriers, 0);
    std::vector<double> nextState(nOrigins + nDestinations, 0.0);
    std::vector<double> spotRatesTmp(nSpotCarriers * nLanes, 0.0);
    std::vector<double> x(n, 0.0);
    std::vector<double> xI(nOrigins, 0.0);
    std::vector<double> xJ(nDestinations, 0.0);
    int nextStateIdx = 0;
    HighsBasis lastBasis;
    bool hasBasis = false;

    const int volumeRow = nStrategicCarriers + nSpotCarriers + nOrigins + nDestinations;
    const int ompChunk = (order == StateOrder::Hilbert) ? chunkSize : 1;

    try {
#pragma omp for schedule(dynamic, ompChunk)
      for (int i = 0; i < nSdx; ++i) {
        const auto& stateIdx = stateIndices[i];
        // Reset HiGHS internal basis to force a cold start for this state.
        // Without this, threadHighs retains its last optimal basis across
        // states even when hasBasis==false, causing traversal-order-dependent
        // LP tie-breaking in degenerate problems.
        threadHighs.setBasis();
        hasBasis = false;

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
                for (int ldx = 0; ldx < nLanes; ++ldx) {
                  spotRatesTmp[k3dx * nLanes + ldx] = spotRateSupport[spotRateIdx[k3dx]];
                }
              }
              updateSpotRates(threadHighs, colIdx, spotRatesTmp, nServices, nSpotCarriers, nLanes);

              if (hasBasis) {
                threadHighs.setBasis(lastBasis);
              }
              threadHighs.run();

              if (threadHighs.getModelStatus() == HighsModelStatus::kOptimal) {
                lastBasis = threadHighs.getBasis();
                hasBasis = true;
                double objval = threadHighs.getObjectiveValue();
                const auto& sol = threadHighs.getSolution();

                for (int k1dx = 0; k1dx < n; ++k1dx) {
                  x[k1dx] = sol.col_value[colIdx[k1dx]];
                }

                std::fill(xI.begin(), xI.end(), 0.0);
                for (int k1dx = 0; k1dx < nOrigins; ++k1dx) {
                  for (const auto& ldx : fromOrigin[k1dx]) {
                    xI[k1dx] += x[ldx - 1];
                  }
                }

                std::fill(xJ.begin(), xJ.end(), 0.0);
                for (int k1dx = 0; k1dx < nDestinations; ++k1dx) {
                  for (const auto& ldx : toDestination[k1dx]) {
                    xJ[k1dx] += x[ldx - 1];
                  }
                }

                for (int k1dx = 0; k1dx < nOrigins; ++k1dx) {
                  fdx[k1dx] = inflowIdx[k1dx];
                  nextState[k1dx] =
                      std::max<double>(std::min<double>(stateSupport[stateIdx[k1dx]] +
                                                            flowSupport[inflowIdx[k1dx]] - xI[k1dx],
                                                        storageLimit),
                                       0.0);
                }

                for (int k2 = 0; k2 < nDdx; ++k2) {
                  const std::vector<int>& outflowIdx = outflowIndices[k2];
                  for (int k2dx = 0; k2dx < nDestinations; ++k2dx) {
                    fdx[nOrigins + k2dx] = outflowIdx[k2dx];
                    nextState[nOrigins + k2dx] =
                        storageLimit +
                        std::min<double>(
                            std::max<double>(extendedStateSupport[stateIdx[nOrigins + k2dx]] -
                                                 flowSupport[outflowIdx[k2dx]] + xJ[k2dx],
                                             -storageLimit),
                            storageLimit);
                  }
                  nextStateIdx =
                      std::inner_product(nextState.begin(), nextState.end(), stateKeys.begin(), 0);

                  int kdx = std::inner_product(fdx.begin(), fdx.end(), flowKeys.begin(), 0);
                  const int ij = lexPos[i] * nAdx + j;
                  double pb = scnpb[kdx];
                  sumW[ij] += pb;
                  sumCost[ij] += pb * objval;
                  sumV[ij] += pb * V_next[nextStateIdx];
                  hasFeas[ij] = 1;
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
  // Output arrays are full-grid sized so R can always index by lex state index.
  // Unsampled states retain NA (V_t, Q_t) or 0 (pi_*).
  Rcpp::NumericVector Q_t(nSdxFull * nAdx, NA_REAL);
  Rcpp::NumericVector V_t(nSdxFull, NA_REAL);
  Rcpp::NumericVector pi_star_t(nSdxFull * nAdx, 0.0);
  Rcpp::NumericVector pi_rand_t(nSdxFull * nAdx, 0.0);

  for (int i = 0; i < nSdx; ++i) {
    const auto& stateIdx = stateIndices[i];
    const int lexI = lexPos[i];

    double holdCost = 0.0;
    for (int k = 0; k < nOrigins; ++k) {
      holdCost += stateSupport[stateIdx[k]] * alpha[k];
    }
    for (int k = 0; k < nDestinations; ++k) {
      double sj = extendedStateSupport[stateIdx[nOrigins + k]];
      holdCost += std::max(sj, 0.0) * alpha[nOrigins + k];
      holdCost -= std::min(sj, 0.0) * alpha[nOrigins + nDestinations + k];
    }

    double Qmax{};
    double Qmin{};
    bool hasAnyFeasible{false};
    for (int j = 0; j < nAdx; ++j) {
      int ij = lexI * nAdx + j;
      if (hasFeas[ij] == 0) {
        continue;
      }
      double q = -holdCost - sumCost[ij] + sumV[ij];
      Q_t[ij] = q;
      if (!hasAnyFeasible) { Qmax = q; Qmin = q; hasAnyFeasible = true; }
      else { Qmax = std::max(Qmax, q); Qmin = std::min(Qmin, q); }
    }

    if (hasAnyFeasible) {
      V_t[lexI] = Qmax;
    }

    double randSum = 0.0;
    for (int j = 0; j < nAdx; ++j) {
      int ij = lexI * nAdx + j;
      if (!Rcpp::NumericVector::is_na(Q_t[ij])) {
        randSum += Q_t[ij] - Qmin + 1.0;
      }
    }
    for (int j = 0; j < nAdx; ++j) {
      int ij = lexI * nAdx + j;
      if (Rcpp::NumericVector::is_na(Q_t[ij])) {
        continue;
      }
      pi_star_t[ij] = (Q_t[ij] == Qmax) ? 1.0 : 0.0;
      pi_rand_t[ij] = (randSum > 0.0) ? (Q_t[ij] - Qmin + 1.0) / randSum : 0.0;
    }
  }

  return Rcpp::List::create(Rcpp::Named("V_t") = V_t, Rcpp::Named("Q_t") = Q_t,
                            Rcpp::Named("pi_star_t") = pi_star_t,
                            Rcpp::Named("pi_rand_t") = pi_rand_t);
}

// JSON wrapper — parses file then delegates to the canonical impl.
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Rcpp::List bellmanUpdateCx(const std::string& jsonFile, int t,
                           const std::vector<double>& stateSupport,
                           const std::vector<double>& flowSupport, const std::vector<double>& scnpb,
                           const std::vector<double>& alpha, const std::vector<double>& V_next,
                           int numThreads = 8, const std::string& traversalOrder = "lexicographic",
                           int chunkSize = 32, SEXP stateSubset = R_NilValue) {
  StateOrder order =
      (traversalOrder == "hilbert") ? StateOrder::Hilbert : StateOrder::Lexicographic;
  std::vector<int> subset;
  if (stateSubset != R_NilValue && Rf_length(stateSubset) > 0) {
    const int* p = INTEGER(stateSubset);
    const int n = Rf_length(stateSubset);
    subset.resize(n);
    std::transform(p, p + n, subset.begin(), [](int x) { return x - 1; });
  }
  return bellmanUpdateImpl(parseProblemData(jsonFile), t, stateSupport, flowSupport, scnpb, alpha,
                           V_next, numThreads, order, chunkSize, subset);
}

// XPtr wrapper — dereferences pointer then delegates to the canonical impl.
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Rcpp::List bellmanUpdatePtr(SEXP problem_ptr, int t, const std::vector<double>& stateSupport,
                            const std::vector<double>& flowSupport,
                            const std::vector<double>& scnpb, const std::vector<double>& alpha,
                            const std::vector<double>& V_next, int numThreads = 8,
                            const std::string& traversalOrder = "lexicographic", int chunkSize = 32,
                            SEXP stateSubset = R_NilValue) {
  Rcpp::XPtr<ProblemData> pd(problem_ptr);
  if (!pd) {
    Rcpp::stop("problem_ptr is NULL or expired");
  }
  StateOrder order =
      (traversalOrder == "hilbert") ? StateOrder::Hilbert : StateOrder::Lexicographic;
  std::vector<int> subset;
  if (stateSubset != R_NilValue && Rf_length(stateSubset) > 0) {
    const int* p = INTEGER(stateSubset);
    const int n = Rf_length(stateSubset);
    subset.resize(n);
    std::transform(p, p + n, subset.begin(), [](int x) { return x - 1; });
  }
  return bellmanUpdateImpl(*pd, t, stateSupport, flowSupport, scnpb, alpha, V_next, numThreads,
                           order, chunkSize, subset);
}

// ─────────────────────────────────────────────────────────────────────────────
// simulateStepImpl — LP-accurate one-step forward simulation for RTDP Phase 2.
//
// Given a (state, scenario) pair, solves nAdx LPs (one per action) with warm-
// starting, picks the greedy-optimal action under V_next, and returns the
// resulting transport cost, next state index, and chosen action.
//
// Scenario decomposition from flat kdx (0-based):
//   k1 = kdx % nQdx          — index into inflowIndices[]
//   k2 = (kdx / nQdx) % nDdx — index into outflowIndices[]
//   k3 = kdx / (nQdx * nDdx) — index into spotRateIndices[]
// This follows directly from the flowKeys mixed-radix encoding in dp_config.
// ─────────────────────────────────────────────────────────────────────────────
[[nodiscard]] static Rcpp::List
simulateStepImpl(const ProblemData& pd, int t, int state_idx, int scenario_kdx,
                 const std::vector<double>& stateSupport,
                 const std::vector<double>& flowSupport,
                 const std::vector<double>& alpha,
                 const std::vector<double>& V_next) {

  const int nInventoryLevels = pd.R + 1;
  const int nFlowLevels      = pd.nQ;
  const int nSpotRates       = pd.nW;
  const int nOrigins         = pd.nI;
  const int nDestinations    = pd.nJ;
  const int nStrategicCarriers = pd.nCS;
  const int nSpotCarriers    = pd.nCO;
  const int nServices        = pd.nL_;
  const int nLanes           = pd.nL;
  const int nActions         = pd.R + 1;
  const int storageLimit     = pd.R;

  const auto& fromOrigin    = pd.fromOrigin;
  const auto& toDestination = pd.toDestination;
  const auto& Cb            = pd.Cb;
  const auto& Co            = pd.Co;
  const auto& bids          = pd.bids;
  const auto& lanes         = pd.lanes;
  const auto& winnerKeys    = pd.winnerKeys;
  const auto& spotRates     = pd.spotRates;
  const auto& spotRateSupport = pd.spotRateSupport;
  const auto& stateKeys     = pd.stateKeys;
  const auto& winners       = pd.winners;
  const auto& CTb           = pd.CTb;
  const auto& carrierIdx    = pd.carrierIdx;

  const int n = nServices + nSpotCarriers * nLanes;
  const std::vector<int> nWarehouses = {nOrigins, nDestinations};
  const std::vector<int> limits(nOrigins + nDestinations, storageLimit);

  const std::vector<double> extendedStateSupport = utils::mirrorAndNegateVector(stateSupport);

  // Build state and flow index sets (same as bellmanUpdateImpl)
  std::vector<std::vector<int>> stateSupportStack =
      stackStateIdxVectors(nInventoryLevels, nOrigins, nDestinations);
  const std::vector<std::vector<int>> stateIndices = CartesianProductIntSTL(stateSupportStack);

  std::vector<int> flowIdxSingle = utils::createIndexVector(nFlowLevels);
  std::vector<std::vector<int>> inflowIdxStack;
  utils::appendIndexVectors(inflowIdxStack, flowIdxSingle, nOrigins);
  const std::vector<std::vector<int>> inflowIndices = CartesianProductIntSTL(inflowIdxStack);

  std::vector<std::vector<int>> outflowIdxStack;
  utils::appendIndexVectors(outflowIdxStack, flowIdxSingle, nDestinations);
  const std::vector<std::vector<int>> outflowIndices = CartesianProductIntSTL(outflowIdxStack);

  std::vector<int> spotRateIdx = utils::createIndexVector(nSpotRates);
  std::vector<std::vector<int>> spotRateIdxStack;
  utils::appendIndexVectors(spotRateIdxStack, spotRateIdx, nSpotCarriers);
  const std::vector<std::vector<int>> spotRateIndices = CartesianProductIntSTL(spotRateIdxStack);

  const int nQdx = static_cast<int>(inflowIndices.size());
  const int nDdx = static_cast<int>(outflowIndices.size());
  const int nAdx = nActions;

  // Decompose flat scenario index → (k1, k2, k3)
  const int k1 = scenario_kdx % nQdx;
  const int k2 = (scenario_kdx / nQdx) % nDdx;
  const int k3 = scenario_kdx / (nQdx * nDdx);

  const auto& stateIdx    = stateIndices[state_idx];
  const auto& inflowIdx   = inflowIndices[k1];
  const auto& outflowIdx  = outflowIndices[k2];
  const auto& srIdx       = spotRateIndices[k3];

  // ── Build LP for period t ──────────────────────────────────────────────────
  Highs highs;
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("threads", 1);

  std::vector<int> colIdx = createTransportVars(highs, n, winnerKeys, winners, bids, lanes,
                                                CTb, spotRates[t], nSpotCarriers, nLanes);
  try {
    addCapacityConstraints(highs, colIdx, winnerKeys, winners, bids, carrierIdx,
                           nSpotCarriers, nLanes, Cb[t], Co[t]);
    addStorageLimitConstraints(highs, colIdx, winnerKeys, winners, bids, lanes,
                               nSpotCarriers, nLanes, nWarehouses, limits);
    addVolumeConstraint(highs, colIdx, n, 0); // placeholder; overridden per action below
  } catch (const std::exception& e) {
    Rcpp::stop("simulateStep: error building LP: %s", e.what());
  }

  const int volumeRow = nStrategicCarriers + nSpotCarriers + nOrigins + nDestinations;

  // Destination storage bounds (state-dependent)
  for (int idx = 0; idx < nDestinations; ++idx) {
    const int p = nStrategicCarriers + nSpotCarriers + nOrigins + idx;
    highs.changeRowBounds(p, -kHighsInf,
                          storageLimit - extendedStateSupport[stateIdx[nOrigins + idx]]);
  }

  // Inflow bounds (scenario k1)
  for (int d = 0; d < nOrigins; ++d) {
    const int p = nStrategicCarriers + nSpotCarriers + d;
    highs.changeRowBounds(p, -kHighsInf,
                          stateSupport[stateIdx[d]] + flowSupport[inflowIdx[d]]);
  }

  // Spot rates (scenario k3)
  std::vector<double> spotRatesTmp(nSpotCarriers * nLanes, 0.0);
  for (int d = 0; d < nSpotCarriers; ++d) {
    for (int l = 0; l < nLanes; ++l) {
      spotRatesTmp[d * nLanes + l] = spotRateSupport[srIdx[d]];
    }
  }
  updateSpotRates(highs, colIdx, spotRatesTmp, nServices, nSpotCarriers, nLanes);

  // Hold cost for this state
  double holdCost = 0.0;
  for (int k = 0; k < nOrigins; ++k) {
    holdCost += stateSupport[stateIdx[k]] * alpha[k];
  }
  for (int k = 0; k < nDestinations; ++k) {
    const double sj = extendedStateSupport[stateIdx[nOrigins + k]];
    holdCost += std::max(sj, 0.0) * alpha[nOrigins + k];
    holdCost -= std::min(sj, 0.0) * alpha[nOrigins + nDestinations + k];
  }

  // ── Solve LP for each action; pick greedy-optimal ─────────────────────────
  double Qstar = std::numeric_limits<double>::lowest();
  int    j_star    = 0;
  double cost_star = 0.0;
  int    next_star = state_idx;

  HighsBasis lastBasis;
  bool hasBasis = false;

  std::vector<double> x(n, 0.0);
  std::vector<double> xI(nOrigins, 0.0);
  std::vector<double> xJ(nDestinations, 0.0);
  std::vector<double> nextState(nOrigins + nDestinations, 0.0);

  for (int j = 0; j < nAdx; ++j) {
    highs.changeRowBounds(volumeRow, static_cast<double>(j), static_cast<double>(j));
    if (hasBasis) {
      highs.setBasis(lastBasis);
    }
    highs.run();

    if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
      continue;
    }
    lastBasis = highs.getBasis();
    hasBasis  = true;

    const double cost_j = highs.getObjectiveValue();
    const auto&  sol    = highs.getSolution();

    for (int d = 0; d < n; ++d) {
      x[d] = sol.col_value[colIdx[d]];
    }

    std::fill(xI.begin(), xI.end(), 0.0);
    for (int d = 0; d < nOrigins; ++d) {
      for (const auto& l : fromOrigin[d]) {
        xI[d] += x[l - 1];
      }
    }

    std::fill(xJ.begin(), xJ.end(), 0.0);
    for (int d = 0; d < nDestinations; ++d) {
      for (const auto& l : toDestination[d]) {
        xJ[d] += x[l - 1];
      }
    }

    // Compute next state under (action j, inflow k1, outflow k2)
    for (int d = 0; d < nOrigins; ++d) {
      nextState[d] = std::max(
          std::min(stateSupport[stateIdx[d]] + flowSupport[inflowIdx[d]] - xI[d],
                   static_cast<double>(storageLimit)),
          0.0);
    }
    for (int d = 0; d < nDestinations; ++d) {
      nextState[nOrigins + d] =
          storageLimit +
          std::min(std::max(extendedStateSupport[stateIdx[nOrigins + d]] -
                                flowSupport[outflowIdx[d]] + xJ[d],
                            static_cast<double>(-storageLimit)),
                   static_cast<double>(storageLimit));
    }
    const int nextStateIdx =
        std::inner_product(nextState.begin(), nextState.end(), stateKeys.begin(), 0);

    const double Q_j = -holdCost - cost_j + V_next[nextStateIdx];
    if (Q_j > Qstar) {
      Qstar     = Q_j;
      j_star    = j;
      cost_star = cost_j;
      next_star = nextStateIdx;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("cost")           = cost_star,
      Rcpp::Named("next_state_idx") = next_star,  // 0-based
      Rcpp::Named("action_idx")     = j_star       // 0-based
  );
}

// XPtr wrapper for simulateStepImpl.
// state_idx and scenario_kdx are 1-based from R; converted to 0-based internally.
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Rcpp::List simulateStepPtr(SEXP problem_ptr, int t, int state_idx, int scenario_kdx,
                           const std::vector<double>& stateSupport,
                           const std::vector<double>& flowSupport,
                           const std::vector<double>& alpha,
                           const std::vector<double>& V_next) {
  Rcpp::XPtr<ProblemData> pd(problem_ptr);
  if (!pd) {
    Rcpp::stop("problem_ptr is NULL or expired");
  }
  return simulateStepImpl(*pd, t, state_idx - 1, scenario_kdx - 1,
                          stateSupport, flowSupport, alpha, V_next);
}

// Create a HiGHS model (replaces createGRBmodel)
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
SEXP createHIGHSmodel() {
  auto* highs = new Highs();
  highs->setOptionValue("output_flag", false);
  return Rcpp::XPtr<Highs>(highs, true);
}

template <typename T> std::unordered_map<std::string, std::vector<T>> rListToMap(Rcpp::List rlist) {
  std::unordered_map<std::string, std::vector<T>> result;

  Rcpp::CharacterVector names = rlist.names();

  for (int i = 0; i < static_cast<int>(rlist.size()); ++i) {
    auto name = Rcpp::as<std::string>(names[i]);
    Rcpp::Vector<Rcpp::traits::r_sexptype_traits<T>::rtype> vec = rlist[i];
    result[name] = Rcpp::as<std::vector<T>>(vec);
  }

  return result;
}

// Overload for scalars
template <typename T> std::unordered_map<std::string, T> rListToMapScalar(Rcpp::List rlist) {
  std::unordered_map<std::string, T> result;

  Rcpp::CharacterVector names = rlist.names();

  for (int i = 0; i < static_cast<int>(rlist.size()); ++i) {
    auto name = Rcpp::as<std::string>(names[i]);
    T value = Rcpp::as<T>(rlist[i]);
    result[name] = value;
  }

  return result;
}

// Function to expose createTransportVars to R
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
SEXP createTransportVarsCx(SEXP model_ptr, int n, const std::vector<std::string>& winnerKeys,
                           const Rcpp::List& winners, const std::vector<std::vector<int>>& bids,
                           const std::vector<std::vector<int>>& lanes,
                           const Rcpp::List& contractRates, const std::vector<double>& spotRates,
                           int nSpotCarriers, int nLanes) {

  Rcpp::XPtr<Highs> highs(model_ptr);
  if (!highs) {
    Rcpp::stop("model_ptr is NULL or expired");
  }

  // Check if model already has variables of size n
  int numVars = static_cast<int>(highs->getLp().col_cost_.size());
  if (numVars == n) {
    Rcpp::Rcout << "Model already has an objective vector of size " << n << '\n';
    // Return existing column indices (0..n-1)
    auto* colIdx = new std::vector<int>(n);
    std::iota(colIdx->begin(), colIdx->end(), 0);
    return Rcpp::XPtr<std::vector<int>>(colIdx, true);
  }

  std::unordered_map<std::string, std::vector<int>> winnersCx = rListToMap<int>(winners);
  std::unordered_map<std::string, std::vector<double>> contractRatesCx =
      rListToMap<double>(contractRates);

  auto* colIdx =
      new std::vector<int>(createTransportVars(*highs, n, winnerKeys, winnersCx, bids, lanes,
                                               contractRatesCx, spotRates, nSpotCarriers, nLanes));

  return Rcpp::XPtr<std::vector<int>>(colIdx, true);
}

// Function to add capacity constraints
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void addCapacityConstraintsCx(SEXP model_ptr, SEXP transport_ptr,
                              const std::vector<std::string>& winnerKeys, const Rcpp::List& winners,
                              const std::vector<std::vector<int>>& bids,
                              const Rcpp::List& carrierIdx,
                              const std::vector<int>& strategicCapacities,
                              const std::vector<int>& spotCapacities, int nSpotSources,
                              int nLanes) {

  Rcpp::XPtr<Highs> highs(model_ptr);
  Rcpp::XPtr<std::vector<int>> colIdx(transport_ptr);
  if (!highs) {
    Rcpp::stop("model_ptr is NULL or expired");
  }
  if (!colIdx) {
    Rcpp::stop("transport_ptr is NULL or expired");
  }

  std::unordered_map<std::string, std::vector<int>> winnersCx = rListToMap<int>(winners);
  std::unordered_map<std::string, int> carrierIdxCx = rListToMapScalar<int>(carrierIdx);

  addCapacityConstraints(*highs, *colIdx, winnerKeys, winnersCx, bids, carrierIdxCx, nSpotSources,
                         nLanes, strategicCapacities, spotCapacities);
}

// Function to add storage limit constraints
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void addStorageLimitConstraintsCx(SEXP model_ptr, SEXP transport_ptr,
                                  const std::vector<std::string>& winnerKeys,
                                  const Rcpp::List& winners,
                                  const std::vector<std::vector<int>>& bids,
                                  const std::vector<std::vector<int>>& lanes, int nSpotCarriers,
                                  int nLanes, const std::vector<int>& nWarehouses,
                                  const std::vector<int>& limits) {

  Rcpp::XPtr<Highs> highs(model_ptr);
  Rcpp::XPtr<std::vector<int>> colIdx(transport_ptr);
  if (!highs) {
    Rcpp::stop("model_ptr is NULL or expired");
  }
  if (!colIdx) {
    Rcpp::stop("transport_ptr is NULL or expired");
  }

  std::unordered_map<std::string, std::vector<int>> winnersCx = rListToMap<int>(winners);

  addStorageLimitConstraints(*highs, *colIdx, winnerKeys, winnersCx, bids, lanes, nSpotCarriers,
                             nLanes, nWarehouses, limits);
}

// Function to add volume constraint
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void addVolumeConstraintCx(SEXP model_ptr, SEXP transport_ptr, int n, double At) {
  Rcpp::XPtr<Highs> highs(model_ptr);
  Rcpp::XPtr<std::vector<int>> colIdx(transport_ptr);
  if (!highs) {
    Rcpp::stop("model_ptr is NULL or expired");
  }
  if (!colIdx) {
    Rcpp::stop("transport_ptr is NULL or expired");
  }

  addVolumeConstraint(*highs, *colIdx, n, At);
}

// Function to expose printing of objective vector to R
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void printObjectiveVectorCx(SEXP model_ptr) {
  Rcpp::XPtr<Highs> highs(model_ptr);
  if (!highs) {
    Rcpp::stop("model_ptr is NULL or expired");
  }
  printObjectiveVector(*highs);
}

// Function to expose printing of constraints to R
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void printConstraintsCx(SEXP model_ptr) {
  Rcpp::XPtr<Highs> highs(model_ptr);
  if (!highs) {
    Rcpp::stop("model_ptr is NULL or expired");
  }
  printConstraints(*highs);
}

// Solve a fully-built HiGHS model (built via createHIGHSmodel +
// createTransportVarsCx + addCapacity/StorageLimit/VolumeConstraintCx)
// and return objective value, solution vector, and status.
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Rcpp::List solveModelCx(SEXP model_ptr, SEXP transport_ptr) {
  Rcpp::XPtr<Highs> highs(model_ptr);
  Rcpp::XPtr<std::vector<int>> colIdx(transport_ptr);
  if (!highs) {
    Rcpp::stop("model_ptr is NULL or expired");
  }
  if (!colIdx) {
    Rcpp::stop("transport_ptr is NULL or expired");
  }

  const int n = static_cast<int>(colIdx->size());

  double objval{0.0};
  std::string status{"INITIAL"};
  Rcpp::NumericVector x(n, 0.0);

  try {
    highs->run();

    if (highs->getModelStatus() == HighsModelStatus::kOptimal) {
      status = "OPTIMAL";
      objval = highs->getObjectiveValue();
      const auto& sol = highs->getSolution();
      for (int i = 0; i < n; ++i) {
        x[i] = sol.col_value[(*colIdx)[i]];
      }
    } else if (highs->getModelStatus() == HighsModelStatus::kInfeasible) {
      status = "INFEASIBLE";
      Rcpp::Rcout << "Model is infeasible" << '\n';
    } else if (highs->getModelStatus() == HighsModelStatus::kUnbounded) {
      status = "UNBOUNDED";
      Rcpp::Rcout << "Model is unbounded" << '\n';
    } else {
      status = "OTHER";
      Rcpp::Rcout << "Optimization ended with non-optimal status" << '\n';
    }
  } catch (const std::exception& e) {
    Rcpp::Rcout << "Error: " << e.what() << '\n';
  } catch (...) {
    Rcpp::Rcout << "Exception during optimization" << '\n';
  }

  return Rcpp::List::create(Rcpp::_["objval"] = objval,
                            Rcpp::_["x"]      = x,
                            Rcpp::_["status"] = status);
}

// Optimize the model and return objective and decision in R list
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Rcpp::List optimizeModelFromJSON(const std::string& jsonFile, int t,
                                 const std::vector<double>& spotRates,
                                 const std::vector<int>& storage_limits, int volume) {

  const ProblemData pd = parseProblemData(jsonFile);

  const int n = pd.nL_ + pd.nCO * pd.nL;
  const std::vector<int> nWarehouses = {pd.nI, pd.nJ};

  double objval{0.0};
  std::string status{"INITIAL"};
  Rcpp::NumericVector x(n, 0.0);

  try {
    Highs highs;
    highs.setOptionValue("output_flag", false);

    std::vector<int> colIdx = createTransportVars(highs, n, pd.winnerKeys, pd.winners, pd.bids,
                                                  pd.lanes, pd.CTb, spotRates, pd.nCO, pd.nL);

    addCapacityConstraints(highs, colIdx, pd.winnerKeys, pd.winners, pd.bids, pd.carrierIdx, pd.nCO,
                           pd.nL, pd.Cb[t], pd.Co[t]);
    addStorageLimitConstraints(highs, colIdx, pd.winnerKeys, pd.winners, pd.bids, pd.lanes, pd.nCO,
                               pd.nL, nWarehouses, storage_limits);
    addVolumeConstraint(highs, colIdx, n, volume);

    highs.run();

    if (highs.getModelStatus() == HighsModelStatus::kOptimal) {
      status = "OPTIMAL";
      objval = highs.getObjectiveValue();
      const auto& sol = highs.getSolution();
      for (int i = 0; i < n; ++i) {
        x[i] = sol.col_value[colIdx[i]];
      }
    } else if (highs.getModelStatus() == HighsModelStatus::kInfeasible) {
      status = "INFEASIBLE";
      Rcpp::Rcout << "Model is infeasible" << '\n';
    } else if (highs.getModelStatus() == HighsModelStatus::kUnbounded) {
      status = "UNBOUNDED";
      Rcpp::Rcout << "Model is unbounded" << '\n';
    } else {
      status = "OTHER";
      Rcpp::Rcout << "Optimization ended with non-optimal status" << '\n';
    }
  } catch (const std::exception& e) {
    Rcpp::Rcout << "Error: " << e.what() << '\n';
  } catch (...) {
    Rcpp::Rcout << "Exception during optimization" << '\n';
  }

  return Rcpp::List::create(Rcpp::_["objval"] = objval,
                            Rcpp::_["x"]      = x,
                            Rcpp::_["status"] = status);
}

// Function to update the state index
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
std::vector<int> updateStateIdx(const std::vector<int>& stateIdx, const std::vector<int>& inflowIdx,
                                const std::vector<std::vector<int>>& outflowIndices,
                                const std::vector<double>& stateSupport,
                                const std::vector<double>& extendedStateSupport,
                                const std::vector<double>& flowSupport,
                                const std::vector<double>& xI, const std::vector<double>& xJ,
                                double storageLimit, const std::vector<int>& stateKeys,
                                int nOrigins, int nDestinations) {

  std::vector<int> nextStateIdx(outflowIndices.size(), 0);
  std::vector<double> nextState(nOrigins + nDestinations, 0.0);

  for (int i = 0; i < nOrigins; ++i) {
    nextState[i] = std::max(
        std::min(stateSupport[stateIdx[i]] + flowSupport[inflowIdx[i]] - xI[i], storageLimit), 0.0);
  }

  for (size_t k = 0; k < outflowIndices.size(); ++k) {
    const std::vector<int>& outflowIdx = outflowIndices[k];
    for (int j = 0; j < nDestinations; ++j) {
      nextState[nOrigins + j] =
          storageLimit + std::min(std::max(extendedStateSupport[stateIdx[nOrigins + j]] -
                                               flowSupport[outflowIdx[j]] + xJ[j],
                                           -storageLimit),
                                  storageLimit);
    }
    nextStateIdx[k] = std::inner_product(nextState.begin(), nextState.end(), stateKeys.begin(), 0);
  }

  return nextStateIdx;
}

// RAII guard that redirects C-level stdout to /dev/null for its lifetime.
// Restores the original fd in the destructor — safe even if an exception is
// thrown. C++ code should prefer this over the exported begin/end pair.
struct StdoutSuppressor {
  int saved_fd{-1};

  StdoutSuppressor() {
    fflush(stdout);
    saved_fd = dup(1); // NOLINT(cppcoreguidelines-prefer-member-initializer)
    if (saved_fd == -1) {
      return;
    }
    int devnull = open("/dev/null", O_WRONLY);
    if (devnull == -1) {
      close(saved_fd);
      saved_fd = -1;
      return;
    }
    if (dup2(devnull, 1) == -1) {
      close(devnull);
      close(saved_fd);
      saved_fd = -1;
      return;
    }
    close(devnull);
  }

  ~StdoutSuppressor() {
    if (saved_fd == -1) {
      return;
    }
    fflush(stdout);
    dup2(saved_fd, 1);
    close(saved_fd);
  }

  StdoutSuppressor(const StdoutSuppressor&) = delete;
  StdoutSuppressor& operator=(const StdoutSuppressor&) = delete;
  StdoutSuppressor(StdoutSuppressor&&) = delete;
  StdoutSuppressor& operator=(StdoutSuppressor&&) = delete;
};

// begin_suppress_stdout / end_suppress_stdout are exported to R for use in
// optimal_assignment (R/immediate_cost.R). They delegate to StdoutSuppressor
// so the POSIX logic lives in one place. C++ code should use
// StdoutSuppressor directly.
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
int begin_suppress_stdout() {
  StdoutSuppressor s;
  int fd = s.saved_fd;
  s.saved_fd = -1; // disarm destructor — caller (R) takes ownership of fd
  return fd;
}

//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
void end_suppress_stdout(int saved_fd) {
  if (saved_fd == -1) {
    return;
  }
  fflush(stdout);
  dup2(saved_fd, 1);
  close(saved_fd);
}

#ifdef STANDALONE_BUILD
int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
    return (EXIT_FAILURE);
  }

  std::ifstream file(argv[1]);
  nlohmann::json input;
  file >> input;

  const std::vector<double> flowSupport = input["Q"]["vals"].get<std::vector<double>>();
  const int nOrigins = input["nI"][0];
  const int nDestinations = input["nJ"][0];
  const int nStrategicSources = input["nL_"][0];
  const int nSpotCarriers = input["nCO"][0];
  const int nL = input["nL"][0];
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
    std::vector<int> colIdx = createTransportVars(highs, n, winnerKeys, winners, bids, lanes, CTb,
                                                  spotRates[0], nSpotCarriers, nL);

    addCapacityConstraints(highs, colIdx, winnerKeys, winners, bids, carrierIdx, nSpotCarriers, nL,
                           Cb[0], Co[0]);
    addStorageLimitConstraints(highs, colIdx, winnerKeys, winners, bids, lanes, nSpotCarriers, nL,
                               nWarehouses, limits);
    addVolumeConstraint(highs, colIdx, n, 10);

  } catch (const std::exception& e) {
    cout << "Error: " << e.what() << endl;
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}
#endif // STANDALONE_BUILD
