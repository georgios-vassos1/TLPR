# TLPR

**T**ransport **L**ane **P**ortfolio **R**eservation — an R package implementing stochastic dynamic programming for the Carrier Selection and Slot Allocation Problem (CSSAP) in drayage / container transport.

The methodology and numerical results are described in:

> arxiv:2505.01808

---

## Overview

The CSSAP models a shipper who must decide, each period, how many container slots to reserve across a mix of contracted strategic carriers and spot-market carriers. Demand and supply are stochastic. The package solves the problem via exact stochastic DP (Bellman backward induction) and evaluates capacity-reservation policies against a multi-period LP formulation.

**Core pipeline:**

```
generate_instance()          # build a problem instance (auction, routing, stochastics)
    ↓
dp_config()                  # populate env with state/action/scenario index arrays
    ↓
rolling_dp_cx()              # backward induction (C++/OpenMP, memory-efficient)
  — or —
computeEnvironmentCx()       # materialise full transit table (R+C++)
dynamic_programming()        # backward induction over transit table (R)
    ↓
system_transition()          # simulate a policy on sampled scenario paths
    ↓
multiperiod_expansion()      # build multi-period LP
solve_lp()                   # solve with HiGHS
```

---

## Installation

```r
# From source (requires C++17 and a working Rcpp toolchain)
install.packages("remotes")
remotes::install_local("path/to/TLPR")

# Or with devtools
devtools::install("path/to/TLPR")
```

**System requirements:**

| Dependency | Notes |
|---|---|
| C++17 compiler | clang++ (macOS) or g++ (Linux) |
| [HiGHS](https://highs.dev) | via the `highs` R package |
| OpenMP | for parallel LP solves in `computeEnvironmentCx` / `bellmanUpdateCx` |
| nlohmann/json | header-only; bundled or via Homebrew (`brew install nlohmann-json`) |

On Apple Silicon, set the include path in `src/Makevars`:

```makefile
PKG_CXXFLAGS += -I/opt/homebrew/opt/nlohmann-json/include
```

---

## Quick start

```r
library(TLPR)

# 1. Generate a 1×1 instance (1 origin, 1 destination, tau=4 periods)
env <- generate_instance(nI = 1, nJ = 1, tau = 4, seed = 42)

# 2. Configure DP index structures
dp_config(env)

# 3. Solve via memory-efficient rolling Bellman updates (C++)
dp <- rolling_dp_cx(env, path_to_json, numThreads = 8L)

# 4. Inspect value function at initial state S0 = (0, 8)
i0 <- which(env$stateSupport[env$Sdx[, 1L]] == 0L &
            env$extendedStateSupport[env$Sdx[, 2L]] == 8L)
cat("V[t=1, S0] =", dp$V[1L, i0])
```

---

## Demo scripts

Run the numbered pipeline in order; each script caches its output for the next.

| Script | Description |
|---|---|
| `demo/01_transit.R` | Compute and cache the full transit matrix (Cx + Rx cross-check) |
| `demo/02_dp.R` | Backward induction; stochastic policy simulation |
| `demo/03_capacity_opt.R` | Optimise carrier capacities; reproduce Table 2 / Table 3 |
| `demo/04_simulation.R` | Compare myopic, random, heuristic, and optimal policies |
| `demo/05_portfolio.R` | Evaluate portfolio contract over action grid |
| `demo/06_dp_scalability.R` | Benchmark DP across state-space sizes and topologies |
| `demo/07_rolling_dp.R` | Compare `dynamic_programming` vs `rolling_dp_cx` (value cross-check) |
| `demo/reproduce_paper.R` | End-to-end reproduction of all paper assertions |

Shared instance configuration lives in `demo/config/instance1x1_4.R`.

---

## Key functions

### Instance generation

| Function | Description |
|---|---|
| `generate_instance()` | Generate a complete TLPR instance; optionally write to JSON |
| `generate_cssap()` | Generate auction structure and routing (lower level) |
| `dp_config()` | Populate `env` with `Sdx`, `Adx`, `scndx`, `scnpb`, `flowKeys`, etc. |

### Dynamic programming

| Function | Description |
|---|---|
| `rolling_dp_cx()` | Memory-efficient DP: one `bellmanUpdateCx` call per period; O(nSdx × nAdx) working memory |
| `dynamic_programming()` | R-level backward induction over a pre-materialised transit table |
| `computeEnvironmentCx()` | C++/OpenMP transit table computation (all periods, parallel LP solves) |
| `bellmanUpdateCx()` | Single-period Bellman update in C++ (called by `rolling_dp_cx`) |
| `system_transition()` | Simulate a DP policy on sampled scenario paths |

### LP / assignment

| Function | Description |
|---|---|
| `solve_lp()` | Solve a fully-specified LP/MIP via HiGHS; Gurobi-format model list |
| `optimal_assignment()` | Single-period optimal slot assignment |
| `heuristic_assignment()` | Greedy feasibility-respecting heuristic |
| `create_model()` | Build the single-period constraint matrix |
| `multiperiod_expansion()` | Expand single-period LP to multi-period block-diagonal form |

### Simulation and evaluation

| Function | Description |
|---|---|
| `simulate_system()` | Multi-policy simulation with exogenous Q/D sampling |
| `eval_portfolio()` | State-dependent portfolio contract evaluation |
| `eval_stateless_portfolio()` | Stateless (no inventory) portfolio evaluation |
| `post_hoc_simulation()` | Recover inventory and allocation trajectories from LP solution |

---

## Architecture notes

**Transit matrix layout** (`tau × nSdx × nAdx × nScen × 6`):

| Column | Content |
|---|---|
| `[0]` | next state index (1-based) |
| `[1]` | transport cost |
| `[2]` | state index (1-based) |
| `[3]` | action index (1-based) |
| `[4]` | scenario index (1-based) |
| `[5]` | period (1-based) |

**`rolling_dp_cx` vs `dynamic_programming`:**
Prefer `rolling_dp_cx` for large instances — it avoids materialising the full `tau × nSdx × nAdx × nScen × 6` transit table. Both paths implement identical Bellman aggregation (raw `scnpb` partial sums over feasible scenarios; infeasible scenarios contribute 0 probability mass, penalising risky actions).

**`solve_lp` interface:** accepts a Gurobi-style model list (`A`, `obj`, `rhs`, `sense`, `modelsense`, `vtype`, `lb`, `ub`). Returns `list(x, objval, status)`.

**`dp_config` must be called before DP:** it populates `env$Sdx`, `env$Adx`, `env$scndx`, `env$scnpb`, `env$flowKeys`, etc.

**OpenMP:** never call `omp_set_num_threads()` in package code — it is process-global. All parallel regions use the `num_threads()` clause instead.

---

## Testing

```r
# Run the full test suite (177 tests)
devtools::test()
# or
Rscript -e "testthat::test_local('.')"
```

| Test file | Coverage |
|---|---|
| `test-dp.R` | DP results, Bellman convergence, `rolling_dp_cx` vs `dynamic_programming` |
| `test-immediate-cost.R` | `h.t`, `heuristic_assignment`, `capacitated_random_assignment` |
| `test-utils.R` | `CartesianProduct*`, `rmvnorm`, `convertListToMap` |
| `test-cssap-gen.R` | `generate_cssap`, `dp_config`, `init_env` |

---

## Known gotchas

- `NAMESPACE` is not fully `roxygen2`-managed — `Rcpp::compileAttributes()` regenerates only the `RcppExports.*` portion. After removing a `[[Rcpp::export]]` annotation, manually remove the corresponding `export()` line from `NAMESPACE`.
- `begin_suppress_stdout` / `end_suppress_stdout` redirect fd 1 to `/dev/null` to silence the HiGHS banner. They return an `int` fd; `end_suppress_stdout(-1)` is a no-op.
- `CartesianProductIntParallel` and `CartesianProductIntParallelxLB` use `std::thread` / `std::async` (not OpenMP); thread count is a direct parameter.
- Instance JSON files must be written via `generate_instance(path = ...)` after any change to problem parameters — C++ functions read directly from the JSON at call time.
