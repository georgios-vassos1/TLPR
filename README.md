# TLPR

**T**ransport **L**ane **P**ortfolio **R**eservation — an R package implementing stochastic dynamic programming for the Carrier Selection and Slot Allocation Problem (CSSAP) in drayage / container transport.

The methodology and numerical results are described in:

> arxiv:2505.01808

---

## Overview

The CSSAP models a shipper who must decide, each period, how many container slots to reserve across a mix of contracted strategic carriers and spot-market carriers. Demand and supply are stochastic. The package solves the problem via exact stochastic DP (Bellman backward induction) and a family of approximate DP methods that scale to instances where exact DP is intractable.

**Core pipeline:**

```
generate_instance()          # build a problem instance (auction, routing, stochastics)
    ↓
dp_config()                  # populate env with state/action/scenario index arrays
    ↓
rolling_dp_ptr()             # exact backward induction via XPtr (C++/OpenMP)
  — or —
rolling_dp_vfa()             # separable VFA approximation (sampled states, 1-pass)
  — or —
rolling_dp_rtdp_p2()         # RTDP Phase 2/3: trajectory-guided ADP (scales to |S|>50k)
  — or —
rolling_dp_rtdp_p2_heur()    # same as above but with Lagrangian Bellman backups (GPU-ready)
    ↓
system_transition()          # simulate a policy on sampled scenario paths
    ↓
multiperiod_expansion()      # build multi-period LP
solve_lp()                   # solve with HiGHS
```

**Choosing a solver:**

| Instance size | nSdx | Recommended |
|---|---|---|
| Small | < 3k | `rolling_dp_ptr` (exact, fast) |
| Medium | 3k – 50k | `rolling_dp_vfa` or `rolling_dp_rtdp_p2` |
| Large | > 50k | `rolling_dp_rtdp_p2` or `rolling_dp_rtdp_p2_heur` (< 0.4% cost gap, 1.3–1.5× CPU speedup) |

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
| OpenMP | for parallel LP solves in `bellmanUpdatePtr` / `computeEnvironmentCx` |
| nlohmann/json | header-only; bundled or via Homebrew (`brew install nlohmann-json`) |

On Apple Silicon, set the include path in `src/Makevars`:

```makefile
PKG_CXXFLAGS += -I/opt/homebrew/opt/nlohmann-json/include
```

---

## Quick start

```r
library(TLPR)

# 1. Generate a 2×2 instance (2 origins, 2 destinations, tau=4 periods)
tmpf <- tempfile(fileext = ".json")
generate_instance(nI = 2, nJ = 2, tau = 4, seed = 42, path = tmpf)
env <- new.env()
jsonlite::fromJSON(tmpf) |> list2env(envir = env)
env$from_i <- lapply(env$from_i, as.integer)
env$to_j   <- lapply(env$to_j,   as.integer)
create_model(env)
dp_config(env)

# 2a. Exact DP (feasible when nSdx < ~3k)
dp <- rolling_dp_ptr(env, tmpf, numThreads = 8L)
cat("V[t=1, s=1] =", dp$V[1L, 1L], "\n")

# 2b. RTDP Phase 2/3 (large instances; convergence-based stopping)
r2 <- rolling_dp_rtdp_p2(env, tmpf,
                           n_iter = 30L, n_traj = 50L,
                           tol = 1e-2, patience = 3L,
                           numThreads = 8L, seed = 42L)
cat("Converged in", r2$n_done, "iterations\n")
cat("V[t=1, s=1] ≈", r2$V_approx[1L, 1L], "\n")
```

---

## SDDP pipeline (paper: arXiv:2505.01808)

The demo scripts 20–32 implement the SDDP-based framework described in the paper. They depend on the [MSTP](../MSTP) R package for instance generation and the [SDDP.jl](https://github.com/odow/SDDP.jl) Julia library (accessed via JuliaCall) for training and simulation.

**Core pipeline:**

```
MSTP::generate_instance()       # build MSTP instance (carriers, lanes, cost coefficients,
                                #   lambda-calibrated capacities, spot pricing)
    ↓
MSTP::mstp_config()             # populate Julia-side HyperParams struct
    ↓
MSTP::train_model()             # SDDP.jl training: single-cut Benders, HiGHS solver,
                                #   batched crash-resilient loop with cut checkpointing
    ↓
MSTP::simulate_model()          # out-of-bag evaluation: Gaussian copula draws,
                                #   per-trial SDDP cost
    ↓
MSTP::clairvoyant_lp()          # per-trial perfect-foresight LP lower bound (HiGHS)
    ↓
regret = (SDDP_i - LP*_i) / |LP*_i|   # clairvoyant-LP regret metric
```

**Capacity optimization extension (demo/28, demo/32):**

```
MSTP::capacity_duals()          # extract SDDP stage dual variables (shadow prices)
    ↓
stochastic gradient ∇f = v - E[μ]     # reservation cost minus expected dual
    ↓
projected gradient descent on x        # normalised step, box projection to [0, x̄]
    ↓
MSTP::simulate_model()          # re-evaluate at updated capacities
```

The LP proxy baseline (`demo/32`) optimises capacities via L-BFGS-B on the deterministic mean-scenario LP (via `optim(..., method="L-BFGS-B")` with exact LP duals from `highs_solve()`), serving as a comparison point for the SDDP dual gradient.

**Uncertainty model:**

Container inflows and outflows are drawn from a Gaussian copula with Poisson marginals. The correlation matrix has a two-block structure (within-origins, within-destinations, cross-block):

```r
# Block structure of the (nI + nJ) × (nI + nJ) correlation matrix
Σ = [ ρ_within * J + (1-ρ_within) * I,   ρ_cross * J  ]
    [ ρ_cross * J,                         ρ_within * J + (1-ρ_within) * I ]
```

Scenario draws use the NORTA method: Gaussian samples from `mvtnorm::rmvnorm(Σ)` transformed to Poisson quantiles via `qpois(pnorm(...))`.

---

## Demo scripts (SDDP era, 20–32)

| Script | Description |
|---|---|
| `demo/20_vss.R` | Gain of multi-stage recourse: SDDP vs myopic single-stage policy |
| `demo/21_multi_topology_regret.R` | Clairvoyant-LP regret across topologies (1×1 through 6×6) |
| `demo/24_sddp_rtdp_apples_benchmark.R` | Apples-to-apples cost comparison: SDDP vs RTDP on shared OOB trajectories |
| `demo/25_sensitivity_analysis.R` | Parameter sweeps: λ ∈ {200,700,2000}, ρ_cross ∈ {0,0.2,0.4}, spot mult ∈ {1,2,4} |
| `demo/26_longer_horizons.R` | SDDP on τ ∈ {12,26,52} horizons (closes gap with Schmiedel 2025) |
| `demo/27_large_instances.R` | Scalability showcase: 20×20×100 and 40×40×100 instances |
| `demo/28_sddp_capacity_opt.R` | SDDP dual projected gradient capacity optimisation, 6×6×20, λ=50 |
| `demo/29_paper_figures.R` | Generate all publication figures from cached RDS results |
| `demo/32_lp_proxy_mstp.R` | Head-to-head: LP proxy (+32.3%) vs SDDP duals (−8.4%) vs default on 6×6×20 |

Results are cached under `demo/results/` as `.rds` files. `demo/29_paper_figures.R` reads these and writes figures to `demo/figs/`.

---

## Demo scripts (exact DP / ADP era, 01–14)

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
| `demo/08_xptr_equivalence.R` | Verify XPtr API (`rolling_dp_ptr`) matches JSON-file counterparts |
| `demo/09_hilbert_benchmark.R` | Lexicographic vs Hilbert traversal correctness and timing |
| `demo/10_vfa_benchmark.R` | Separable VFA approximation quality vs exact DP |
| `demo/11_vfa_analysis.R` | Multi-seed, multi-regime VFA error analysis |
| `demo/12_vfa_scaling.R` | VFA execution time vs instance size; 60-second frontier |
| `demo/13_rtdp_benchmark.R` | RTDP Phase 1/2/3: quality, speedup, convergence, large-instance quality |
| `demo/14_capacity_opt_multi.R` | Capacity optimisation across 2×1/1×2/2×2 topologies; §5.2.2 multi-topology table |
| `demo/09_lagrangian_test.R` | Lagrangian heuristic vs HiGHS at LP level; canonical + 16 generated instances |
| `demo/23_heuristic_rtdp_comparison.R` | Heuristic RTDP vs exact RTDP on 2×2 instances (nSdx=53k); cost gap and speedup |
| `demo/rerun_mc_figures.R` | Recompute Table 3, SAA/regret figures; `FIGURES_ONLY=1` to skip computation |
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

### Exact dynamic programming

| Function | Description |
|---|---|
| `rolling_dp_ptr()` | Exact backward induction via XPtr; one `bellmanUpdatePtr` call per period |
| `rolling_dp_cx()` | Memory-efficient DP via `bellmanUpdateCx`; older JSON-file API |
| `dynamic_programming()` | R-level backward induction over a pre-materialised transit table |
| `computeEnvironmentCx()` | C++/OpenMP transit table computation (all periods, parallel LP solves) |
| `bellmanUpdatePtr()` | Single-period Bellman update (C++); accepts optional `stateSubset` for sampling |
| `bellmanUpdateCx()` | Single-period Bellman update, JSON-file API |
| `bellmanUpdateHeurPtr()` | Heuristic Bellman update via Lagrangian dual ascent; same signature as `bellmanUpdatePtr` |
| `bellmanUpdateHeurCx()` | Heuristic Bellman update, JSON-file API |
| `system_transition()` | Simulate a DP policy on sampled scenario paths |

### Approximate DP (ADP)

| Function | Description |
|---|---|
| `rolling_dp_vfa()` | Separable VFA: OLS-fit piecewise-linear approximation, optional K-state sampling |
| `rolling_dp_rtdp()` | RTDP Phase 1: iterative K-state sampling with VFA cache/blend |
| `rolling_dp_rtdp_p2()` | RTDP Phase 2/3: trajectory-guided ADP; LP cost O(n_traj × τ × n_A), independent of \|S\| |
| `rolling_dp_rtdp_p2_heur()` | Heuristic RTDP Phase 2/3: identical to above but Bellman backups use Lagrangian dual ascent instead of HiGHS; GPU-batchable fixed-iteration arithmetic |
| `simulateStepPtr()` | LP-accurate single-step forward simulation; returns greedy-optimal next state and action |
| `vfa_fit()` | Fit separable VFA from (possibly sparse) cached value observations |
| `vfa_eval()` | Evaluate separable VFA over all states |

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

**`rolling_dp_ptr` vs `rolling_dp_cx`:**
Both implement identical Bellman aggregation. Prefer `rolling_dp_ptr` — it uses the XPtr API, avoids redundant JSON parsing on each call, and accepts a `stateSubset` parameter for sampled updates. `rolling_dp_cx` is retained for backward compatibility.

**ADP scaling property (`rolling_dp_rtdp_p2`):**
Each iteration runs `n_traj` forward trajectories via `simulateStepPtr`. LP cost is O(n_traj × τ × n_A) per iteration — independent of |S|. Only trajectory-visited states receive a Bellman update. Convergence is checked using the mean relative change in V_approx[1,] over a `patience`-iteration window; early stopping saves 40–70% of the iteration budget in practice.

**Lagrangian heuristic Bellman backups (`rolling_dp_rtdp_p2_heur`):**
Replaces HiGHS LP calls in Bellman backups with Lagrangian dual ascent + capacitated greedy recovery. Each LP call becomes a fixed-iteration arithmetic loop (no branch-on-solver-status), making it embarrassingly parallel across all state-action-scenario triples — the same batch structure as a GPU warp. On 2×2 instances (nSdx = 53k) validated on CPU: v_delta convergence curves match exact RTDP, out-of-sample cost gap < 0.4% on tight-regime instances, 1.3–1.5× CPU speedup (larger gains expected on GPU). Forward trajectories intentionally retain exact LP (`simulateStepPtr`) to keep visited-state distributions comparable.

**Separable VFA:**
`V(s) ≈ Σ_k v_k(s_k)` where each `v_k` is a piecewise-linear function of component `k` of the state, fitted by OLS from cached exact Bellman values. The VFA fills unvisited states in the blended `V_next = cache (priority) + VFA (fallback)`.

**`solve_lp` interface:** accepts a Gurobi-style model list (`A`, `obj`, `rhs`, `sense`, `modelsense`, `vtype`, `lb`, `ub`). Returns `list(x, objval, status)`.

**`dp_config` must be called before DP:** it populates `env$Sdx`, `env$Adx`, `env$scndx`, `env$scnpb`, `env$flowKeys`, etc.

**OpenMP:** never call `omp_set_num_threads()` in package code — it is process-global. All parallel regions use the `num_threads()` clause instead.

---

## Testing

```r
# Run the full test suite
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
- `vfa_fit()`, `vfa_eval()`, `rolling_dp_vfa()`, `rolling_dp_rtdp()`, and `rolling_dp_rtdp_p2()` are all exported via `NAMESPACE` and available after `library(TLPR)`; no explicit `source()` is needed.
