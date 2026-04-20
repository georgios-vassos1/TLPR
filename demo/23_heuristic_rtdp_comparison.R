## ============================================================
## 23_heuristic_rtdp_comparison.R
##
## Compares exact RTDP (HiGHS Bellman backups) against heuristic
## RTDP (Lagrangian Bellman backups) on generated medium instances
## that are too large for full exact DP.
##
## Both methods use identical forward trajectories (exact LP via
## simulateStepPtr) so that visited-state distributions are matched.
## Comparison metrics:
##   1. Convergence speed (v_delta vs iteration)
##   2. Cache coverage (states updated per iteration)
##   3. Out-of-sample simulation cost via post_hoc_simulation
##   4. Wall-clock time
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(jsonlite)
})

TLPR_ROOT  <- path.expand("~/drayage/TLPR")
N_THREADS  <- 4L
LAGR_ITER  <- 50L
N_ITER     <- 15L
N_TRAJ     <- 50L
N_SIM      <- 200L   # out-of-sample trajectories for policy evaluation
SEED_RTDP  <- 42L
SEED_SIM   <- 99L

## ── instances ─────────────────────────────────────────────────────────────────
gen_specs <- list(
  list(label = "2x2_bal_nCS3",   nI=2L, nJ=2L, tau=4L, nB=6L,  nCS=3L, nCO=3L, rate=1.0, regime="balanced",  seed=101L),
  list(label = "2x2_tight_nCS3", nI=2L, nJ=2L, tau=4L, nB=6L,  nCS=3L, nCO=3L, rate=1.0, regime="tight",     seed=102L)
)

## ── helper: simulate out-of-sample cost for an RTDP result ───────────────────
sim_cost <- function(rtdp_res, env, json_path, n_sim, seed) {
  set.seed(seed)
  costs <- numeric(n_sim)
  prob  <- loadProblemDataCx(json_path)
  theta_final <- vfa_fit(rtdp_res$V_cache[, 1L], env)   # re-fit from cached values

  for (rep in seq(n_sim)) {
    s   <- sample(env$nSdx, 1L) - 1L
    tot <- 0.0
    for (t in seq(env$tau)) {
      V_next <- vfa_eval(theta_final, env)
      cached <- rtdp_res$V_cache[, t + 1L]
      known  <- !is.na(cached)
      if (any(known)) V_next[known] <- cached[known]

      kdx_r <- sample(env$nScen, 1L, prob = env$scnpb)
      step  <- simulateStepPtr(
        problem_ptr  = prob,
        t            = t - 1L,
        state_idx    = s + 1L,
        scenario_kdx = kdx_r,
        stateSupport = as.double(env$stateSupport),
        flowSupport  = as.double(env$Q$vals),
        alpha        = env$alpha,
        V_next       = V_next
      )
      tot <- tot + step$cost
      s   <- step$next_state_idx
    }
    # add terminal holding cost
    sdx_vec <- env$Sdx[s + 1L, ]
    term <- -sum(c(
      env$stateSupport[sdx_vec[env$I_]],
      pmax(env$extendedStateSupport[sdx_vec[env$nI + env$J_]], 0L),
     -pmin(env$extendedStateSupport[sdx_vec[env$nI + env$J_]], 0L)
    ) * env$alpha)
    costs[rep] <- tot + term
  }
  costs
}

## ── main loop ────────────────────────────────────────────────────────────────
all_results <- list()

for (spec in gen_specs) {
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)

  label_args <- spec[setdiff(names(spec), "label")]
  do.call(generate_instance, c(label_args, list(path = tmp)))

  env_raw <- fromJSON(tmp)
  env     <- new.env(parent = baseenv())
  for (nm in names(env_raw)) assign(nm, env_raw[[nm]], envir = env)
  dp_config(env)

  cat(sprintf("\n══ %s  [nSdx=%d  nAdx=%d  tau=%d]\n",
              spec$label, env$nSdx, env$nAdx, env$tau))

  ## -- exact RTDP ─────────────────────────────────────────────────────────────
  t0 <- proc.time()["elapsed"]
  rtdp_exact <- rolling_dp_rtdp_p2(
    env        = env,
    jsonFile   = tmp,
    n_iter     = N_ITER,
    n_traj     = N_TRAJ,
    s0         = NULL,
    epsilon    = 0.2,
    tol        = 1e-3,
    min_iter   = 5L,
    patience   = 3L,
    numThreads = N_THREADS,
    seed       = SEED_RTDP
  )
  t_exact <- proc.time()["elapsed"] - t0
  cat(sprintf("  Exact RTDP  : %d iters  %.1f s  cache=%d/%d\n",
              rtdp_exact$n_done, t_exact,
              tail(rtdp_exact$cache_coverage, 1L), env$nSdx))

  ## -- heuristic RTDP ─────────────────────────────────────────────────────────
  t0 <- proc.time()["elapsed"]
  rtdp_heur <- rolling_dp_rtdp_p2_heur(
    env        = env,
    jsonFile   = tmp,
    n_iter     = N_ITER,
    n_traj     = N_TRAJ,
    s0         = NULL,
    epsilon    = 0.2,
    tol        = 1e-3,
    min_iter   = 5L,
    patience   = 3L,
    numThreads = N_THREADS,
    lagrIter   = LAGR_ITER,
    seed       = SEED_RTDP
  )
  t_heur <- proc.time()["elapsed"] - t0
  cat(sprintf("  Heuristic   : %d iters  %.1f s  cache=%d/%d  speedup=%.1fx\n",
              rtdp_heur$n_done, t_heur,
              tail(rtdp_heur$cache_coverage, 1L), env$nSdx,
              t_exact / max(t_heur, 0.01)))

  ## -- convergence comparison ─────────────────────────────────────────────────
  ne <- length(rtdp_exact$v_delta)
  nh <- length(rtdp_heur$v_delta)
  cat(sprintf("  v_delta (last 3): exact [%.4f, %.4f, %.4f]  heur [%.4f, %.4f, %.4f]\n",
              rtdp_exact$v_delta[max(1L, ne-2L)],
              rtdp_exact$v_delta[max(1L, ne-1L)],
              rtdp_exact$v_delta[ne],
              rtdp_heur$v_delta[max(1L, nh-2L)],
              rtdp_heur$v_delta[max(1L, nh-1L)],
              rtdp_heur$v_delta[nh]))

  ## -- out-of-sample cost comparison ─────────────────────────────────────────
  cat(sprintf("  Simulating %d trajectories for each policy...\n", N_SIM))

  c_exact <- sim_cost(rtdp_exact, env, tmp, N_SIM, SEED_SIM)
  c_heur  <- sim_cost(rtdp_heur,  env, tmp, N_SIM, SEED_SIM)

  gap_abs <- mean(c_heur) - mean(c_exact)
  gap_rel <- gap_abs / (abs(mean(c_exact)) + 1e-8)

  cat(sprintf("  Sim cost  exact : mean=%.2f  sd=%.2f\n", mean(c_exact), sd(c_exact)))
  cat(sprintf("  Sim cost  heur  : mean=%.2f  sd=%.2f\n", mean(c_heur),  sd(c_heur)))
  cat(sprintf("  Cost gap        : %.4f  (%.4f%%)\n", gap_abs, 100 * gap_rel))

  all_results[[spec$label]] <- list(
    label      = spec$label,
    nSdx       = env$nSdx,
    t_exact    = t_exact,
    t_heur     = t_heur,
    iters_exact = rtdp_exact$n_done,
    iters_heur  = rtdp_heur$n_done,
    cache_exact = tail(rtdp_exact$cache_coverage, 1L),
    cache_heur  = tail(rtdp_heur$cache_coverage,  1L),
    mean_exact  = mean(c_exact),
    mean_heur   = mean(c_heur),
    gap_rel     = 100 * gap_rel
  )
}

## ── global summary ────────────────────────────────────────────────────────────
cat("\n=== SUMMARY ===\n")
cat(sprintf("%-22s  %7s  %7s  %7s  %7s\n",
            "Instance", "t_exact", "t_heur", "speedup", "gap%"))
for (r in all_results) {
  cat(sprintf("%-22s  %7.1f  %7.1f  %7.1fx  %7.4f%%\n",
              r$label, r$t_exact, r$t_heur,
              r$t_exact / max(r$t_heur, 0.01),
              r$gap_rel))
}

cat("\n[23] Heuristic RTDP comparison complete.\n")
