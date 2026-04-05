## ============================================================
## 02_dp.R — Backward induction and policy extraction
##
## Runs the stochastic Bellman recursion over the transit table
## (loaded from cache or computed via 01_transit.R).
## Evaluates the policies by sampling scenarios from the
## empirical distribution (stochastic evaluation).
##
## Outputs (in parent / .GlobalEnv):
##   V        — (tau+1) x nSdx value function
##   Q        — tau x (nSdx*nAdx) action-value matrix
##   pi_star  — deterministic greedy policy
##   pi_rand  — softmax stochastic policy
## ============================================================

suppressPackageStartupMessages(library(TLPR))

## -- Config / dependencies -----------------------------------
.demo_dir <- local({
  ofiles <- Filter(Negate(is.null), lapply(sys.frames(), `[[`, "ofile"))
  if (length(ofiles)) dirname(normalizePath(tail(ofiles, 1L)[[1L]])) else {
    farg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(farg)) dirname(normalizePath(sub("--file=", "", farg[1L]))) else getwd()
  }
})
if (!exists("JSON_PATH")) source(file.path(.demo_dir, "config/instance1x1_4.R"))
if (!exists("transit"))   source(file.path(.demo_dir, "01_transit.R"))

## -- Load cached DP if available -----------------------------
if (file.exists(DP_CACHE)) {
  cat("[02] Loading DP results from cache...\n")
  dp <- readRDS(DP_CACHE)
  V       <- dp$V
  Q       <- dp$Q
  pi_star <- dp$pi_star
  pi_rand <- dp$pi_rand
  rm(dp)
} else {
  cat("[02] Running stochastic backward induction...\n")
  timex <- Sys.time()
  dp    <- dynamic_programming(env, transit)
  elapsed <- as.numeric(Sys.time() - timex, units = "secs")
  cat(sprintf("[02] DP done in %.1f s\n", elapsed))

  V       <- dp$V
  Q       <- dp$Q
  pi_star <- dp$pi_star
  pi_rand <- dp$pi_rand

  saveRDS(dp, DP_CACHE)
  cat(sprintf("[02] DP results saved to: %s\n", DP_CACHE))
}

## -- Report value at S0 = (0, 8) ----------------------------
i0 <- which(
  env$stateSupport[env$Sdx[, 1L]] == S0[1L] &
  env$extendedStateSupport[env$Sdx[, 2L]] == S0[2L])
cat(sprintf("[02] V[t=1, S0=(0,8)] = %.4f\n", V[1L, i0]))

## -- Stochastic simulation: sample N scenario paths ----------
## The DP policy is designed for stochastic evaluation.
## For each run, draw tau scenarios i.i.d. from scnpb.
## tryCatch handles rare infeasible action-scenario combinations.
cat("[02] Stochastic simulation (N=500 sampled scenario paths)...\n")
set.seed(42L)
N <- 500L
run_sim <- function(pi) {
  tryCatch({
    varphidx <- sample(env$nScen, env$tau, replace = TRUE, prob = env$scnpb)
    system_transition(env, transit, pi, varphidx, init_s = S0)$cost
  }, error = function(e) NA_real_)
}
sim_star <- vapply(seq(N), function(.) run_sim(pi_star), numeric(1L))
sim_rand <- vapply(seq(N), function(.) run_sim(pi_rand), numeric(1L))

## Discard infeasible paths
ok_star <- is.finite(sim_star)
ok_rand <- is.finite(sim_rand) & is.finite(sim_star)

cat(sprintf("[02] pi_star: mean=%.2f  sd=%.2f  (feasible: %d/%d)\n",
            mean(sim_star[ok_star]), sd(sim_star[ok_star]), sum(ok_star), N))
cat(sprintf("[02] pi_rand: mean=%.2f  sd=%.2f  (feasible: %d/%d)\n",
            mean(sim_rand[ok_rand]), sd(sim_rand[ok_rand]), sum(ok_rand), N))
cat(sprintf("[02] pi_star vs pi_rand (matched feasible paths): %.1f%% lower cost\n",
            mean(sim_rand[ok_rand] > sim_star[ok_rand]) * 100))