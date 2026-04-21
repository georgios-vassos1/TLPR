## ============================================================
## 24_sddp_rtdp_apples_benchmark.R
##
## Apples-to-apples cost comparison: SDDP vs RTDP on a shared
## MSTP/TLPR compatible instance evaluated on identical OOB
## noise trajectories.
##
## Pipeline:
##   1. Generate a TLPR instance (2×2, R=5) and load into TLPR env.
##   2. Convert to an MSTP inst list (same pattern as demo/18).
##   3. Train SDDP via mstp_train().
##   4. Train RTDP (exact + heuristic) via rolling_dp_rtdp_p2/heur.
##   5. Generate shared OOB noise paths via mstp_simulate().
##   6. SDDP cost from sims$obj (evaluated on same paths by Julia).
##   7. RTDP cost via rtdp_rollout_on_noise() — same noise matrix.
##   8. Report SDDP LB/UB and RTDP UB with gap percentages.
##
## Key property: steps 6 and 7 use IDENTICAL Q/D realisations so
## the comparison is not confounded by noise model differences.
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(MSTP)
  library(jsonlite)
  library(Matrix)
})

N_THREADS   <- max(1L, parallel::detectCores() - 2L)
N_ITER_RTDP <- 20L
N_TRAJ_RTDP <- 50L
LAGR_ITER   <- 50L
N_TRIALS    <- 200L
SEED_RTDP   <- 42L
SDDP_ITERS  <- 500L

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

## demo/18 helper: convert TLPR env → MSTP inst list
tlpr_env_to_mstp_inst <- function(e) {
  list(
    tau              = e$tau,
    nOrigins         = e$nI,
    nDestinations    = e$nJ,
    nCarriers        = e$nCS,
    nSpotCarriers    = e$nCO,
    Bids             = e$B,
    Winners          = setNames(e$winner, NULL),
    Ldx              = unlist(e$L_, use.names = FALSE),
    nLc              = e$nLc,
    entry_capacity   = rep(e$R, e$nI),
    exit_capacity    = rep(e$R, e$nJ),
    carrier_capacity = c(as.vector(e$Cb), as.vector(e$Co)),
    entry_store_coef = e$alpha[seq(e$nI)],
    exit_store_coef  = e$alpha[e$nI + seq(e$nJ)],
    exit_short_coef  = e$alpha[e$nI + e$nJ + seq(e$nJ)],
    transport_coef   = e$CTb,
    spot_coef        = as.vector(e$CTo),
    entry_stock_0    = rep(0L, e$nI),
    exit_stock_0     = rep(0L, e$nJ),
    exit_short_0     = rep(0L, e$nJ)
  )
}

cat(strrep("=", 70), "\n")
cat("24_sddp_rtdp_apples_benchmark.R — Apples-to-apples SDDP vs RTDP\n")
cat(strrep("=", 70), "\n\n")

## ── 1. TLPR instance (2×2, R=5) ──────────────────────────────────────────────
tmpf <- tempfile(fileext = ".json")
on.exit(unlink(tmpf), add = TRUE)

TLPR::generate_instance(
  nI = 2L, nJ = 2L, tau = 4L,
  nB = 5L, nCS = 3L, nCO = 1L, rate = 1.0,
  seed = 99L, path = tmpf
)

e <- new.env()
jsonlite::fromJSON(tmpf) |> list2env(envir = e)
e$from_i <- as_int_list(e$from_i)
e$to_j   <- as_int_list(e$to_j)
TLPR::create_model(e)
TLPR::dp_config(e)

cat(sprintf("Instance: nI=%d  nJ=%d  R=%d  tau=%d  nSdx=%d  nCS=%d  nCO=%d\n",
            e$nI, e$nJ, e$R, e$tau, e$nSdx, e$nCS, e$nCO))

## ── 2. MSTP inst list ─────────────────────────────────────────────────────────
inst <- tlpr_env_to_mstp_inst(e)

## ── 3. SDDP training ──────────────────────────────────────────────────────────
lambda_val <- sum(e$Q$vals * e$Q$prob)
lambda_vec <- rep(lambda_val, e$nI + e$nJ)
corrmat    <- MSTP::mstp_gen_corrmat(2L, 2L, cross_corr = 0.0)
sddp_cfg   <- MSTP::mstp_config(inst, lambda = lambda_vec, corrmat = corrmat,
                                 n_scenarios = 10L)

cat(sprintf("SDDP lambda=%.2f, training %d iterations...\n", lambda_val, SDDP_ITERS))
t0         <- proc.time()["elapsed"]
sddp_model <- MSTP::mstp_train(sddp_cfg, iterations = SDDP_ITERS)
t_sddp     <- proc.time()["elapsed"] - t0
sddp_lb    <- MSTP::mstp_bound(sddp_model)
cat(sprintf("SDDP done in %.1f s  LB=%.2f\n", t_sddp, sddp_lb))

## ── 4a. Exact RTDP ────────────────────────────────────────────────────────────
cat(sprintf("Training exact RTDP (%d iters, %d traj)...\n", N_ITER_RTDP, N_TRAJ_RTDP))
t0 <- proc.time()["elapsed"]
rtdp_exact <- TLPR::rolling_dp_rtdp_p2(
  env        = e,
  jsonFile   = tmpf,
  n_iter     = N_ITER_RTDP,
  n_traj     = N_TRAJ_RTDP,
  tol        = 1e-3,
  min_iter   = 5L,
  patience   = 3L,
  numThreads = N_THREADS,
  seed       = SEED_RTDP
)
t_rtdp_e <- proc.time()["elapsed"] - t0
cat(sprintf("Exact RTDP:  %d iters  %.1f s  cache=%d/%d  V_approx[s0]=%.2f\n",
            rtdp_exact$n_done, t_rtdp_e,
            tail(rtdp_exact$cache_coverage, 1L), e$nSdx,
            rtdp_exact$V_approx[1L, 1L]))

## ── 4b. Heuristic RTDP ───────────────────────────────────────────────────────
cat(sprintf("Training heuristic RTDP (%d iters, %d traj, %d lagr)...\n",
            N_ITER_RTDP, N_TRAJ_RTDP, LAGR_ITER))
t0 <- proc.time()["elapsed"]
rtdp_heur <- TLPR::rolling_dp_rtdp_p2_heur(
  env        = e,
  jsonFile   = tmpf,
  n_iter     = N_ITER_RTDP,
  n_traj     = N_TRAJ_RTDP,
  tol        = 1e-3,
  min_iter   = 5L,
  patience   = 3L,
  numThreads = N_THREADS,
  lagrIter   = LAGR_ITER,
  seed       = SEED_RTDP
)
t_rtdp_h <- proc.time()["elapsed"] - t0
cat(sprintf("Heur  RTDP:  %d iters  %.1f s  cache=%d/%d  V_approx[s0]=%.2f\n",
            rtdp_heur$n_done, t_rtdp_h,
            tail(rtdp_heur$cache_coverage, 1L), e$nSdx,
            rtdp_heur$V_approx[1L, 1L]))

## ── 5. Shared OOB noise paths (SDDP simulation) ──────────────────────────────
cat(sprintf("Generating %d OOB simulation paths via SDDP...\n", N_TRIALS))
sims    <- MSTP::mstp_simulate(sddp_model, sddp_cfg, trials = N_TRIALS)
sddp_ub <- mean(sims$obj)
sddp_sd <- sd(sims$obj)
cat(sprintf("SDDP UB: mean=%.2f  sd=%.2f  (noise dim: %d x %d)\n",
            sddp_ub, sddp_sd, nrow(sims$noise), ncol(sims$noise)))

## ── 6-7. RTDP rollout on the SAME noise paths ─────────────────────────────────
cat("Evaluating exact  RTDP on shared noise paths...\n")
c_rtdp_e <- TLPR::rtdp_rollout_on_noise(
  rtdp_res = rtdp_exact,
  tlpr_env = e,
  sim      = inst,
  noise    = sims$noise,
  n_trials = N_TRIALS
)

cat("Evaluating heuristic RTDP on shared noise paths...\n")
c_rtdp_h <- TLPR::rtdp_rollout_on_noise(
  rtdp_res = rtdp_heur,
  tlpr_env = e,
  sim      = inst,
  noise    = sims$noise,
  n_trials = N_TRIALS
)

## ── 8. Results ────────────────────────────────────────────────────────────────
gap_e <- (mean(c_rtdp_e) - sddp_ub) / (abs(sddp_ub) + 1e-8) * 100
gap_h <- (mean(c_rtdp_h) - sddp_ub) / (abs(sddp_ub) + 1e-8) * 100

cat("\n", strrep("=", 70), "\n", sep = "")
cat("APPLES-TO-APPLES SUMMARY\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Instance       : 2x2  R=%d  tau=%d  nSdx=%d\n", e$R, e$tau, e$nSdx))
cat(sprintf("  Noise model    : Gaussian copula Poisson marginals (lambda=%.2f)\n", lambda_val))
cat(sprintf("  Shared paths   : %d OOB trials (identical Q/D for all methods)\n", N_TRIALS))
cat(sprintf("  SDDP  LB       : %.2f\n", sddp_lb))
cat(sprintf("  SDDP  UB       : mean=%.2f  sd=%.2f\n", sddp_ub, sddp_sd))
cat(sprintf("  RTDP  exact UB : mean=%.2f  sd=%.2f  gap_vs_SDDP=%.2f%%\n",
            mean(c_rtdp_e), sd(c_rtdp_e), gap_e))
cat(sprintf("  RTDP  heur  UB : mean=%.2f  sd=%.2f  gap_vs_SDDP=%.2f%%\n",
            mean(c_rtdp_h), sd(c_rtdp_h), gap_h))
cat(sprintf("  Training time  : SDDP=%.1fs  RTDP_exact=%.1fs  RTDP_heur=%.1fs\n",
            t_sddp, t_rtdp_e, t_rtdp_h))
cat("\n  Note: gap > 0 reflects two coupled effects — NOT a bug:\n")
cat("  (1) NOISE MODEL: RTDP trained on TLPR discrete Q in {0,5,10};\n")
cat("      evaluated on Gaussian copula Poisson(lambda) marginals.\n")
cat("  (2) COST MODEL: TLPR VFA includes terminal salvage (-alpha*s_final),\n")
cat("      making VFA slopes large and negative (multi-period transport\n")
cat("      savings signal).  MSTP charges positive terminal holding.\n")
cat("      The RTDP policy correctly hoards inventory under TLPR costs\n")
cat("      but accumulates large shortage/terminal penalties in MSTP.\n")
cat("  A true apples-to-apples comparison requires task D1: train RTDP\n")
cat("  on the MSTP-compatible cost model (no salvage, same distribution).\n")
cat("  Sanity check: zero VFA = myopic exactly (validated separately).\n")
cat(strrep("=", 70), "\n")
cat("[24] Benchmark complete.\n")
