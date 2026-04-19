## ============================================================
## 18_sddp_vs_rtdp.R — SDDP vs RTDP head-to-head (B4)
##
## Three methods on the same 2×2 network instance:
##   1. Exact DP  — provably optimal baseline, enumeration-based
##   2. RTDP      — trajectory-guided ADP, Hilbert ordering
##   3. SDDP      — Benders cuts via SDDP.jl, copula uncertainty
##
## Instance: 2×2, R=5, tau=4  (nSdx=4,356; exact DP ~6.6 s)
## Initial state: s0 = all-zero inventory
##
## Quality metric:
##   Exact DP / RTDP: V[1, s0]  — estimated expected total cost from s0
##   SDDP:           LB = SDDP.calculate_bound()  (dual lower bound)
##                   UB = mean(mstp_simulate()$obj) (simulation upper bound)
##
## Note on uncertainty models:
##   RTDP / Exact DP use the TLPR discrete Q/D distributions (finite support).
##   SDDP uses a Gaussian copula with Poisson(lambda) marginals, where
##   lambda = weighted mean of Q$vals.  The two models are consistent in
##   expectation but differ in tail behaviour.
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(MSTP)
})

N_THREADS  <- max(1L, parallel::detectCores() - 2L)
SDDP_ITER  <- 500L    # SDDP training iterations (increase for final paper run)
SDDP_TRIALS <- 500L   # OOB simulation trials

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

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
    transport_coef   = unlist(e$CTb_list, use.names = FALSE),
    spot_coef        = as.vector(e$CTo),
    entry_stock_0    = rep(0L, e$nI),
    exit_stock_0     = rep(0L, e$nJ),
    exit_short_0     = rep(0L, e$nJ)
  )
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("18_sddp_vs_rtdp.R — SDDP vs RTDP head-to-head (B4)\n")
cat("Instance: 2x2, R=5, tau=4, seed=42\n")
cat(strrep("=", 70), "\n\n")

## ── Generate TLPR instance ────────────────────────────────────────────────────
tmpf <- tempfile(fileext = ".json")
TLPR::generate_instance(
  nI = 2L, nJ = 2L, tau = 4L,
  nB = 5L, nCS = 10L, nCO = 1L, rate = 0.5,
  seed = 42L, path = tmpf
)

e <- new.env()
jsonlite::fromJSON(tmpf) |> list2env(envir = e)
e$from_i <- as_int_list(e$from_i)
e$to_j   <- as_int_list(e$to_j)
create_model(e)
dp_config(e)

# Initial state s0 = all-zero inventory (lex index 1)
s0_idx <- 1L
cat(sprintf("nI=%d  nJ=%d  R=%d  tau=%d  nSdx=%d  nCS=%d  nCO=%d\n",
            e$nI, e$nJ, e$R, e$tau, e$nSdx, e$nCS, e$nCO))
cat(sprintf("Initial state s0=all-zero: lex index %d\n", s0_idx))

## ── Method 1: Exact DP ────────────────────────────────────────────────────────
cat("\n--- Exact DP ---\n")
t0_dp  <- proc.time()[["elapsed"]]
res_dp <- rolling_dp_ptr(e, tmpf, numThreads = N_THREADS, traversalOrder = "hilbert")
t_dp   <- proc.time()[["elapsed"]] - t0_dp

V_dp_s0 <- res_dp$V[1L, s0_idx]
cat(sprintf("  Time: %.2fs  |  V*(s0): %.4f\n", t_dp, V_dp_s0))

## ── Method 2: RTDP ───────────────────────────────────────────────────────────
cat("\n--- RTDP (Phase 2, n_iter=30, n_traj=100, start from s0) ---\n")
t0_rtdp  <- proc.time()[["elapsed"]]
res_rtdp <- rolling_dp_rtdp_p2(
  e, tmpf,
  n_iter     = 30L,
  n_traj     = 100L,
  s0         = s0_idx,
  epsilon    = 0.2,
  tol        = 1e-3,
  numThreads = N_THREADS,
  seed       = 42L
)
t_rtdp <- proc.time()[["elapsed"]] - t0_rtdp

V_rtdp_s0  <- res_rtdp$V_approx[1L, s0_idx]
pi_agr      <- mean(res_dp$pi_star[1L, ] == res_rtdp$pi_star[1L, ], na.rm = TRUE) * 100
rel_err_s0  <- 100 * abs(V_rtdp_s0 - V_dp_s0) / (abs(V_dp_s0) + 1e-8)
# Mean relative error across visited (non-NA) states
ok          <- is.finite(res_dp$V[1L,]) & is.finite(res_rtdp$V_approx[1L,])
rel_err_all <- 100 * mean(abs(res_rtdp$V_approx[1L, ok] - res_dp$V[1L, ok]) /
                          (abs(res_dp$V[1L, ok]) + 1e-8))
cat(sprintf("  Time: %.2fs  |  V_approx(s0): %.4f  |  rel err@s0: %.2f%%  mean: %.2f%%\n",
            t_rtdp, V_rtdp_s0, rel_err_s0, rel_err_all))
cat(sprintf("  pi_star agreement vs Exact DP (t=1): %.1f%%\n", pi_agr))

## ── Method 3: SDDP ───────────────────────────────────────────────────────────
cat(sprintf("\n--- SDDP (%d iterations, %d OOB trials) ---\n", SDDP_ITER, SDDP_TRIALS))

# Lambda: match TLPR Q distribution mean broadcast to all nI+nJ dimensions
lambda_val <- sum(e$Q$vals * e$Q$prob)
lambda_vec <- rep(lambda_val, e$nI + e$nJ)
cat(sprintf("  lambda=%.3f (mean of TLPR Q distribution) x %d dims\n",
            lambda_val, length(lambda_vec)))

inst   <- tlpr_env_to_mstp_inst(e)
config <- MSTP::mstp_config(inst, lambda = lambda_vec, n_scenarios = 10L)

t0_sddp  <- proc.time()[["elapsed"]]
model    <- MSTP::mstp_train(config, iterations = SDDP_ITER)
t_train  <- proc.time()[["elapsed"]] - t0_sddp

lb <- MSTP::mstp_bound(model)
cat(sprintf("  Train time: %.2fs  |  SDDP lower bound (LB): %.4f\n", t_train, lb))

t0_sim  <- proc.time()[["elapsed"]]
sims    <- MSTP::mstp_simulate(model, config, trials = SDDP_TRIALS)
t_sim   <- proc.time()[["elapsed"]] - t0_sim

ub_mean <- mean(sims$obj)
ub_sd   <- sd(sims$obj)
gap_pct <- 100 * (ub_mean - lb) / (abs(lb) + 1e-8)
cat(sprintf("  Sim time:   %.2fs  |  UB (sim mean ± sd): %.4f ± %.4f  |  gap: %.2f%%\n",
            t_sim, ub_mean, ub_sd, gap_pct))

## ── Summary table ─────────────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("%-22s  %8s  %12s  %12s  %9s\n",
            "Method", "Time (s)", "V(s0) / LB", "UB (sim)", "Err@s0%"))
cat(strrep("-", 70), "\n")
cat(sprintf("%-22s  %8.2f  %12.4f  %12s  %9s\n",
            "Exact DP", t_dp, V_dp_s0, "—", "0.00%"))
cat(sprintf("%-22s  %8.2f  %12.4f  %12s  %9.2f%%\n",
            "RTDP (s0-seeded)", t_rtdp, V_rtdp_s0, "—", rel_err_s0))
cat(sprintf("%-22s  %8.2f  %12.4f  %12.4f  %9s\n",
            sprintf("SDDP (%di)", SDDP_ITER), t_train, lb, ub_mean, "n/a*"))
cat(strrep("-", 70), "\n")
cat(sprintf("RTDP mean rel err (all states visited): %.2f%%\n", rel_err_all))
cat(sprintf("RTDP pi_star agreement vs Exact DP:     %.1f%%\n", pi_agr))
cat(sprintf("SDDP LB-UB gap: %.2f%%  |  Threads: %d\n", gap_pct, N_THREADS))
cat("* SDDP uses Poisson(lambda) copula; TLPR methods use discrete Q/D.\n")
cat("  Cost scales are not directly comparable across noise models.\n")

unlink(tmpf)
invisible(gc(verbose = FALSE))
