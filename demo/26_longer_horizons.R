## ============================================================
## 26_longer_horizons.R — Longer-horizon SDDP runs (P1.3)
##
## Trains SDDP on 6×6×20 instances with planning horizons
##   τ ∈ {12, 26, 52}
## to close the horizon gap with Schmiedel (2025) who tests
## up to τ = 52 periods.
##
## Instances generated fresh via MSTP::generate_instance():
##   nOrigins=6, nDestinations=6, nCarriers=20, nSpotCarriers=20
##
## Metrics per horizon:
##   - LB  : Benders lower bound after training
##   - UB  : mean cost over 200 OOB trials
##   - Gap%: (UB − LB) / |LB| × 100
##   - LB/period, UB/period: normalised per planning period
##   - Training and simulation runtimes
##
## Addresses: G6 (competitive runtimes on larger instances),
##            comparison to Schmiedel Table 3 (τ ∈ {12,26,52})
## ============================================================

suppressPackageStartupMessages({
  library(MSTP)
})

## ── Parameters ──────────────────────────────────────────────────────────────

SEED        <- 42L
SDDP_ITER   <- 1000L   # same iteration budget for all horizons
N_TRIALS    <- 200L
N_SCENARIOS <- 10L

LAMBDA_VAL <- 700.0   # medium demand rate (≈ Schmiedel mean demand)
CROSS_CORR <- 0.4

TAU_VALS <- c(12L, 26L, 52L)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("26_longer_horizons.R — τ ∈ {12, 26, 52} on 6×6×20\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("SDDP iters=%d  |  OOB trials=%d  |  λ=%.0f  |  ρ_cross=%.1f\n\n",
            SDDP_ITER, N_TRIALS, LAMBDA_VAL, CROSS_CORR))

## ── Main loop ───────────────────────────────────────────────────────────────

results <- vector("list", length(TAU_VALS))

for (ti in seq_along(TAU_VALS)) {
  tau <- TAU_VALS[ti]

  cat(sprintf("\n%s\n[τ = %d]\n%s\n", strrep("=", 60), tau, strrep("=", 60)))

  ## Generate instance with this horizon
  inst <- MSTP::generate_instance(
    tau           = tau,
    nOrigins      = 6L,
    nDestinations = 6L,
    nCarriers     = 20L,
    nSpotCarriers = 20L,
    nBids         = 10L,
    seed          = SEED
  )

  nOD        <- inst$nOrigins + inst$nDestinations
  lambda_vec <- rep(LAMBDA_VAL, nOD)
  set.seed(SEED)
  corrmat    <- MSTP::gen_corrmat(n_blocks   = 2L,
                                   block_size = inst$nOrigins,
                                   cross_corr = CROSS_CORR)

  cat(sprintf("  nI=%d  nJ=%d  nCS=%d  nCSO=%d  tau=%d\n",
              inst$nOrigins, inst$nDestinations,
              inst$nCarriers, inst$nSpotCarriers, tau))
  cat(sprintf("  Corrmat: %d×%d  |  λ=%.0f  |  Iters=%d\n",
              nrow(corrmat), ncol(corrmat), LAMBDA_VAL, SDDP_ITER))

  config <- MSTP::mstp_config(inst,
                               lambda      = lambda_vec,
                               corrmat     = corrmat,
                               n_scenarios = N_SCENARIOS)

  cat(sprintf("  Training SDDP (%d iterations)...\n", SDDP_ITER))
  t0    <- proc.time()[["elapsed"]]
  model <- MSTP::mstp_train(config, iterations = SDDP_ITER)
  t_tr  <- proc.time()[["elapsed"]] - t0

  lb <- MSTP::mstp_bound(model)
  cat(sprintf("  Train: %.1fs  |  LB = %.2f\n", t_tr, lb))

  cat(sprintf("  Simulating %d OOB trials...\n", N_TRIALS))
  t0   <- proc.time()[["elapsed"]]
  sims <- MSTP::mstp_simulate(model, config, trials = N_TRIALS)
  t_si <- proc.time()[["elapsed"]] - t0

  ub      <- mean(sims$obj)
  ub_sd   <- sd(sims$obj)
  gap_pct <- 100 * (ub - lb) / (abs(lb) + 1e-8)

  cat(sprintf("  Sim:   %.1fs  |  UB = %.2f ± %.2f  |  Gap = %.2f%%\n",
              t_si, ub, ub_sd, gap_pct))

  results[[ti]] <- list(
    tau           = tau,
    lb            = lb,
    ub            = ub,
    ub_sd         = ub_sd,
    gap_pct       = gap_pct,
    per_period_lb = lb / tau,
    per_period_ub = ub / tau,
    t_train       = t_tr,
    t_sim         = t_si
  )

  rm(sims, model)
  gc(verbose = FALSE)
}

## ── Summary table ────────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 92), "\n", sep = "")
cat(sprintf("LONGER-HORIZON SDDP  (6×6×20, λ=%.0f, ρ_cross=%.1f, seed=%d)\n",
            LAMBDA_VAL, CROSS_CORR, SEED))
cat(strrep("=", 92), "\n")

hdr <- sprintf("%-5s  %12s  %12s  %10s  %8s  %10s  %10s  %8s  %8s",
               "τ", "LB", "UB", "UB_sd", "Gap%",
               "LB/period", "UB/period", "Train(s)", "Sim(s)")
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")

for (r in results)
  cat(sprintf("%-5d  %12.1f  %12.1f  %10.1f  %7.2f%%  %10.1f  %10.1f  %7.1fs  %7.1fs\n",
              r$tau, r$lb, r$ub, r$ub_sd, r$gap_pct,
              r$per_period_lb, r$per_period_ub, r$t_train, r$t_sim))

cat(strrep("-", nchar(hdr)), "\n")
cat(sprintf("SDDP iters=%d  |  OOB trials=%d  |  Scenarios/iter=%d\n",
            SDDP_ITER, N_TRIALS, N_SCENARIOS))
cat("Schmiedel (2025): τ ∈ {12,26,52}, Normal demand, independent nodes\n")
cat("Gap% = (UB − LB) / |LB| × 100\n")

## ── Horizon scaling note ────────────────────────────────────────────────────

cat("\nHorizon scaling analysis:\n")
if (length(results) >= 2L) {
  r12 <- results[[1]]
  for (rr in results[-1]) {
    ratio_lb  <- rr$per_period_lb / r12$per_period_lb
    ratio_t   <- rr$t_train / r12$t_train
    ratio_tau <- rr$tau / r12$tau
    cat(sprintf("  τ=%2d vs τ=%2d: cost/period ratio=%.2f  train-time ratio=%.2f  (τ ratio=%.2f)\n",
                rr$tau, r12$tau, ratio_lb, ratio_t, ratio_tau))
  }
}
cat("Expected: train-time scales ≈ linearly with τ (more stages = more LP solves).\n")
cat("Cost/period: stable if stationary dynamics, increasing if congestion builds.\n")

## ── LaTeX snippet ─────────────────────────────────────────────────────────

cat("\n%% LaTeX snippet for horizon comparison table (§6.2)\n")
cat(paste0("$\\tau$ & LB & LB/period & UB & UB/period",
           " & Gap\\% & Train (s) \\\\\n"))
cat("\\hline\n")
for (r in results)
  cat(sprintf(paste0("$%d$ & $%.0f$ & $%.0f$ & $%.0f \\pm %.0f$",
                     " & $%.0f$ & $%.2f\\%%$ & $%.0f$ \\\\\n"),
              r$tau, r$lb, r$per_period_lb,
              r$ub, r$ub_sd, r$per_period_ub,
              r$gap_pct, r$t_train))

invisible(gc(verbose = FALSE))
