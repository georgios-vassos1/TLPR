## ============================================================
## 25_sensitivity_analysis.R — Systematic sensitivity analysis (P1.2)
##
## Sweeps key uncertainty model parameters on 6×6×20 instances
## (tau=12, nOrigins=6, nDestinations=6, nCarriers=20):
##
##   Sweep 1 — demand rate λ ∈ {200, 700, 2000}  (ρ_cross = 0.4)
##   Sweep 2 — cross-block correlation ρ_cross ∈ {0.0, 0.2, 0.4}  (λ = 700)
##   Sweep 3 — spot rate multiplier m ∈ {1, 2, 4}  (λ = 700, ρ_cross = 0.4)
##             spot_coef scaled by m relative to base U[3,9]; tests whether
##             spot price risk materially affects policy cost and gap%.
##
## Primary metric: SDDP gap% = (UB − LB) / |LB| × 100
##   UB = mean cost over 200 OOB trials
##   LB = Benders lower bound after 500 SDDP iterations
##
## Parallels Schmiedel (2025) §5.2 sensitivity design:
##   their σ ∈ {60,350,700} maps to our λ ∈ {200,700,2000}
##   their ρ ∈ {0,0.4,0.8}  maps to our ρ_cross ∈ {0.0,0.2,0.4}
##
## Addresses: G1 (Poisson marginals), G2 (correlated demand), G4 (spot risk), R5
## ============================================================

suppressPackageStartupMessages({
  library(MSTP)
})

## ── Parameters ──────────────────────────────────────────────────────────────

SEED        <- 42L
SDDP_ITER   <- 500L
N_TRIALS    <- 200L
N_SCENARIOS <- 10L
N_INST      <- 3L   # number of instances averaged per configuration

## Parameter grids
LAMBDA_VALS     <- c(200, 700, 2000)    # Poisson demand intensity
CROSS_CORR_VALS <- c(0.0, 0.2, 0.4)    # cross-block copula correlation
SPOT_MULT_VALS  <- c(1.0, 2.0, 4.0)    # spot cost multiplier (1x = base U[3,9])

## Fixed values used when the other axes are swept
FIXED_CROSS      <- 0.0
FIXED_LAMBDA     <- 700.0
FIXED_SPOT_MULT  <- 1.0                 # base spot cost level

## ── Generate instances ──────────────────────────────────────────────────────

cat("\n", strrep("=", 70), "\n", sep = "")
cat("25_sensitivity_analysis.R — λ and ρ sweeps on 6×6×20\n")
cat(strrep("=", 70), "\n\n")

instances <- lapply(seq_len(N_INST), function(i)
  MSTP::generate_instance(tau=12L, nOrigins=6L, nDestinations=6L,
                           nCarriers=20L, seed=SEED + i - 1L))
cat(sprintf("Generated %d instance(s) (6×6×20, τ=12, seeds %d–%d)\n\n",
            N_INST, SEED, SEED + N_INST - 1L))

## ── Helper: train + simulate one configuration ──────────────────────────────

run_config <- function(inst, lambda_val, cross_corr, spot_mult = 1.0) {
  inst_run <- inst
  inst_run$spot_coef <- inst$spot_coef * spot_mult

  lambda_vec <- rep(as.numeric(lambda_val),
                    inst$nOrigins + inst$nDestinations)
  corrmat <- MSTP::gen_corrmat(
    n_blocks   = 2L,
    block_size = inst$nOrigins,
    cross_corr = cross_corr
  )

  config <- MSTP::mstp_config(inst_run,
                               lambda      = lambda_vec,
                               corrmat     = corrmat,
                               n_scenarios = N_SCENARIOS)

  t0    <- proc.time()[["elapsed"]]
  model <- MSTP::mstp_train(config, iterations = SDDP_ITER)
  t_tr  <- proc.time()[["elapsed"]] - t0

  lb   <- MSTP::mstp_bound(model)

  t0   <- proc.time()[["elapsed"]]
  sims <- MSTP::mstp_simulate(model, config, trials = N_TRIALS)
  t_si <- proc.time()[["elapsed"]] - t0

  ub      <- mean(sims$obj)
  ub_sd   <- sd(sims$obj)
  gap_pct <- 100 * (ub - lb) / (abs(lb) + 1e-8)

  list(lb = lb, ub = ub, ub_sd = ub_sd, gap_pct = gap_pct,
       t_train = t_tr, t_sim = t_si)
}

aggregate_runs <- function(runs) {
  fields <- c("lb", "ub", "ub_sd", "gap_pct", "t_train")
  mat    <- do.call(rbind, lapply(runs, function(r) unlist(r[fields])))
  as.list(colMeans(mat))
}

## ── Sweep 1: λ variation  (ρ_cross = FIXED_CROSS) ──────────────────────────

cat(strrep("=", 70), "\n")
cat(sprintf("SWEEP 1: λ variation  (ρ_cross = %.1f, n_inst = %d)\n",
            FIXED_CROSS, N_INST))
cat(strrep("=", 70), "\n\n")

res_lambda <- vector("list", length(LAMBDA_VALS))

for (li in seq_along(LAMBDA_VALS)) {
  lv <- LAMBDA_VALS[li]
  cat(sprintf("  λ = %4d ... ", lv))
  flush.console()

  runs <- lapply(instances, run_config,
                 lambda_val = lv,
                 cross_corr = FIXED_CROSS)

  agg <- aggregate_runs(runs)
  res_lambda[[li]] <- c(list(lambda = lv, cross_corr = FIXED_CROSS), agg)

  cat(sprintf("LB=%.0f  UB=%.0f±%.0f  Gap=%.2f%%  Train=%.1fs\n",
              agg$lb, agg$ub, agg$ub_sd, agg$gap_pct, agg$t_train))
}

## ── Sweep 2: ρ_cross variation  (λ = FIXED_LAMBDA) ─────────────────────────

cat(sprintf("\n%s\n", strrep("=", 70)))
cat(sprintf("SWEEP 2: ρ_cross variation  (λ = %.0f, n_inst = %d)\n",
            FIXED_LAMBDA, N_INST))
cat(strrep("=", 70), "\n\n")

res_corr <- vector("list", length(CROSS_CORR_VALS))

for (ci in seq_along(CROSS_CORR_VALS)) {
  cc <- CROSS_CORR_VALS[ci]
  cat(sprintf("  ρ_cross = %.1f ... ", cc))
  flush.console()

  runs <- lapply(instances, run_config,
                 lambda_val = FIXED_LAMBDA,
                 cross_corr = cc)

  agg <- aggregate_runs(runs)
  res_corr[[ci]] <- c(list(lambda = FIXED_LAMBDA, cross_corr = cc), agg)

  cat(sprintf("LB=%.0f  UB=%.0f±%.0f  Gap=%.2f%%  Train=%.1fs\n",
              agg$lb, agg$ub, agg$ub_sd, agg$gap_pct, agg$t_train))
}

## ── Sweep 3: spot multiplier variation  (λ = FIXED_LAMBDA, ρ = FIXED_CROSS) ─

cat(sprintf("\n%s\n", strrep("=", 70)))
cat(sprintf("SWEEP 3: spot rate multiplier  (λ = %.0f, ρ_cross = %.1f, n_inst = %d)\n",
            FIXED_LAMBDA, FIXED_CROSS, N_INST))
cat(sprintf("  base spot_coef ~ U[3,9]; multiplier scales mean spot rate\n"))
cat(sprintf("  mean contracted ≈ 7; spot means at 1x=6, 2x=12, 4x=24\n"))
cat(strrep("=", 70), "\n\n")

res_spot <- vector("list", length(SPOT_MULT_VALS))

for (si in seq_along(SPOT_MULT_VALS)) {
  sm <- SPOT_MULT_VALS[si]
  cat(sprintf("  spot_mult = %.1f (mean spot ≈ %.0f) ... ", sm, sm * 6.0))
  flush.console()

  runs <- lapply(instances, run_config,
                 lambda_val = FIXED_LAMBDA,
                 cross_corr = FIXED_CROSS,
                 spot_mult  = sm)

  agg <- aggregate_runs(runs)
  res_spot[[si]] <- c(list(spot_mult = sm, lambda = FIXED_LAMBDA,
                            cross_corr = FIXED_CROSS), agg)

  cat(sprintf("LB=%.0f  UB=%.0f±%.0f  Gap=%.2f%%  Train=%.1fs\n",
              agg$lb, agg$ub, agg$ub_sd, agg$gap_pct, agg$t_train))
}

## ── Summary tables ──────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 82), "\n", sep = "")
cat(sprintf("SENSITIVITY: λ  (ρ_cross=%.1f, averaged over %d instances)\n",
            FIXED_CROSS, N_INST))
cat(strrep("=", 82), "\n")
cat(sprintf("%-8s  %10s  %10s  %10s  %8s  %8s\n",
            "λ", "LB", "UB", "UB_sd", "Gap%", "Train(s)"))
cat(strrep("-", 82), "\n")
for (r in res_lambda)
  cat(sprintf("%-8d  %10.1f  %10.1f  %10.1f  %7.2f%%  %7.1fs\n",
              r$lambda, r$lb, r$ub, r$ub_sd, r$gap_pct, r$t_train))
cat(strrep("=", 82), "\n")

cat("\n\n", strrep("=", 82), "\n", sep = "")
cat(sprintf("SENSITIVITY: ρ_cross  (λ=%.0f, averaged over %d instances)\n",
            FIXED_LAMBDA, N_INST))
cat(strrep("=", 82), "\n")
cat(sprintf("%-10s  %10s  %10s  %10s  %8s  %8s\n",
            "ρ_cross", "LB", "UB", "UB_sd", "Gap%", "Train(s)"))
cat(strrep("-", 82), "\n")
for (r in res_corr)
  cat(sprintf("%-10.1f  %10.1f  %10.1f  %10.1f  %7.2f%%  %7.1fs\n",
              r$cross_corr, r$lb, r$ub, r$ub_sd, r$gap_pct, r$t_train))
cat(strrep("=", 82), "\n")

cat("\n\n", strrep("=", 82), "\n", sep = "")
cat(sprintf("SENSITIVITY: spot multiplier  (λ=%.0f, ρ_cross=%.1f, averaged over %d instances)\n",
            FIXED_LAMBDA, FIXED_CROSS, N_INST))
cat(strrep("=", 82), "\n")
cat(sprintf("%-10s  %10s  %10s  %10s  %10s  %8s  %8s\n",
            "spot_mult", "mean_spot", "LB", "UB", "UB_sd", "Gap%", "Train(s)"))
cat(strrep("-", 82), "\n")
for (r in res_spot)
  cat(sprintf("%-10.1f  %10.1f  %10.1f  %10.1f  %10.1f  %7.2f%%  %7.1fs\n",
              r$spot_mult, r$spot_mult * 6.0,
              r$lb, r$ub, r$ub_sd, r$gap_pct, r$t_train))
cat(strrep("=", 82), "\n")

cat(sprintf("\nInstance: 6×6×20 (tau=12)  |  SDDP iters=%d  |  OOB trials=%d  |  n_inst=%d\n",
            SDDP_ITER, N_TRIALS, N_INST))
cat("Copula: Gaussian 2-block corrmat, Poisson(λ) marginals\n")
cat("Gap% = (UB − LB) / |LB| × 100   (smaller = tighter bound)\n")

## ── LaTeX snippet ────────────────────────────────────────────────────────────

cat("\n%% LaTeX snippet for sensitivity table (§6.4)\n\n")

cat("%% Sweep 1: λ variation\n")
cat("$\\lambda$ & LB & UB ($\\pm$ sd) & Gap\\% & Train (s) \\\\\n")
cat("\\hline\n")
for (r in res_lambda)
  cat(sprintf("$%d$ & $%.0f$ & $%.0f \\pm %.0f$ & $%.2f\\%%$ & $%.0f$ \\\\\n",
              r$lambda, r$lb, r$ub, r$ub_sd, r$gap_pct, r$t_train))

cat("\n%% Sweep 2: ρ_cross variation\n")
cat("$\\rho_{\\rm cross}$ & LB & UB ($\\pm$ sd) & Gap\\% & Train (s) \\\\\n")
cat("\\hline\n")
for (r in res_corr)
  cat(sprintf("$%.1f$ & $%.0f$ & $%.0f \\pm %.0f$ & $%.2f\\%%$ & $%.0f$ \\\\\n",
              r$cross_corr, r$lb, r$ub, r$ub_sd, r$gap_pct, r$t_train))

cat("\n%% Sweep 3: spot multiplier variation\n")
cat("$m$ & $\\bar{c}_{\\rm spot}$ & LB & UB ($\\pm$ sd) & Gap\\% & Train (s) \\\\\n")
cat("\\hline\n")
for (r in res_spot)
  cat(sprintf("$%.0f\\times$ & $%.0f$ & $%.0f$ & $%.0f \\pm %.0f$ & $%.2f\\%%$ & $%.0f$ \\\\\n",
              r$spot_mult, r$spot_mult * 6.0,
              r$lb, r$ub, r$ub_sd, r$gap_pct, r$t_train))

## ── Persist results ───────────────────────────────────────────────────────────
dir.create("demo/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(lambda = res_lambda, corr = res_corr, spot = res_spot),
        "demo/results/25_results.rds")
cat("\nResults saved to demo/results/25_results.rds\n")

invisible(gc(verbose = FALSE))
