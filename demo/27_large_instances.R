## ============================================================
## 27_large_instances.R — Large-instance SDDP showcase (P1.4)
##
## Demonstrates SDDP scalability on instances that exceed
## Schmiedel (2025) in network width and carrier count:
##
##   20×20×400  (nI=20, nJ=20, nCS=400, nCSO=400, τ=12)
##   40×40×400  (nI=40, nJ=40, nCS=400, nCSO=400, τ=12)
##
## Schmiedel's largest: 10×50, 3 strategic + ∞ spot, τ ∈ {12,26,52}.
## Our instances have 133× more contracted carriers and larger
## networks, while keeping the stochastic spot rate model.
##
## Instances generated fresh via MSTP::generate_instance().
## Reports: LB, UB, gap%, train and simulation runtimes.
##
## Addresses: G6 (HiGHS + SDDP.jl competitive runtimes),
##            R4 (instances > Schmiedel), R5 (multi-instance study)
## ============================================================

suppressPackageStartupMessages({
  library(MSTP)
})

## ── Parameters ──────────────────────────────────────────────────────────────

SEED        <- 42L
N_TRIALS    <- 100L
N_SCENARIOS <- 10L
LAMBDA_VAL  <- 700.0
CROSS_CORR  <- 0.4

## Instance configurations: topology parameters and SDDP iteration budget
INSTANCE_CONFIGS <- list(
  list(label="6×6×20",   nI=6L,  nJ=6L,  nCS=20L,  nIter=1000L),
  list(label="20×20×400", nI=20L, nJ=20L, nCS=400L, nIter=500L),
  list(label="40×40×400", nI=40L, nJ=40L, nCS=400L, nIter=300L)
)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("27_large_instances.R — Large-instance SDDP showcase\n")
cat(strrep("=", 70), "\n\n")
cat("Schmiedel (2025) Table 2 reference:\n")
cat("  S: 5×25,  3 strategic + ∞ spot,  τ ∈ {12,26,52}\n")
cat("  M: 10×25, 3 strategic + ∞ spot,  τ ∈ {12,26,52}\n")
cat("  L: 10×50, 3 strategic + ∞ spot,  τ ∈ {12,26,52}\n\n")
cat(sprintf("λ=%.0f  |  ρ_cross=%.1f  |  OOB trials=%d\n\n",
            LAMBDA_VAL, CROSS_CORR, N_TRIALS))

## ── Main loop ───────────────────────────────────────────────────────────────

results <- vector("list", length(INSTANCE_CONFIGS))

for (ci in seq_along(INSTANCE_CONFIGS)) {
  cfg <- INSTANCE_CONFIGS[[ci]]

  cat(sprintf("\n%s\n[%s]\n%s\n",
              strrep("=", 60), cfg$label, strrep("=", 60)))

  inst <- MSTP::generate_instance(tau=12L, nOrigins=cfg$nI, nDestinations=cfg$nJ,
                                   nCarriers=cfg$nCS, seed=SEED)

  nI  <- inst$nOrigins
  nJ  <- inst$nDestinations
  nCS <- inst$nCarriers
  tau <- inst$tau
  nOD <- nI + nJ

  cat(sprintf("  nI=%d  nJ=%d  nCS=%d  τ=%d\n", nI, nJ, nCS, tau))

  lambda_vec <- rep(LAMBDA_VAL, nOD)
  corrmat    <- MSTP::gen_corrmat(n_blocks   = 2L,
                                   block_size = nI,
                                   cross_corr = CROSS_CORR)

  cat(sprintf("  Corrmat: %d×%d  |  SDDP iters=%d\n",
              nrow(corrmat), ncol(corrmat), cfg$nIter))

  config <- MSTP::mstp_config(inst,
                               lambda      = lambda_vec,
                               corrmat     = corrmat,
                               n_scenarios = N_SCENARIOS)

  cat(sprintf("  Training SDDP (%d iterations)...\n", cfg$nIter))
  t0    <- proc.time()[["elapsed"]]
  model <- MSTP::mstp_train(config, iterations = cfg$nIter)
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

  cat("  State-space |S|: astronomically large (SDDP does not enumerate)\n")

  results[[ci]] <- list(
    label    = cfg$label,
    nI       = nI,
    nJ       = nJ,
    nCS      = nCS,
    tau      = tau,
    nIter    = cfg$nIter,
    lb       = lb,
    ub       = ub,
    ub_sd    = ub_sd,
    gap_pct  = gap_pct,
    t_train  = t_tr,
    t_sim    = t_si
  )

  rm(sims, model, inst)
  gc(verbose = FALSE)
}

## ── Summary table ────────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 96), "\n", sep = "")
cat("LARGE-INSTANCE SDDP SHOWCASE\n")
cat(strrep("=", 96), "\n")

hdr <- sprintf("%-14s  %4s  %4s  %6s  %4s  %5s  %12s  %12s  %8s  %8s  %8s",
               "Instance", "nI", "nJ", "nCS", "τ", "Iters",
               "LB", "UB", "Gap%", "Train(s)", "Sim(s)")
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")

for (r in results) {
  if (is.null(r)) next
  cat(sprintf("%-14s  %4d  %4d  %6d  %4d  %5d  %12.1f  %12.1f  %7.2f%%  %7.1fs  %7.1fs\n",
              r$label, r$nI, r$nJ, r$nCS, r$tau, r$nIter,
              r$lb, r$ub, r$gap_pct, r$t_train, r$t_sim))
}

cat(strrep("-", nchar(hdr)), "\n")
cat(sprintf("λ=%.0f  |  ρ_cross=%.1f  |  OOB trials=%d  |  Scenarios/iter=%d\n",
            LAMBDA_VAL, CROSS_CORR, N_TRIALS, N_SCENARIOS))
cat("\nKey comparisons vs Schmiedel (2025):\n")
cat("  Network width: our 40×40 > Schmiedel's 10×50 in carrier lanes\n")
cat("  Carrier count: our 400 contracted vs Schmiedel's 3 strategic\n")
cat("  Spot model: our stochastic finite-capacity vs deterministic infinite\n")
cat("  Solver: HiGHS (open-source) vs Gurobi\n")
cat("Gap% = (UB − LB) / |LB| × 100\n")

## ── LaTeX snippet ─────────────────────────────────────────────────────────

cat("\n%% LaTeX snippet for large-instance table (§6.2)\n")
cat(paste0("Instance & $n_I$ & $n_J$ & $n_{\\rm CS}$ & $\\tau$ & Iters",
           " & LB & UB & Gap\\% & Train (s) \\\\\n"))
cat("\\hline\n")
for (r in results) {
  if (is.null(r)) next
  cat(sprintf(paste0("%s & %d & %d & %d & %d & %d",
                     " & $%.0f$ & $%.0f \\pm %.0f$ & $%.2f\\%%$ & $%.0f$ \\\\\n"),
              r$label, r$nI, r$nJ, r$nCS, r$tau, r$nIter,
              r$lb, r$ub, r$ub_sd, r$gap_pct, r$t_train))
}

invisible(gc(verbose = FALSE))
