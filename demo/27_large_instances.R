## ============================================================
## 27_large_instances.R — Large-instance SDDP showcase (P1.4)
##
## Demonstrates SDDP scalability on instances that exceed
## Schmiedel (2025) in network width and carrier count:
##
##   20×20×100  (nI=20, nJ=20, nCS=100, τ=12)
##   40×40×100  (nI=40, nJ=40, nCS=100, τ=12)
##
## Schmiedel's largest: 10×50, 3 strategic + ∞ spot, τ ∈ {12,26,52}.
## Our instances have 33× more contracted carriers and larger
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
N_TRIALS    <- 50L    # reduced for large instances (simulation is also expensive)
N_SCENARIOS <- 10L
LAMBDA_VAL  <- 700.0
CROSS_CORR  <- 0.0

## Instance configurations: topology parameters and SDDP iteration budget.
## Large instances (20×20×100, 40×40×100) use reduced iteration counts because
## each SDDP iteration is expensive at this scale.
## 50 and 20 iterations are sufficient to demonstrate convergence and measure
## per-iteration runtime; gap% reflects early-stage convergence, not asymptotic.
INSTANCE_CONFIGS <- list(
  list(label="6×6×20",   nI=6L,  nJ=6L,  nCS=20L,  nIter=500L),
  list(label="20×20×100", nI=20L, nJ=20L, nCS=100L, nIter=50L),
  list(label="40×40×100", nI=40L, nJ=40L, nCS=100L, nIter=20L)
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
  ## Retry on JuliaCall numerical errors (Julia restarts between attempts)
  inst_result <- NULL
  for (.attempt in 1:6L) {
    inst_result <- tryCatch({
      t0    <- proc.time()[["elapsed"]]
      model <- MSTP::mstp_train(config, iterations = cfg$nIter)
      t_tr  <- proc.time()[["elapsed"]] - t0

      lb <- MSTP::mstp_bound(model)

      t0   <- proc.time()[["elapsed"]]
      sims <- MSTP::mstp_simulate(model, config, trials = N_TRIALS)
      t_si <- proc.time()[["elapsed"]] - t0

      list(model = model, lb = lb, sims = sims, t_tr = t_tr, t_si = t_si)
    }, error = function(e) {
      cat(sprintf("  [attempt %d failed: %s — retrying]\n", .attempt, conditionMessage(e)))
      NULL
    })
    if (!is.null(inst_result)) break
  }
  if (is.null(inst_result)) stop(cfg$label, " failed after 6 attempts")
  model  <- inst_result$model
  lb     <- inst_result$lb
  sims   <- inst_result$sims
  t_tr   <- inst_result$t_tr
  t_si   <- inst_result$t_si

  cat(sprintf("  Train: %.1fs  |  LB = %.2f\n", t_tr, lb))
  cat(sprintf("  Simulating %d OOB trials...\n", N_TRIALS))
  cat(sprintf("  Sim:   %.1fs\n", t_si))

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
cat("  Carrier count: our 100 contracted vs Schmiedel's 3 strategic\n")
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

## ── Persist results ───────────────────────────────────────────────────────────
dir.create("demo/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, "demo/results/27_results.rds")
cat("\nResults saved to demo/results/27_results.rds\n")

invisible(gc(verbose = FALSE))
