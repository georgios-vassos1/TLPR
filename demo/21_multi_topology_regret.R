## ============================================================
## 21_multi_topology_regret.R — Multi-topology SDDP regret (B1)
##
## Trains SDDP on multiple network topologies, simulates OOB
## trials, and computes per-trial regret against the clairvoyant
## (perfect-foresight) LP:
##
##   regret_i = (SDDP_cost_i - LP*_i) / |LP*_i|
##
## Topologies (same as demo/14):
##   2x1 R=5    (|S| =    396)
##   1x2 R=5    (|S| =    726)
##   2x2 R=5    (|S| =  4,356)
##   2x2 R=10   (|S| = 53,361)
##
## Outputs:
##   - Regret summary table (mean, sd, quantiles per topology)
##   - Q-Q quantiles for cross-topology regret comparison
##   - LaTeX snippet for paper §6
##
## Addresses: G1 (optimality benchmark), G2 (correlated demand),
##            R5 (multiple instances)
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(MSTP)
})

## ── Parameters ──────────────────────────────────────────────────────────────

SEED        <- 42L
SDDP_ITER   <- 1000L
N_TRIALS    <- 500L
N_SCENARIOS <- 10L

CONFIGS <- list(
  list(nI = 2L, nJ = 1L, rate = 0.5, label = "2x1 R=5"),
  list(nI = 1L, nJ = 2L, rate = 0.5, label = "1x2 R=5"),
  list(nI = 2L, nJ = 2L, rate = 0.5, label = "2x2 R=5"),
  list(nI = 2L, nJ = 2L, rate = 1.0, label = "2x2 R=10")
)

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

## ── Adapter (CTb-fix, matches demo/18-20) ───────────────────────────────────

tlpr_env_to_mstp_inst <- function(e) {
  # as.array() ensures length-1 vectors stay as Julia Vector, not scalar
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
    entry_capacity   = as.array(rep(e$R, e$nI)),
    exit_capacity    = as.array(rep(e$R, e$nJ)),
    carrier_capacity = c(as.vector(e$Cb), as.vector(e$Co)),
    entry_store_coef = as.array(e$alpha[seq(e$nI)]),
    exit_store_coef  = as.array(e$alpha[e$nI + seq(e$nJ)]),
    exit_short_coef  = as.array(e$alpha[e$nI + e$nJ + seq(e$nJ)]),
    transport_coef   = e$CTb,
    spot_coef        = as.vector(e$CTo),
    entry_stock_0    = as.array(rep(0L, e$nI)),
    exit_stock_0     = as.array(rep(0L, e$nJ)),
    exit_short_0     = as.array(rep(0L, e$nJ))
  )
}

## ── Clairvoyant LP (generalised from demo/19) ──────────────────────────────
## Builds a deterministic multi-period LP with realised noise.
## Fixes applied: alpha interleaving (nJ>=2), terminal holding zeroed, s0=0.

solve_clairvoyant <- function(e, noise_mat, trial_idx) {
  rows  <- (trial_idx - 1L) * e$tau + seq_len(e$tau)
  Q_vec <- c(t(noise_mat[rows, seq(e$nI)]))
  D_vec <- c(t(noise_mat[rows, e$nI + seq(e$nJ)]))

  # Interleave exit coefficients to match LP variable layout
  alpha_il <- c(
    e$alpha[seq(e$nI)],
    c(rbind(e$alpha[e$nI + seq(e$nJ)], e$alpha[e$nI + e$nJ + seq(e$nJ)]))
  )

  alpha_orig <- e$alpha
  e$alpha <- alpha_il

  ccx  <- carrier_capacity_padded(e)
  tlx  <- transition_logic(e, q = Q_vec[seq(e$nI)], d = D_vec[seq(e$nJ)])
  slx  <- storage_limits(e,  q = Q_vec[seq(e$nI)])
  obj_ <- c(alpha_il, e$CTb, e$CTo[1L, ], alpha_il)

  A   <- rbind(ccx$A, tlx$A, slx$A)
  rhs <- c(ccx$rhs, tlx$rhs, slx$rhs)
  sns <- c(ccx$sense, tlx$sense, slx$sense)

  mp <- multiperiod_expansion(e, Q_vec, D_vec, A, obj_, rhs, sns)
  e$alpha <- alpha_orig

  mp$modelsense <- "min"
  mp$vtype      <- rep("C", ncol(mp$A))

  # Zero out terminal holding costs (SDDP charges tau, LP charges tau+1)
  s0_off <- e$nI + 2L * e$nJ
  n_col  <- ncol(mp$A)
  mp$obj[seq(n_col - s0_off + 1L, n_col)] <- 0

  # Pin s0 = 0
  A_s0   <- cbind(Matrix::Diagonal(s0_off),
                   Matrix::Matrix(0L, nrow = s0_off, ncol = n_col - s0_off))
  mp$A     <- rbind(mp$A, A_s0)
  mp$rhs   <- c(mp$rhs,   numeric(s0_off))
  mp$sense <- c(mp$sense, rep("=", s0_off))

  opt <- solve_lp(mp)
  if (is.null(opt$objval) || !is.finite(opt$objval)) return(NA_real_)
  opt$objval
}

## ── Main loop: one topology at a time ───────────────────────────────────────

cat("\n", strrep("=", 70), "\n", sep = "")
cat("21_multi_topology_regret.R — Multi-topology SDDP regret (B1)\n")
cat(strrep("=", 70), "\n\n")

results <- vector("list", length(CONFIGS))

for (ci in seq_along(CONFIGS)) {
  cfg <- CONFIGS[[ci]]
  cat(sprintf("\n%s\n[%s]  nI=%d nJ=%d rate=%.1f\n%s\n",
              strrep("=", 60), cfg$label, cfg$nI, cfg$nJ, cfg$rate,
              strrep("=", 60)))

  ## Generate TLPR instance
  tmpf <- tempfile(fileext = ".json")
  TLPR::generate_instance(
    nI = cfg$nI, nJ = cfg$nJ, tau = 4L,
    nB = 5L, nCS = 10L, nCO = 1L, rate = cfg$rate,
    seed = SEED, path = tmpf
  )

  e <- new.env()
  jsonlite::fromJSON(tmpf) |> list2env(envir = e)
  e$from_i <- as_int_list(e$from_i)
  e$to_j   <- as_int_list(e$to_j)
  create_model(e)
  dp_config(e)

  cat(sprintf("  |S|=%s  nCS=%d  nCO=%d  tau=%d  R=%d\n",
              format(e$nSdx, big.mark = ","), e$nCS, e$nCO, e$tau, e$R))

  ## Build MSTP config with Gaussian copula
  lambda_val <- sum(e$Q$vals * e$Q$prob)
  lambda_vec <- rep(lambda_val, cfg$nI + cfg$nJ)
  nOD        <- cfg$nI + cfg$nJ
  corrmat    <- MSTP::gen_corrmat(n_blocks = 2L,
                                   block_size = max(cfg$nI, cfg$nJ),
                                   cross_corr = 0.0)
  if (nrow(corrmat) > nOD) corrmat <- corrmat[1:nOD, 1:nOD]

  inst   <- tlpr_env_to_mstp_inst(e)
  config <- MSTP::mstp_config(inst, lambda = lambda_vec, corrmat = corrmat,
                               n_scenarios = N_SCENARIOS)

  ## Train SDDP
  cat(sprintf("  Training SDDP (%d iterations)...\n", SDDP_ITER))
  t0    <- proc.time()[["elapsed"]]
  model <- MSTP::mstp_train(config, iterations = SDDP_ITER)
  t_train <- proc.time()[["elapsed"]] - t0

  lb <- MSTP::mstp_bound(model)
  cat(sprintf("  Train: %.1fs  |  LB: %.2f\n", t_train, lb))

  ## Simulate OOB
  cat(sprintf("  Simulating %d OOB trials...\n", N_TRIALS))
  t0   <- proc.time()[["elapsed"]]
  sims <- MSTP::mstp_simulate(model, config, trials = N_TRIALS)
  t_sim <- proc.time()[["elapsed"]] - t0

  ub_mean <- mean(sims$obj)
  ub_sd   <- sd(sims$obj)
  gap_pct <- 100 * (ub_mean - lb) / (abs(lb) + 1e-8)
  cat(sprintf("  Sim:   %.1fs  |  UB: %.2f +/- %.2f  |  Gap: %.2f%%\n",
              t_sim, ub_mean, ub_sd, gap_pct))

  ## Compute per-trial regret
  cat(sprintf("  Computing regret (%d clairvoyant LPs)... ", N_TRIALS))
  flush.console()
  t0_reg <- proc.time()[["elapsed"]]
  regret <- numeric(N_TRIALS)

  for (i in seq_len(N_TRIALS)) {
    lp_opt    <- solve_clairvoyant(e, sims$noise, i)
    regret[i] <- (sims$obj[i] - lp_opt) / (abs(lp_opt) + 1e-8)
  }
  t_reg <- proc.time()[["elapsed"]] - t0_reg
  cat(sprintf("done (%.1fs)\n", t_reg))

  ok <- is.finite(regret)
  q  <- quantile(regret[ok] * 100, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

  results[[ci]] <- list(
    label    = cfg$label,
    nI       = cfg$nI,
    nJ       = cfg$nJ,
    nSdx     = e$nSdx,
    t_train  = t_train,
    t_sim    = t_sim,
    t_regret = t_reg,
    lb       = lb,
    ub_mean  = ub_mean,
    ub_sd    = ub_sd,
    gap_pct  = gap_pct,
    regret   = regret[ok],
    n_ok     = sum(ok),
    quantiles = q
  )

  cat(sprintf("  Regret: mean=%.3f%%  sd=%.3f%%  median=%.3f%%  (n=%d/%d)\n",
              100 * mean(regret[ok]), 100 * sd(regret[ok]),
              100 * median(regret[ok]), sum(ok), N_TRIALS))
  cat(sprintf("  Quantiles: min=%.3f%%  Q5=%.3f%%  Q95=%.3f%%  max=%.3f%%\n",
              q[1], q[2], q[6], q[7]))

  unlink(tmpf)
  gc(verbose = FALSE)
}

## ── Summary table ────────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 90), "\n", sep = "")
cat("MULTI-TOPOLOGY SDDP REGRET\n")
cat(strrep("=", 90), "\n\n")

hdr <- sprintf("%-12s  %8s  %7s  %10s  %10s  %10s  %10s  %10s",
               "Topology", "|S|", "Gap%", "Regret%", "Regret sd", "Median%",
               "Q5%", "Q95%")
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")

for (res in results) {
  cat(sprintf("%-12s  %8s  %6.2f%%  %9.3f%%  %9.3f%%  %9.3f%%  %9.3f%%  %9.3f%%\n",
              res$label,
              format(res$nSdx, big.mark = ","),
              res$gap_pct,
              100 * mean(res$regret),
              100 * sd(res$regret),
              100 * median(res$regret),
              res$quantiles[2],
              res$quantiles[6]))
}

cat(strrep("-", nchar(hdr)), "\n")
cat(sprintf("SDDP iters=%d  |  Trials=%d  |  Scenarios/iter=%d  |  Seed=%d\n",
            SDDP_ITER, N_TRIALS, N_SCENARIOS, SEED))
cat("Copula: Gaussian, 2-block corrmat (cross=0.4), Poisson(lambda) marginals\n")
cat("Regret = (SDDP_cost - LP*) / |LP*|  where LP* = clairvoyant optimal\n")

## ── Regret distribution quantiles (for Q-Q plots) ──────────────────────────

cat("\n\n", strrep("=", 90), "\n", sep = "")
cat("REGRET Q-Q QUANTILES (per topology)\n")
cat(strrep("=", 90), "\n")

qq_probs <- seq(0, 1, by = 0.05)

cat(sprintf("\n%-12s", "Quantile"))
for (res in results) cat(sprintf("  %12s", res$label))
cat("\n")
cat(strrep("-", 12 + length(results) * 14), "\n")

for (p in qq_probs) {
  cat(sprintf("%-12s", sprintf("%.0f%%", p * 100)))
  for (res in results) {
    q <- quantile(res$regret * 100, probs = p)
    cat(sprintf("  %12.3f", q))
  }
  cat("\n")
}

## ── Timing summary ──────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("TIMING SUMMARY\n")
cat(strrep("-", 70), "\n")
cat(sprintf("%-12s  %8s  %8s  %8s  %8s\n",
            "Topology", "Train", "Sim", "Regret", "Total"))
cat(strrep("-", 70), "\n")

for (res in results) {
  total <- res$t_train + res$t_sim + res$t_regret
  cat(sprintf("%-12s  %7.1fs  %7.1fs  %7.1fs  %7.1fs\n",
              res$label, res$t_train, res$t_sim, res$t_regret, total))
}
cat(strrep("-", 70), "\n")

## ── LaTeX table snippet ─────────────────────────────────────────────────────

cat("\n%% LaTeX snippet for multi-topology regret table (§6)\n")
cat("%% Topology & |S| & LB & UB +/- sd & Gap% & Regret% & Regret sd%\n")
for (res in results) {
  topo <- sub("x", "$\\\\times$", sub(" R=(\\d+)", " ($R=\\1$)", res$label))
  cat(sprintf(
    "%s & %s & %.2f & $%.2f \\pm %.0f$ & %.2f\\%% & %.3f\\%% & %.3f\\%% \\\\\n",
    topo,
    format(res$nSdx, big.mark = ","),
    res$lb, res$ub_mean, res$ub_sd, res$gap_pct,
    100 * mean(res$regret), 100 * sd(res$regret)))
}

cat("\n%% Q-Q quantiles (pgfplots / tikz)\n")
cat("%% prob")
for (res in results) cat(sprintf("  %s", res$label))
cat("\n")
for (p in qq_probs) {
  cat(sprintf("%% %.2f", p))
  for (res in results) {
    q <- quantile(res$regret * 100, probs = p)
    cat(sprintf("  %.4f", q))
  }
  cat("\n")
}

## ── Persist results ───────────────────────────────────────────────────────────
dir.create("demo/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, "demo/results/21_results.rds")
cat("\nResults saved to demo/results/21_results.rds\n")

invisible(gc(verbose = FALSE))
