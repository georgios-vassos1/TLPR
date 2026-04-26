## ============================================================
## 19_regret_stopping.R — Regret-based stopping criterion (B5)
##
## Trains SDDP at multiple iteration counts on the same instance,
## simulates OOB trials at each, and computes per-trial regret:
##
##   regret_i = (SDDP_cost_i - LP*_i) / |LP*_i|
##
## where LP*_i is the clairvoyant (perfect-foresight) LP optimal
## for trial i's realised noise.
##
## Goal: show that regret stabilises early under Gaussian copula +
## single-cut SDDP, contradicting Schmiedel2025's finding that
## regret-based stopping "should not be recommended."
##
## Instance: 2×2, R=10, tau=4  (same as demo/18)
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(MSTP)
})

## ── Parameters ────────────────────────────────────────────────────────────────

ITER_COUNTS <- c(100L, 250L, 500L, 1000L, 1500L)
N_TRIALS    <- 500L       # OOB simulation trials per iteration count
N_SCENARIOS <- 10L        # training scenarios per SDDP iteration
SEED        <- 42L

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

## ── Generate TLPR instance ───────────────────────────────────────────────────

cat("\n", strrep("=", 70), "\n", sep = "")
cat("19_regret_stopping.R — Regret vs SDDP iterations (B5)\n")
cat(strrep("=", 70), "\n\n")

tmpf <- tempfile(fileext = ".json")
TLPR::generate_instance(
  nI = 2L, nJ = 2L, tau = 4L,
  nB = 5L, nCS = 10L, nCO = 1L, rate = 1.0,
  seed = SEED, path = tmpf
)

e <- new.env()
jsonlite::fromJSON(tmpf) |> list2env(envir = e)
e$from_i <- as_int_list(e$from_i)
e$to_j   <- as_int_list(e$to_j)
create_model(e)
dp_config(e)

cat(sprintf("Instance: nI=%d  nJ=%d  R=%d  tau=%d  nCS=%d  nCO=%d\n",
            e$nI, e$nJ, e$R, e$tau, e$nCS, e$nCO))

## ── Build MSTP config (reuse adapter from demo/18) ───────────────────────────

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
    transport_coef   = e$CTb,  # carrier-index order, matches move variable layout
    spot_coef        = as.vector(e$CTo),
    entry_stock_0    = rep(0L, e$nI),
    exit_stock_0     = rep(0L, e$nJ),
    exit_short_0     = rep(0L, e$nJ)
  )
}

lambda_val <- sum(e$Q$vals * e$Q$prob)
lambda_vec <- rep(lambda_val, e$nI + e$nJ)
nOD        <- e$nI + e$nJ
corrmat    <- MSTP::gen_corrmat(n_blocks = 2L, block_size = max(e$nI, e$nJ),
                                 cross_corr = 0.0)
if (nrow(corrmat) > nOD) corrmat <- corrmat[1:nOD, 1:nOD]

inst   <- tlpr_env_to_mstp_inst(e)
config <- MSTP::mstp_config(inst, lambda = lambda_vec, corrmat = corrmat,
                             n_scenarios = N_SCENARIOS)

cat(sprintf("Lambda = %.3f × %d dims  |  Corr: 2-block (cross=0.0, O-D independent)\n",
            lambda_val, nOD))
cat(sprintf("Scenarios/iter = %d  |  OOB trials = %d\n", N_SCENARIOS, N_TRIALS))

## ── Clairvoyant LP builder ───────────────────────────────────────────────────
## For each trial, builds a deterministic multi-period LP with the realised
## noise and solves to get the perfect-foresight optimal cost.
## Uses TLPR's constraint builders (same pattern as demo/14).

solve_clairvoyant <- function(noise_mat, trial_idx) {
  rows  <- (trial_idx - 1L) * e$tau + seq_len(e$tau)
  Q_vec <- c(t(noise_mat[rows, seq(e$nI)]))           # length tau*nI
  D_vec <- c(t(noise_mat[rows, e$nI + seq(e$nJ)]))    # length tau*nJ

  # Interleave exit coefficients to match LP variable layout
  # LP vars: [entry(nI), exitp_1, exitm_1, exitp_2, exitm_2, ...]
  # alpha:   [entry_store(nI), exit_store(nJ), exit_short(nJ)]
  alpha_il <- c(
    e$alpha[seq(e$nI)],
    c(rbind(e$alpha[e$nI + seq(e$nJ)], e$alpha[e$nI + e$nJ + seq(e$nJ)]))
  )

  # Temporarily swap alpha for multiperiod_expansion (uses env$alpha internally)
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
  e$alpha <- alpha_orig   # restore

  mp$modelsense <- "min"
  mp$vtype      <- rep("C", ncol(mp$A))

  # Zero out terminal holding costs — SDDP charges holding on .in only
  # (tau charges total), but multiperiod_expansion adds an extra alpha
  # on state_end of the last period (tau+1 charges). Remove it to match.
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

## ── Main loop: train → simulate → regret at each iteration count ─────────────

results <- vector("list", length(ITER_COUNTS))

for (ki in seq_along(ITER_COUNTS)) {
  iters <- ITER_COUNTS[ki]
  cat(sprintf("\n--- SDDP iterations = %d ---\n", iters))

  t0    <- proc.time()[["elapsed"]]
  model <- MSTP::mstp_train(config, iterations = iters)
  t_train <- proc.time()[["elapsed"]] - t0

  lb <- MSTP::mstp_bound(model)
  cat(sprintf("  Train: %.1fs  |  LB: %.2f\n", t_train, lb))

  t0   <- proc.time()[["elapsed"]]
  sims <- MSTP::mstp_simulate(model, config, trials = N_TRIALS)
  t_sim <- proc.time()[["elapsed"]] - t0

  ub_mean <- mean(sims$obj)
  ub_sd   <- sd(sims$obj)
  gap_pct <- 100 * (ub_mean - lb) / (abs(lb) + 1e-8)
  cat(sprintf("  Sim:   %.1fs  |  UB: %.2f ± %.2f  |  Gap: %.2f%%\n",
              t_sim, ub_mean, ub_sd, gap_pct))

  # Compute per-trial regret
  cat(sprintf("  Computing regret (%d clairvoyant LPs)... ", N_TRIALS))
  flush.console()
  t0_reg <- proc.time()[["elapsed"]]
  regret <- numeric(N_TRIALS)

  for (i in seq_len(N_TRIALS)) {
    lp_opt    <- solve_clairvoyant(sims$noise, i)
    regret[i] <- (sims$obj[i] - lp_opt) / (abs(lp_opt) + 1e-8)
  }
  t_reg <- proc.time()[["elapsed"]] - t0_reg
  cat(sprintf("done (%.1fs)\n", t_reg))

  ok <- is.finite(regret)
  results[[ki]] <- list(
    iters    = iters,
    t_train  = t_train,
    t_sim    = t_sim,
    t_regret = t_reg,
    lb       = lb,
    ub_mean  = ub_mean,
    ub_sd    = ub_sd,
    gap_pct  = gap_pct,
    regret   = regret[ok],
    n_ok     = sum(ok)
  )

  cat(sprintf("  Regret: mean=%.4f%%  median=%.4f%%  sd=%.4f%%  (n=%d/%d)\n",
              100 * mean(regret[ok]), 100 * median(regret[ok]),
              100 * sd(regret[ok]), sum(ok), N_TRIALS))
}

## ── Summary table ─────────────────────────────────────────────────────────────

cat("\n", strrep("=", 90), "\n", sep = "")
cat("REGRET STABILISATION ACROSS SDDP ITERATION COUNTS\n")
cat(strrep("=", 90), "\n\n")

hdr <- sprintf("%-8s  %8s  %10s  %12s  %10s  %12s  %12s",
               "Iters", "t_train", "LB", "UB ± sd", "Gap%",
               "Regret%", "Regret sd%")
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")

for (res in results) {
  cat(sprintf("%-8d  %8.1f  %10.2f  %7.2f ± %3.0f  %9.2f%%  %11.4f%%  %11.4f%%\n",
              res$iters, res$t_train, res$lb,
              res$ub_mean, res$ub_sd, res$gap_pct,
              100 * mean(res$regret), 100 * sd(res$regret)))
}

cat(strrep("-", nchar(hdr)), "\n")
cat(sprintf("Instance: 2×2 R=10 tau=4 seed=%d  |  Trials=%d  |  Scenarios/iter=%d\n",
            SEED, N_TRIALS, N_SCENARIOS))
cat("Copula: Gaussian, 2-block corrmat (cross=0.4), Poisson(lambda) marginals\n")
cat("Regret = (SDDP_cost - LP*) / |LP*|  where LP* = clairvoyant optimal\n")

## ── Regret distribution quantiles ────────────────────────────────────────────

cat("\n", strrep("=", 90), "\n", sep = "")
cat("Regret distribution quantiles (%):\n")
cat(strrep("=", 90), "\n")

cat(sprintf("%-8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
            "Iters", "min", "Q1", "median", "Q3", "Q95", "max"))
cat(strrep("-", 60), "\n")

for (res in results) {
  q <- quantile(res$regret * 100, probs = c(0, 0.25, 0.5, 0.75, 0.95, 1))
  cat(sprintf("%-8d  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
              res$iters, q[1], q[2], q[3], q[4], q[5], q[6]))
}

cat(strrep("-", 60), "\n")

## ── LaTeX table snippet ──────────────────────────────────────────────────────

cat("\n% LaTeX snippet for §5.3 / §6\n")
cat("% Iters & t_train & LB & UB±sd & Gap% & Regret% & Regret sd%\n")
for (res in results) {
  cat(sprintf(
    "%d & %.1f & %.2f & $%.2f \\pm %.0f$ & %.2f\\%% & %.4f\\%% & %.4f\\%% \\\\\n",
    res$iters, res$t_train, res$lb,
    res$ub_mean, res$ub_sd, res$gap_pct,
    100 * mean(res$regret), 100 * sd(res$regret)))
}

cat("\nKey finding: if regret mean and sd are stable across iteration counts,\n")
cat("regret-based stopping is viable — contradicting Schmiedel2025 §5.3.\n")
cat("Their regular SDDP showed regret up to 1100%; their slow+forget variant\n")
cat("achieved 0.3-3.5% but dismissed regret stopping due to large negative biases.\n")
cat("Under copula uncertainty + single-cut, we expect well-concentrated regret.\n")

unlink(tmpf)

## ── Persist results ───────────────────────────────────────────────────────────
dir.create("demo/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, "demo/results/19_results.rds")
cat("\nResults saved to demo/results/19_results.rds\n")

invisible(gc(verbose = FALSE))
