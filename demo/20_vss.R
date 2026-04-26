## ============================================================
## 20_vss.R — Value of Stochastic Solution (A4)
##
## Computes gain of recourse: how much better the multi-stage
## SDDP policy is compared to a myopic (greedy single-stage)
## policy on the same OOB noise trajectories.
##
##   gain_i = (myopic_cost_i - SDDP_cost_i) / myopic_cost_i
##
## The myopic policy solves a single-stage LP at each period
## with one-step look-ahead costs (holding + transport), but
## no future cost-to-go.  It uses the same cost coefficients
## as the SDDP model.
##
## Instance: 2x2, R=10, tau=4  (same as demo/18, demo/19)
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(MSTP)
})

SEED        <- 42L
SDDP_ITER   <- 1000L
N_TRIALS    <- 500L
N_SCENARIOS <- 10L

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

## ── Generate TLPR instance ──────────────────────────────────────────────────

cat("\n", strrep("=", 70), "\n", sep = "")
cat("20_vss.R — Value of Stochastic Solution (A4)\n")
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

## ── Adapter (with CTb fix) ──────────────────────────────────────────────────

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

cat(sprintf("Lambda = %.3f x %d dims  |  Corr: 2-block (cross=0.0, O-D independent)\n",
            lambda_val, nOD))

## ── Train + simulate ────────────────────────────────────────────────────────

cat(sprintf("\nTraining SDDP (%d iterations)...\n", SDDP_ITER))
model <- MSTP::mstp_train(config, iterations = SDDP_ITER)
lb    <- MSTP::mstp_bound(model)

cat(sprintf("Simulating %d OOB trials...\n", N_TRIALS))
sims <- MSTP::mstp_simulate(model, config, trials = N_TRIALS)

cat(sprintf("LB = %.2f  |  UB (sim mean) = %.2f +/- %.2f\n\n",
            lb, mean(sims$obj), sd(sims$obj)))

## ── Cost coefficients (match SDDP model) ────────────────────────────────────

entry_store_coef <- e$alpha[seq(e$nI)]
exit_store_coef  <- e$alpha[e$nI + seq(e$nJ)]
exit_short_coef  <- e$alpha[e$nI + e$nJ + seq(e$nJ)]
transport_coef   <- e$CTb
nLanes           <- e$nL_
nSpotLanes       <- e$nI * e$nJ * e$nCO

## ── Myopic single-stage LP ──────────────────────────────────────────────────
## Solves a one-period LP at each stage.  Decision: carrier allocations.
## Objective: current-period holding costs (constant) + transport costs
##            + next-period holding/shortage costs on resulting state.
## This gives the greedy policy an incentive to move containers and
## satisfy demand, but with no multi-period look-ahead.

solve_myopic_stage <- function(t, entry_in, exitp_in, exitm_in, Q_t, D_t) {
  nMoves <- e$nvars
  nVars  <- nMoves + e$nI + 2L * e$nJ  # moves + entry_out + exitp_out + exitm_out

  col_eo <- nMoves + seq_len(e$nI)
  col_xp <- nMoves + e$nI + seq_len(e$nJ)
  col_xm <- nMoves + e$nI + e$nJ + seq_len(e$nJ)

  # Spot cost at stage t: spot_coef is as.vector(CTo), column-major tau x nSpotLanes
  spot_t <- numeric(nSpotLanes)
  for (sl in seq_len(nSpotLanes)) {
    spot_t[sl] <- e$CTo[(sl - 1L) * e$tau + t]
  }

  obj <- c(transport_coef, spot_t,       # move costs
           entry_store_coef,             # entry_out holding
           exit_store_coef,              # exitp_out holding
           exit_short_coef)              # exitm_out shortage

  # Constraints: carrier cap, entry transition, exit capacity, exit balance
  nCons <- e$nCS + e$nCO + e$nI + e$nJ + e$nJ
  is_ <- integer(0L); js_ <- integer(0L); xs_ <- numeric(0L)
  add <- function(i, j, x = 1.0) {
    is_ <<- c(is_, i); js_ <<- c(js_, j); xs_ <<- c(xs_, x)
  }

  # Carrier capacity
  for (k in seq_len(e$nCS))
    for (col in (e$nLc[k] + 1L):e$nLc[k + 1L]) add(k, col)
  for (c_idx in seq_len(e$nCO))
    for (col in nLanes + (c_idx - 1L) * (e$nI * e$nJ) + seq_len(e$nI * e$nJ))
      add(e$nCS + c_idx, col)

  # Entry transition: entry_out[i] + sum(move from i) = entry_in[i] + Q_t[i]
  r0 <- e$nCS + e$nCO
  for (i in seq_len(e$nI)) {
    add(r0 + i, col_eo[i])
    for (col in e$from_i[[i]]) add(r0 + i, col)
  }

  # Exit capacity: sum(move to j) <= exit_cap[j] - exitp_in[j]
  r0 <- e$nCS + e$nCO + e$nI
  for (j in seq_len(e$nJ))
    for (col in e$to_j[[j]]) add(r0 + j, col)

  # Exit balance: exitp_out[j] - exitm_out[j] - sum(move to j) = exitp_in - exitm_in - D_t[j]
  r0 <- e$nCS + e$nCO + e$nI + e$nJ
  for (j in seq_len(e$nJ)) {
    add(r0 + j, col_xp[j],  1.0)
    add(r0 + j, col_xm[j], -1.0)
    for (col in e$to_j[[j]]) add(r0 + j, col, -1.0)
  }

  A <- methods::as(
    Matrix::sparseMatrix(i = is_, j = js_, x = xs_, dims = c(nCons, nVars)),
    "CsparseMatrix"
  )

  sense <- c(rep("<", e$nCS + e$nCO),
             rep("=", e$nI),
             rep("<", e$nJ),
             rep("=", e$nJ))

  rhs <- c(e$Cb[t, seq_len(e$nCS)],
           e$Co[t, seq_len(e$nCO)],
           entry_in + Q_t,
           e$R - exitp_in,
           exitp_in - exitm_in - D_t)

  res <- solve_lp(list(obj = obj, A = A, rhs = rhs,
                       sense = sense, modelsense = "min",
                       vtype = rep("C", nVars)))
  x <- pmax(0, res$x)

  # Current-period holding (constant term, not in LP objective)
  state_cost <- sum(entry_store_coef * entry_in) +
                sum(exit_store_coef  * exitp_in) +
                sum(exit_short_coef  * exitm_in)

  list(cost      = state_cost + res$objval,
       entry_out = x[col_eo],
       exitp_out = x[col_xp],
       exitm_out = x[col_xm])
}

## ── Compute VSS for each trial ──────────────────────────────────────────────

cat(sprintf("Computing myopic policy costs (%d trials x %d stages)...\n",
            N_TRIALS, e$tau))
flush.console()

myopic_cost <- numeric(N_TRIALS)
t0 <- proc.time()[["elapsed"]]

for (i in seq_len(N_TRIALS)) {
  entry_in <- rep(0, e$nI)
  exitp_in <- rep(0, e$nJ)
  exitm_in <- rep(0, e$nJ)
  total    <- 0

  for (t in seq_len(e$tau)) {
    row <- (i - 1L) * e$tau + t
    Q_t <- sims$noise[row, seq(e$nI)]
    D_t <- sims$noise[row, e$nI + seq(e$nJ)]

    s <- solve_myopic_stage(t, entry_in, exitp_in, exitm_in, Q_t, D_t)
    total    <- total + s$cost
    entry_in <- s$entry_out
    exitp_in <- s$exitp_out
    exitm_in <- s$exitm_out
  }
  myopic_cost[i] <- total
}
t_myopic <- proc.time()[["elapsed"]] - t0

gain <- (myopic_cost - sims$obj) / myopic_cost

cat(sprintf("Done (%.1fs)\n\n", t_myopic))

## ── Summary ─────────────────────────────────────────────────────────────────

cat(strrep("=", 70), "\n")
cat("VALUE OF STOCHASTIC SOLUTION (GAIN OF RECOURSE)\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("SDDP cost:   mean = %.2f  sd = %.2f\n", mean(sims$obj), sd(sims$obj)))
cat(sprintf("Myopic cost: mean = %.2f  sd = %.2f\n", mean(myopic_cost), sd(myopic_cost)))
cat(sprintf("\nGain of recourse = (myopic - SDDP) / myopic\n"))
cat(sprintf("  Mean:   %.2f%%\n", 100 * mean(gain)))
cat(sprintf("  Median: %.2f%%\n", 100 * median(gain)))
cat(sprintf("  SD:     %.2f%%\n", 100 * sd(gain)))

q <- quantile(gain * 100, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
cat(sprintf("\nQuantiles (%%): min=%.2f  Q5=%.2f  Q1=%.2f  med=%.2f  Q3=%.2f  Q95=%.2f  max=%.2f\n",
            q[1], q[2], q[3], q[4], q[5], q[6], q[7]))

cat(sprintf("\nTrials: %d  |  SDDP iters: %d  |  Instance: 2x2 R=10 tau=4 seed=%d\n",
            N_TRIALS, SDDP_ITER, SEED))
cat("Myopic = greedy single-stage LP with one-step holding look-ahead\n")
cat("SDDP = multi-stage policy with Benders cuts (copula, single-cut)\n")

## ── LaTeX snippet ───────────────────────────────────────────────────────────

cat(sprintf("\n%% LaTeX snippet for VSS / gain of recourse\n"))
cat(sprintf("%% SDDP mean & Myopic mean & Gain mean & Gain sd\n"))
cat(sprintf("$%.2f$ & $%.2f$ & $%.2f\\%%$ & $%.2f\\%%$ \\\\\n",
            mean(sims$obj), mean(myopic_cost),
            100 * mean(gain), 100 * sd(gain)))

unlink(tmpf)

## ── Persist results ───────────────────────────────────────────────────────────
dir.create("demo/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(sims_obj = sims$obj, myopic_cost = myopic_cost, gain = gain, lb = lb),
        "demo/results/20_results.rds")
cat("\nResults saved to demo/results/20_results.rds\n")

invisible(gc(verbose = FALSE))
