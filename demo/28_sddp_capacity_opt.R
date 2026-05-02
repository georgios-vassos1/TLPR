## ============================================================
## 28_sddp_capacity_opt.R — SDDP capacity optimisation, 6×6×20 (P1.5)
##
## !! SUPERSEDED — NOT USED IN THE PAPER !!
##
## Superseded by demo/31_vocc_capacity.R, which extends this experiment
## with a VoCC comparison (correct copula vs independence) on the same
## MSTP instance and calibration (λ=50).
##
## Part A (LP proxy vs SDDP duals, 2×2 TLPR instance) has been removed.
## It was superseded by demo/28b (itself archived at demo/archived/).
## The relevant scientific content is now in demo/31.
##
## ── What this file contains ─────────────────────────────────────────
## SDDP projected-gradient capacity optimisation on a 6×6×20 MSTP
## instance (τ=12, λ=50). Demonstrates that SDDP capacity duals guide
## meaningful capacity reallocation, achieving ~8% total cost reduction
## over 15 outer iterations. This is the direct predecessor of demo/31.
## ============================================================

suppressPackageStartupMessages({
  library(MSTP)
})

## ── Parameters ───────────────────────────────────────────────────────────────

SEED        <- 42L
N_TRIALS    <- 300L
N_SCENARIOS <- 10L
LAMBDA_VAL  <- 50.0   # reduced from 700: CV=14% vs 3.8%, P(capacity binds)~17% vs ~0%
CROSS_CORR  <- 0.0

N_CAP_ITER     <- 15L
N_CAP_SAMPLES  <- 200L
SDDP_ITER_COLD <- 100L
STEP_SIZE      <- 20.0
CAP_FLOOR      <- 1.0

## ── Header ───────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 70), "\n", sep = "")
cat("28_sddp_capacity_opt.R — SDDP capacity optimisation (6×6×20, τ=12)\n")
cat(strrep("=", 70), "\n\n")
cat("NOTE: superseded by demo/31_vocc_capacity.R\n\n")

## ── Instance ─────────────────────────────────────────────────────────────────

cat("PART B — SDDP dual-based capacity optimisation (6×6×20, τ=12)\n\n")

inst <- MSTP::generate_instance(tau=12L, nOrigins=6L, nDestinations=6L,
                                 nCarriers=20L, seed=SEED, lambda=LAMBDA_VAL)
nI   <- inst$nOrigins
nJ   <- inst$nDestinations
nCS  <- inst$nCarriers
nCSO <- inst$nSpotCarriers
tau  <- inst$tau
nOD  <- nI + nJ

cat(sprintf("Instance: generated (6×6×20, τ=%d, seed=%d)\n", tau, SEED))
cat(sprintf("  nI=%d  nJ=%d  nCS=%d  nCSO=%d  τ=%d\n", nI, nJ, nCS, nCSO, tau))

lambda_vec <- rep(LAMBDA_VAL, nOD)
cat(sprintf("  lambda=%.0f (reduced from 700; higher CV makes capacity binding)\n",
            LAMBDA_VAL))
set.seed(SEED)
corrmat <- MSTP::gen_corrmat(n_blocks   = 2L,
                              block_size = max(nI, nJ),
                              cross_corr = CROSS_CORR)

## Reservation cost: 0.2 × mean transport rate per capacity unit per period.
v_b <- rep(mean(inst$transport_coef) * 0.2, (nCS + nCSO) * tau)

x0_b    <- as.numeric(inst$carrier_capacity)
cap_max <- max(x0_b) * 2.0

## ── Gradient step helper ──────────────────────────────────────────────────────

train_step <- function(x, n_iter = SDDP_ITER_COLD) {
  inst_x  <- MSTP::mstp_update_capacity(inst, x)
  config  <- MSTP::mstp_config(inst_x,
                                lambda      = lambda_vec,
                                corrmat     = corrmat,
                                n_scenarios = N_SCENARIOS)
  model   <- MSTP::mstp_train(config, n_iter)
  lb    <- MSTP::mstp_bound(model)
  duals <- MSTP::mstp_capacity_duals(model, config, n_samples = N_CAP_SAMPLES)
  grad  <- duals + v_b
  list(model = model, config = config, lb = lb,
       obj   = lb + sum(v_b * x),
       grad  = grad, grad_norm = sqrt(sum(grad^2)))
}

conv_hist <- vector("list", N_CAP_ITER + 1L)

## ── Step 1: Cold start at x0 ─────────────────────────────────────────────────

cat(sprintf("Step 1: SDDP at x0 (%d iters)...\n", SDDP_ITER_COLD))
t0_b  <- proc.time()[["elapsed"]]
step  <- train_step(x0_b, n_iter = SDDP_ITER_COLD)
t_s   <- proc.time()[["elapsed"]] - t0_b
conv_hist[[1]] <- list(iter=0L, obj=step$obj, lb=step$lb,
                        grad_norm=step$grad_norm, x_change=NA_real_,
                        t=t_s, warm=FALSE, n_iter=SDDP_ITER_COLD)
cat(sprintf("  iter  0 | LB=%11.0f | obj=%13.0f | |∇|=%8.3f | Δt=%.1fs [cold]\n",
            step$lb, step$obj, step$grad_norm, t_s))
x <- x0_b

## ── Step 2: Projected gradient ───────────────────────────────────────────────

cat(sprintf("\nStep 2: projected gradient (%d steps, normalised step=%.2f/sqrt(k), cold SDDP=%d iters)\n",
            N_CAP_ITER, STEP_SIZE, SDDP_ITER_COLD))

t0_opt <- proc.time()[["elapsed"]]

for (k in seq_len(N_CAP_ITER)) {
  alpha_k  <- STEP_SIZE / sqrt(k)
  g_hat    <- step$grad / (step$grad_norm + 1e-12)
  x_new    <- round(pmax(CAP_FLOOR, pmin(cap_max, x - alpha_k * g_hat)))
  x_change <- sqrt(sum((x_new - x)^2))

  t0_s  <- proc.time()[["elapsed"]]
  step  <- train_step(x_new, n_iter = SDDP_ITER_COLD)
  t_s   <- proc.time()[["elapsed"]] - t0_s

  conv_hist[[k + 1L]] <- list(iter=k, obj=step$obj, lb=step$lb,
                               grad_norm=step$grad_norm, x_change=x_change,
                               t=t_s, warm=FALSE, n_iter=SDDP_ITER_COLD)
  cat(sprintf("  iter %2d | LB=%11.0f | obj=%13.0f | |∇|=%8.3f | |Δx|=%7.1f | Δt=%.1fs [cold]\n",
              k, step$lb, step$obj, step$grad_norm, x_change, t_s))

  x <- x_new
  step$model <- NULL; gc(verbose = FALSE)
}

t_opt    <- proc.time()[["elapsed"]] - t0_opt
x_star_b <- x

## ── Step 3: Baseline UB at x0 ────────────────────────────────────────────────

cat(sprintf("\nStep 3: full simulation at x0 (%d iters × %d trials)...\n",
            SDDP_ITER_COLD, N_TRIALS))
t0_b    <- proc.time()[["elapsed"]]
config0 <- MSTP::mstp_config(inst, lambda=lambda_vec, corrmat=corrmat,
                              n_scenarios=N_SCENARIOS)
model0  <- MSTP::mstp_train(config0, iterations = SDDP_ITER_COLD)
lb0     <- MSTP::mstp_bound(model0)
sims0   <- MSTP::mstp_simulate(model0, config0, trials = N_TRIALS)
t_base  <- proc.time()[["elapsed"]] - t0_b
ub0     <- mean(sims0$obj); ub0_sd <- sd(sims0$obj)
ub0_with_v <- ub0 + sum(v_b * x0_b)
cat(sprintf("  LB=%.0f  UB=%.0f±%.0f  UB+v·x0=%.0f  (%.1fs)\n",
            lb0, ub0, ub0_sd, ub0_with_v, t_base))
rm(sims0, model0); gc(verbose = FALSE)

## ── Step 4: Validate at x* ───────────────────────────────────────────────────

cat(sprintf("\nStep 4: SDDP at x* (%d iters × %d trials)...\n",
            SDDP_ITER_COLD, N_TRIALS))
last_step  <- conv_hist[[N_CAP_ITER + 1L]]

tryCatch({
  t0_b      <- proc.time()[["elapsed"]]
  inst_star <- MSTP::mstp_update_capacity(inst, x_star_b)
  config_s  <- MSTP::mstp_config(inst_star, lambda=lambda_vec, corrmat=corrmat,
                                  n_scenarios=N_SCENARIOS)
  model_s   <- MSTP::mstp_train(config_s, SDDP_ITER_COLD)
  lb_star   <- MSTP::mstp_bound(model_s)
  sims_s    <- MSTP::mstp_simulate(model_s, config_s, trials = N_TRIALS)
  t_star    <- proc.time()[["elapsed"]] - t0_b
  ub_star   <- mean(sims_s$obj); ub_star_sd <- sd(sims_s$obj)
  ubstar_with_v <- ub_star + sum(v_b * x_star_b)
  cat(sprintf("  LB=%.0f  UB=%.0f±%.0f  UB+v·x*=%.0f  (%.1fs)\n",
              lb_star, ub_star, ub_star_sd, ubstar_with_v, t_star))
  rm(sims_s, model_s, inst_star); gc(verbose = FALSE)
}, error = function(e) {
  cat(sprintf("  [!] Step 4 evaluation failed: %s\n", conditionMessage(e)))
  cat(sprintf("  Using last gradient-step LB+v·x proxy: %.0f\n", last_step$obj))
  lb_star        <<- last_step$lb
  ub_star        <<- NA_real_
  ub_star_sd     <<- NA_real_
  ubstar_with_v  <<- last_step$obj
})

## ── Summary ───────────────────────────────────────────────────────────────────

sddp_b_reduction <- (1.0 - ubstar_with_v / ub0_with_v) * 100.0

cat("\n\n", strrep("=", 80), "\n", sep = "")
cat("CONVERGENCE HISTORY\n")
cat(strrep("=", 80), "\n")
chdr <- sprintf("%-6s  %-5s  %6s  %13s  %13s  %9s  %8s  %6s",
                "iter", "start", "iters", "LB", "obj (LB+v·x)", "|∇|", "|Δx|", "t(s)")
cat(chdr, "\n")
cat(strrep("-", nchar(chdr)), "\n")
for (h in conv_hist) {
  cat(sprintf("%-6d  %-5s  %6d  %13.0f  %13.0f  %9.3f  %8s  %5.1fs\n",
              h$iter,
              if (h$warm) "warm" else "cold",
              h$n_iter,
              h$lb, h$obj, h$grad_norm,
              if (is.na(h$x_change)) "      —" else sprintf("%8.1f", h$x_change),
              h$t))
}
cat(strrep("-", nchar(chdr)), "\n")

cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("SDDP DUAL-BASED CAPACITY OPTIMISATION (6×6×20)\n")
cat(strrep("=", 70), "\n")
cat(sprintf("%-20s  %14s  %14s  %8s\n",
            "Metric", "x0 (original)", "x* (optimised)", "Change"))
cat(strrep("-", 62), "\n")
cat(sprintf("%-20s  %14.0f  %14.0f  %7.1f%%\n",
            "SDDP LB", lb0, lb_star,
            100*(lb_star - lb0) / (abs(lb0) + 1e-8)))
cat(sprintf("%-20s  %14.0f  %14.0f  %7.1f%%\n",
            "SDDP UB (mean)", ub0, ub_star,
            100*(ub_star - ub0) / (abs(ub0) + 1e-8)))
cat(sprintf("%-20s  %14.0f  %14.0f  %7.1f%%\n",
            "UB + v·x", ub0_with_v, ubstar_with_v, sddp_b_reduction))
cat(strrep("-", 62), "\n")
cat(sprintf("Cold steps: %d × %d iters  |  Dual samples: %d\n",
            N_CAP_ITER + 1L, SDDP_ITER_COLD, N_CAP_SAMPLES))
cat(sprintf("Opt wall time: %.1fs\n", t_opt))

## ── Persist ───────────────────────────────────────────────────────────────────

dir.create("demo/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(
    conv_hist     = conv_hist,
    lb0           = lb0,
    ub0           = ub0,
    ub0_with_v    = ub0_with_v,
    lb_star       = lb_star,
    ub_star       = ub_star,
    ubstar_with_v = ubstar_with_v,
    x0            = x0_b,
    x_star        = x_star_b,
    sddp_reduction = sddp_b_reduction
  ),
  "demo/results/28_results.rds"
)
cat("\nResults saved to demo/results/28_results.rds\n")

invisible(gc(verbose = FALSE))
