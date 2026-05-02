## ============================================================
## 32_lp_proxy_mstp.R — LP proxy vs SDDP duals on MSTP instance (G5)
##
## Head-to-head comparison on the 6×6×20 MSTP instance (τ=12, λ=50):
##
##   x0         — default carrier capacities (calibrated to 80% utilisation)
##   x*_LP      — L-BFGS-B on deterministic mean-scenario LP
##   x*_SDDP    — normalised projected gradient from SDDP capacity duals
##                (same 15-iteration budget as demo/28)
##
## LP proxy: solve the multi-period LP with demands fixed at their Poisson
## means (inflow = outflow = λ = 50 per node).  L-BFGS-B optimises x using
## the exact LP dual w.r.t. capacity constraints.  This recovers a capacity
## vector calibrated for the modal scenario only — ignoring demand variance.
##
## All three are evaluated by fresh SDDP (100 iters × 300 OOB trials).
## ============================================================

suppressPackageStartupMessages({
  library(MSTP)
  library(highs)
  library(Matrix)
})

## ── Parameters ────────────────────────────────────────────────────────────────

SEED          <- 42L
LAMBDA_VAL    <- 50.0
CROSS_CORR    <- 0.0
N_SCENARIOS   <- 10L
SDDP_ITER     <- 100L
N_TRIALS      <- 300L
LP_MAXIT      <- 300L

N_CAP_ITER    <- 15L
N_CAP_SAMPLES <- 200L
STEP_SIZE     <- 20.0
CAP_FLOOR     <- 1.0

## ── Instance ──────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 70), "\n", sep = "")
cat("32_lp_proxy_mstp.R — LP proxy vs SDDP duals (6×6×20, τ=12, λ=50)\n")
cat(strrep("=", 70), "\n\n")

inst <- MSTP::generate_instance(tau=12L, nOrigins=6L, nDestinations=6L,
                                 nCarriers=20L, seed=SEED, lambda=LAMBDA_VAL)

nI   <- inst$nOrigins;   nJ  <- inst$nDestinations
nCS  <- inst$nCarriers;  nCO <- inst$nSpotCarriers
tau  <- inst$tau;        nOD <- nI + nJ

## Derived indices (replicating HyperParams constructor in types.jl)
Ldx       <- inst$Ldx          # OD lane index for each contracted lane position
nLc       <- inst$nLc          # cumulative lane counts: length nCS+1
nLanes    <- nLc[length(nLc)]  # total contracted lanes
nSpotLanes <- nI * nJ * nCO    # total spot lanes

cat(sprintf("Instance: nI=%d  nJ=%d  nCS=%d  nCO=%d  τ=%d\n", nI, nJ, nCS, nCO, tau))
cat(sprintf("          nLanes=%d  nSpotLanes=%d\n", nLanes, nSpotLanes))

## Carrier lane index blocks
CarrierIdx_contr <- lapply(seq_len(nCS), function(k) seq(nLc[k] + 1L, nLc[k + 1L]))
CarrierIdx_spot  <- lapply(seq_len(nCO),
                           function(k) seq((k - 1L) * nI * nJ + 1L, k * nI * nJ))

## Origin / destination lane membership (for contracted and spot lanes)
spotLdx  <- rep(seq_len(nI * nJ), nCO)

origin_lanes <- function(i) sapply(seq_len(nJ), function(j) (j - 1L) * nI + i)
dest_lanes   <- function(j) seq.int((j - 1L) * nI + 1L, j * nI)

from_SG <- lapply(seq_len(nI), function(i) which(Ldx    %in% origin_lanes(i)))
from_SP <- lapply(seq_len(nI), function(i) which(spotLdx %in% origin_lanes(i)) + nLanes)
to_SG   <- lapply(seq_len(nJ), function(j) which(Ldx    %in% dest_lanes(j)))
to_SP   <- lapply(seq_len(nJ), function(j) which(spotLdx %in% dest_lanes(j)) + nLanes)

## Reservation cost: 0.2 × mean contracted transport rate
v_b    <- rep(mean(inst$transport_coef) * 0.2, (nCS + nCO) * tau)
x0_b   <- as.numeric(inst$carrier_capacity)   # length = (nCS+nCO)*tau
cap_max <- max(x0_b) * 2.0
nK <- nCS + nCO

cat(sprintf("  x0 range: [%.1f, %.1f]  v=%.4f  v·x0=%.1f\n",
            min(x0_b), max(x0_b), v_b[1], sum(v_b * x0_b)))

## Lambda / corrmat for SDDP
lambda_vec <- rep(LAMBDA_VAL, nOD)
set.seed(SEED)
corrmat <- MSTP::gen_corrmat(n_blocks=2L, block_size=max(nI, nJ),
                              cross_corr=CROSS_CORR)

## ── LP proxy: build multi-period deterministic LP ─────────────────────────────
##
## Variables per period t (n_pp total):
##   entry[1:nI]              — entry inventory at END of t
##   exitp[1:nJ]              — exit positive stock at END of t
##   exitm[1:nJ]              — exit shortage at END of t
##   move[1:nLanes]            — contracted carrier allocations in t
##   move_sp[1:nSpotLanes]     — spot carrier allocations in t
##
## Capacity parameters x are embedded as constraint RHS (updated per solve).

n_pp  <- nI + 2L * nJ + nLanes + nSpotLanes
nVars <- tau * n_pp

## Variable index functions (1-based)
var_e  <- function(i, t) (t - 1L) * n_pp + i
var_xp <- function(j, t) (t - 1L) * n_pp + nI + j
var_xm <- function(j, t) (t - 1L) * n_pp + nI + nJ + j
var_m  <- function(l, t) (t - 1L) * n_pp + nI + 2L * nJ + l
var_sp <- function(l, t) (t - 1L) * n_pp + nI + 2L * nJ + nLanes + l

## Spot cost index: spot_coef[(k-1)*nI*nJ*tau + (l-1)*tau + t]
nIJ <- nI * nJ
spot_cost_idx <- function(k, l, t) (k - 1L) * nIJ * tau + (l - 1L) * tau + t

## ── Objective ─────────────────────────────────────────────────────────────────
## Holding + shortage costs on .in (previous-period state); transport/spot on moves.

obj_lp <- numeric(nVars)

for (t in seq_len(tau)) {
  ## Holding costs on .in = state at START of period = end of t-1
  if (t > 1L) {
    for (i in seq_len(nI)) obj_lp[var_e(i,  t - 1L)] <- obj_lp[var_e(i,  t - 1L)] + inst$entry_store_coef[i]
    for (j in seq_len(nJ)) obj_lp[var_xp(j, t - 1L)] <- obj_lp[var_xp(j, t - 1L)] + inst$exit_store_coef[j]
    for (j in seq_len(nJ)) obj_lp[var_xm(j, t - 1L)] <- obj_lp[var_xm(j, t - 1L)] + inst$exit_short_coef[j]
  }
  ## Transport costs
  for (l in seq_len(nLanes))
    obj_lp[var_m(l, t)] <- inst$transport_coef[l]
  ## Spot costs
  for (k in seq_len(nCO)) for (l in seq_len(nIJ))
    obj_lp[var_sp((k - 1L) * nIJ + l, t)] <- inst$spot_coef[spot_cost_idx(k, l, t)]
}

## ── Constraints (all except carrier capacity) ─────────────────────────────────
## We build sparse row format: three vectors i_row, j_col, val + sense + rhs.

row_list <- list()
add_row <- function(vars, vals, sense, rhs) {
  row_list[[length(row_list) + 1L]] <<- list(v = vars, a = vals, s = sense, b = rhs)
}

for (t in seq_len(tau)) {
  q_t <- LAMBDA_VAL   # mean inflow at each origin
  d_t <- LAMBDA_VAL   # mean outflow at each destination

  ## Entry balance: entry(i,t) - entry(i,t-1) + sum(moves_out) = q[i]
  for (i in seq_len(nI)) {
    vars <- var_e(i, t); vals <- 1.0
    if (t > 1L) { vars <- c(vars, var_e(i, t - 1L)); vals <- c(vals, -1.0) }
    sg  <- sapply(from_SG[[i]], var_m,  t = t)
    spv <- sapply(from_SP[[i]] - nLanes, var_sp, t = t)
    if (length(sg))  { vars <- c(vars, sg);  vals <- c(vals,  rep(1.0, length(sg)))  }
    if (length(spv)) { vars <- c(vars, spv); vals <- c(vals,  rep(1.0, length(spv))) }
    add_row(vars, vals, "=", q_t)
  }

  ## Exit balance: exitp(j,t) - exitm(j,t) - exitp(j,t-1) + exitm(j,t-1) - sum(inbound) = -d[j]
  for (j in seq_len(nJ)) {
    vars <- c(var_xp(j, t), var_xm(j, t)); vals <- c(1.0, -1.0)
    if (t > 1L) { vars <- c(vars, var_xp(j, t - 1L), var_xm(j, t - 1L)); vals <- c(vals, -1.0, 1.0) }
    sg  <- sapply(to_SG[[j]], var_m,  t = t)
    spv <- sapply(to_SP[[j]] - nLanes, var_sp, t = t)
    if (length(sg))  { vars <- c(vars, sg);  vals <- c(vals, rep(-1.0, length(sg)))  }
    if (length(spv)) { vars <- c(vars, spv); vals <- c(vals, rep(-1.0, length(spv))) }
    add_row(vars, vals, "=", -d_t)
  }

  ## Exit storage capacity: exitp(j,t-1) + sum(inbound at t) <= exit_capacity[j]
  if (t > 1L) {
    for (j in seq_len(nJ)) {
      vars <- var_xp(j, t - 1L); vals <- 1.0
      sg  <- sapply(to_SG[[j]], var_m,  t = t)
      spv <- sapply(to_SP[[j]] - nLanes, var_sp, t = t)
      if (length(sg))  { vars <- c(vars, sg);  vals <- c(vals, rep(1.0, length(sg)))  }
      if (length(spv)) { vars <- c(vars, spv); vals <- c(vals, rep(1.0, length(spv))) }
      add_row(vars, vals, "<=", inst$exit_capacity[j])
    }
  }
}

## Convert row list to sparse matrix format for highs
n_rows_base <- length(row_list)
cat(sprintf("  LP: %d vars  |  %d base rows (+%d cap rows per solve)\n",
            nVars, n_rows_base, nK * tau))

## Capacity constraint row indices (for dual extraction) — appended after base rows
cap_row_meta <- vector("list", nK * tau)
cap_idx <- 0L
for (t in seq_len(tau)) {
  for (k in seq_len(nCS)) {
    cap_idx <- cap_idx + 1L
    cidx <- CarrierIdx_contr[[k]]
    cap_row_meta[[cap_idx]] <- list(k = k, t = t,
                                    vars = sapply(cidx, var_m, t = t))
  }
  for (k in seq_len(nCO)) {
    cap_idx <- cap_idx + 1L
    kk <- nCS + k
    cidx <- CarrierIdx_spot[[k]]
    cap_row_meta[[cap_idx]] <- list(k = kk, t = t,
                                    vars = sapply(cidx, var_sp, t = t))
  }
}
n_cap_rows <- length(cap_row_meta)

## Build sparse constraint matrix (base rows only; cap rows appended per solve)
all_rows     <- c(row_list, cap_row_meta)
n_total_rows <- length(all_rows)

## Pre-build triplet for the FIXED part of the constraint matrix
trip_i <- integer(); trip_j <- integer(); trip_v <- double()
lhs_v  <- double();  rhs_v  <- double()

for (ri in seq_len(n_total_rows)) {
  r <- all_rows[[ri]]
  trip_i <- c(trip_i, rep(ri, length(r$v)))
  trip_j <- c(trip_j, r$v)
  trip_v <- c(trip_v, r$a)
  if (ri <= n_rows_base) {
    if (r$s == "=")  { lhs_v <- c(lhs_v, r$b);  rhs_v <- c(rhs_v, r$b)  }
    if (r$s == "<=") { lhs_v <- c(lhs_v, -Inf); rhs_v <- c(rhs_v, r$b)  }
    if (r$s == ">=") { lhs_v <- c(lhs_v, r$b);  rhs_v <- c(rhs_v, Inf)  }
  } else {
    lhs_v <- c(lhs_v, -Inf)
    rhs_v <- c(rhs_v, 0.0)   # placeholder; updated per solve
  }
}

A_base <- sparseMatrix(i = trip_i, j = trip_j, x = trip_v,
                       dims = c(n_total_rows, nVars))

## Variable bounds
lb_v <- rep(0.0, nVars)
ub_v <- rep(Inf, nVars)

## ── LP solve function with dual extraction ────────────────────────────────────
## highs_solve returns solver_msg$row_dual[i] = dual for constraint i.
## For ≤ constraints in minimization, dual ≤ 0 (shadow price = ∂obj/∂rhs).
## grad_x[k,t] = row_dual[cap_constraint_of_(k,t)] + v_b[(k-1)*tau+t]

solve_lp_dual <- function(x_cap) {
  ## Update capacity constraint RHS
  rhs_local <- rhs_v
  for (ci in seq_len(n_cap_rows)) {
    k <- cap_row_meta[[ci]]$k
    t <- cap_row_meta[[ci]]$t
    rhs_local[n_rows_base + ci] <- x_cap[(k - 1L) * tau + t]
  }

  sol <- highs_solve(
    L     = obj_lp,
    lower = lb_v, upper = ub_v,
    A     = A_base,
    lhs   = lhs_v, rhs = rhs_local,
    control = highs_control(
      output_flag          = FALSE,
      presolve             = "on",
      simplex_strategy     = 4L,
      simplex_scale_strategy = 4L
    )
  )

  if (sol$status_message != "Optimal") {
    warning("LP status: ", sol$status_message)
    return(NULL)
  }

  ## Duals for capacity rows (shadow prices = ∂obj/∂rhs, ≤ 0 for ≤ constraints)
  row_duals <- sol$solver_msg$row_dual
  duals_cap <- row_duals[n_rows_base + seq_len(n_cap_rows)]

  ## Assemble into (nK * tau) vector indexed (k-1)*tau + t
  duals_x <- numeric(nK * tau)
  for (ci in seq_len(n_cap_rows)) {
    k <- cap_row_meta[[ci]]$k
    t <- cap_row_meta[[ci]]$t
    duals_x[(k - 1L) * tau + t] <- duals_cap[ci]
  }

  list(obj = sol$objective_value, duals_x = duals_x)
}

## ── Step 1: LP solve at x0 ────────────────────────────────────────────────────

cat("\nStep 1: LP proxy solve at x0...\n")
t0 <- proc.time()[["elapsed"]]
lp0 <- solve_lp_dual(x0_b)
t_lp0 <- proc.time()[["elapsed"]] - t0
if (is.null(lp0)) stop("LP solve at x0 failed")
cat(sprintf("  LP(x0) = %.2f  |  v·x0 = %.1f  |  total = %.1f  (%.2fs)\n",
            lp0$obj, sum(v_b * x0_b), lp0$obj + sum(v_b * x0_b), t_lp0))
cat(sprintf("  |∇| at x0 = %.3f\n", sqrt(sum((lp0$duals_x + v_b)^2))))

## ── Step 2: L-BFGS-B over x ───────────────────────────────────────────────────

cat(sprintf("\nStep 2: L-BFGS-B optimisation (maxit=%d)...\n", LP_MAXIT))

dual_cache <- new.env(parent = emptyenv())
dual_cache$x    <- NULL
dual_cache$result <- NULL

f_lp <- function(x) {
  res <- solve_lp_dual(x)
  if (is.null(res)) return(Inf)
  dual_cache$x      <- x
  dual_cache$result <- res
  res$obj + sum(v_b * x)
}

gr_lp <- function(x) {
  if (!isTRUE(all.equal(dual_cache$x, x))) {
    res <- solve_lp_dual(x)
    if (is.null(res)) return(rep(0, length(x)))
    dual_cache$x      <- x
    dual_cache$result <- res
  }
  dual_cache$result$duals_x + v_b
}

t0_lp <- proc.time()[["elapsed"]]
opt_lp <- optim(x0_b, f_lp, gr = gr_lp,
                method  = "L-BFGS-B",
                lower   = CAP_FLOOR, upper = cap_max,
                control = list(maxit = LP_MAXIT, factr = 1e3, pgtol = 1e-9))
t_lp <- proc.time()[["elapsed"]] - t0_lp

x_lp      <- round(pmax(CAP_FLOOR, opt_lp$par))
lp_opt    <- f_lp(x_lp)
cat(sprintf("  LP obj(x*LP) + v·x*LP = %.1f  (%.1f%% vs x0, %.1fs, %d evals)\n",
            lp_opt,
            100 * (1 - lp_opt / (lp0$obj + sum(v_b * x0_b))),
            t_lp, opt_lp$counts[1]))
cat(sprintf("  x*LP range: [%.1f, %.1f]  mean=%.2f\n",
            min(x_lp), max(x_lp), mean(x_lp)))

## ── SDDP evaluation helper ────────────────────────────────────────────────────

eval_sddp <- function(x_vec, label) {
  cat(sprintf("\n  Evaluating %s (%d SDDP iters × %d trials)...\n",
              label, SDDP_ITER, N_TRIALS))
  inst_x <- MSTP::mstp_update_capacity(inst, x_vec)
  config  <- MSTP::mstp_config(inst_x, lambda=lambda_vec, corrmat=corrmat,
                                n_scenarios=N_SCENARIOS)
  for (.attempt in seq_len(6L)) {
    result <- tryCatch({
      t0    <- proc.time()[["elapsed"]]
      model <- MSTP::mstp_train(config, iterations=SDDP_ITER)
      t_tr  <- proc.time()[["elapsed"]] - t0
      lb    <- MSTP::mstp_bound(model)
      sims  <- MSTP::mstp_simulate(model, config, trials=N_TRIALS)
      ub    <- mean(sims$obj); ub_sd <- sd(sims$obj)
      total <- ub + sum(v_b * x_vec)
      cat(sprintf("    LB=%.0f  UB=%.0f±%.0f  UB+v·x=%.0f  (%.1fs)\n",
                  lb, ub, ub_sd, total, t_tr))
      rm(sims, model); gc(verbose=FALSE)
      list(lb=lb, ub=ub, ub_sd=ub_sd, total=total)
    }, error = function(e) {
      cat(sprintf("    [attempt %d failed: %s]\n", .attempt, conditionMessage(e)))
      NULL
    })
    if (!is.null(result)) return(result)
  }
  stop("eval_sddp failed after 6 attempts for ", label)
}

## ── Step 3: Evaluate x0 ──────────────────────────────────────────────────────

cat("\nStep 3: SDDP evaluation at x0\n")
r0 <- eval_sddp(x0_b, "x0 (default)")

## ── Step 4: Evaluate x*_LP ───────────────────────────────────────────────────

cat("\nStep 4: SDDP evaluation at x*_LP\n")
r_lp <- eval_sddp(x_lp, "x*_LP (LP proxy)")

## ── Step 5: SDDP projected gradient from x0 ──────────────────────────────────

cat(sprintf("\nStep 5: SDDP projected gradient (%d iters, α₀=%.1f/√k)...\n",
            N_CAP_ITER, STEP_SIZE))

train_step <- function(x) {
  inst_x  <- MSTP::mstp_update_capacity(inst, x)
  config  <- MSTP::mstp_config(inst_x, lambda=lambda_vec, corrmat=corrmat,
                                n_scenarios=N_SCENARIOS)
  for (.attempt in seq_len(6L)) {
    result <- tryCatch({
      model  <- MSTP::mstp_train(config, SDDP_ITER)
      lb     <- MSTP::mstp_bound(model)
      duals  <- MSTP::mstp_capacity_duals(model, config, n_samples=N_CAP_SAMPLES)
      grad   <- duals + v_b
      list(lb=lb, obj=lb + sum(v_b * x), grad=grad,
           grad_norm=sqrt(sum(grad^2)))
    }, error = function(e) { NULL })
    if (!is.null(result)) return(result)
  }
  stop("train_step failed")
}

step <- train_step(x0_b)
x    <- x0_b
conv_hist <- vector("list", N_CAP_ITER + 1L)
conv_hist[[1]] <- list(iter=0L, obj=step$obj, lb=step$lb,
                        grad_norm=step$grad_norm, x_change=NA_real_)
cat(sprintf("  iter  0 | LB=%9.0f | obj=%9.0f | |∇|=%7.3f\n",
            step$lb, step$obj, step$grad_norm))

for (k in seq_len(N_CAP_ITER)) {
  alpha_k  <- STEP_SIZE / sqrt(k)
  g_hat    <- step$grad / (step$grad_norm + 1e-12)
  x_new    <- round(pmax(CAP_FLOOR, pmin(cap_max, x - alpha_k * g_hat)))
  x_change <- sqrt(sum((x_new - x)^2))

  t0s  <- proc.time()[["elapsed"]]
  step <- train_step(x_new)
  ts   <- proc.time()[["elapsed"]] - t0s

  conv_hist[[k + 1L]] <- list(iter=k, obj=step$obj, lb=step$lb,
                               grad_norm=step$grad_norm, x_change=x_change)
  cat(sprintf("  iter %2d | LB=%9.0f | obj=%9.0f | |∇|=%7.3f | |Δx|=%6.1f | %.1fs\n",
              k, step$lb, step$obj, step$grad_norm, x_change, ts))
  x <- x_new
  gc(verbose=FALSE)
}
x_sddp <- x

cat("\nStep 6: SDDP evaluation at x*_SDDP\n")
r_sddp <- eval_sddp(x_sddp, "x*_SDDP (stochastic duals)")

## ── Summary table ─────────────────────────────────────────────────────────────

red <- function(a, b) sprintf("%+.1f%%", 100 * (b - a) / abs(a))

cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("LP PROXY vs SDDP DUALS — 6×6×20, τ=12, λ=50\n")
cat(strrep("=", 70), "\n")
cat(sprintf("%-24s  %10s  %10s  %10s  %8s\n",
            "Method", "SDDP UB", "UB_sd", "UB + v·x", "vs x0"))
cat(strrep("-", 70), "\n")
for (r in list(list(label="x0 (default)",          res=r0),
               list(label="x*_LP (LP proxy)",       res=r_lp),
               list(label="x*_SDDP (SDDP duals)",   res=r_sddp))) {
  vs <- if (identical(r$res, r0)) "  —" else red(r0$total, r$res$total)
  cat(sprintf("%-24s  %10.0f  %10.0f  %10.0f  %8s\n",
              r$label, r$res$ub, r$res$ub_sd, r$res$total, vs))
}
cat(strrep("-", 70), "\n")

lp_adv <- 100 * (r_lp$total - r_sddp$total) / abs(r0$total)
cat(sprintf("SDDP advantage over LP proxy: %.1f pp\n", lp_adv))

## ── LaTeX snippet ─────────────────────────────────────────────────────────────

cat("\n%% LaTeX snippet for capacity opt comparison table (§6.5)\n")
cat("Method & SDDP UB & $v\\cdot x$ & UB $+ v\\cdot x$ & vs.\\ $x^0$ \\\\\n")
cat("\\hline\n")
for (r in list(list(label="$x^0$ (default)",                  res=r0,    x=x0_b),
               list(label="$x^*_{\\mathrm{LP}}$ (LP proxy)",  res=r_lp,  x=x_lp),
               list(label="$x^*_{\\mathrm{SDDP}}$ (SDDP duals)", res=r_sddp, x=x_sddp))) {
  vx   <- sum(v_b * r$x)
  pct  <- if (identical(r$res, r0)) "---" else sprintf("%+.1f\\%%", 100*(r$res$total - r0$total)/abs(r0$total))
  cat(sprintf("%s & $%.0f$ & $%.0f$ & $%.0f$ & $%s$ \\\\\n",
              r$label, r$res$ub, vx, r$res$total, pct))
}

## ── Persist ───────────────────────────────────────────────────────────────────

dir.create("demo/results", showWarnings=FALSE, recursive=TRUE)
saveRDS(list(r0=r0, r_lp=r_lp, r_sddp=r_sddp,
             x0=x0_b, x_lp=x_lp, x_sddp=x_sddp,
             lp_adv=lp_adv, conv_hist=conv_hist),
        "demo/results/32_results.rds")
cat("\nResults saved to demo/results/32_results.rds\n")

invisible(gc(verbose=FALSE))
