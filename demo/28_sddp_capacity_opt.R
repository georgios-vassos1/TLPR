## ============================================================
## 28_sddp_capacity_opt.R — SDDP-based capacity optimisation (P1.5)
##
## Part A — Small topologies (2×1, 1×2, 2×2 R=5, 2×2 R=10):
##   LP dual gradient (demo/14 approach) → x*; then SDDP validates
##   that the LP-guided x* also reduces stochastic cost.
##
##   Step 1: Generate TLPR instance
##   Step 2: Build multi-period LP at mode scenario
##   Step 3: L-BFGS-B with LP dual gradient → x*
##   Step 4: Train SDDP at x0 and x* (same corrmat)
##   Step 5: Report LP and SDDP cost reductions
##
## Part B — 6×6×20 instance (tau=12):
##   SDDP capacity duals → first-stage optimisation (no LP proxy).
##
##   Step 1: Load pre-generated MSTP instance
##   Step 2: Train SDDP at x0; extract capacity duals → ∂V/∂x
##   Step 3: L-BFGS-B with SDDP dual gradient + reservation cost → x*
##   Step 4: Retrain SDDP at x*; compare UB_SDDP(x0) vs UB_SDDP(x*)
##
## Addresses: G5 (strategic–tactical capacity optimisation)
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(MSTP)
})

## ── Shared parameters ────────────────────────────────────────────────────────

ALPHA_SCALE <- 3.0
SEED        <- 42L
SDDP_ITER   <- 1000L
N_TRIALS    <- 300L
N_SCENARIOS <- 10L
LAMBDA_VAL  <- 700.0
CROSS_CORR  <- 0.4

## Part A configs (TLPR instances, tau=4, LP proxy)
CONFIGS <- list(
  list(nI = 2L, nJ = 1L, nCS = 10L, rate = 0.5, label = "2×1, R=5"),
  list(nI = 1L, nJ = 2L, nCS = 10L, rate = 0.5, label = "1×2, R=5"),
  list(nI = 2L, nJ = 2L, nCS = 10L, rate = 0.5, label = "2×2, R=5"),
  list(nI = 2L, nJ = 2L, nCS = 10L, rate = 1.0, label = "2×2, R=10")
)

## Part B parameters (MSTP instance, projected gradient + warm-start)
N_CAP_ITER      <- 10L   # projected-gradient outer iterations
N_CAP_SAMPLES   <- 200L  # trajectories for dual averaging per outer step
SDDP_ITER_COLD  <- 500L  # SDDP iterations for cold start (iter 0)
SDDP_ITER_WARM  <- 100L  # SDDP iterations for warm-start steps (iters 1+)
STEP_SIZE       <- 3.0   # α₀ for diminishing step αₖ = α₀ / √k

## ── Helpers ──────────────────────────────────────────────────────────────────

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

## Adapter: TLPR env → MSTP instance list for mstp_config()
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

## Build MSTP instance with updated Cb/Co (after capacity optimisation)
update_capacity <- function(e, x_star) {
  x_mat <- matrix(x_star, nrow = e$tau, ncol = e$nCS + e$nCO, byrow = TRUE)
  e_new <- list2env(as.list(e), parent = emptyenv())
  e_new$Cb <- x_mat[, seq(e$nCS), drop = FALSE]
  e_new$Co <- x_mat[, e$nCS + seq(e$nCO), drop = FALSE]
  e_new
}

## Train SDDP and simulate; corrmat must be pre-generated for x0/x* consistency.
run_sddp <- function(env, corrmat) {
  lambda_vec <- rep(LAMBDA_VAL, env$nI + env$nJ)

  inst   <- tlpr_env_to_mstp_inst(env)
  config <- MSTP::mstp_config(inst,
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

  list(lb    = lb,
       ub    = mean(sims$obj),
       ub_sd = sd(sims$obj),
       t_tr  = t_tr,
       t_si  = t_si)
}

## ── Header ────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 70), "\n", sep = "")
cat("28_sddp_capacity_opt.R — SDDP-based capacity optimisation (P1.5)\n")
cat(strrep("=", 70), "\n\n")

## ══════════════════════════════════════════════════════════════════════════════
## PART A — Small topologies: LP proxy + SDDP validation
## ══════════════════════════════════════════════════════════════════════════════

cat(strrep("=", 70), "\n")
cat("PART A — LP proxy + SDDP validation (small topologies)\n")
cat(strrep("=", 70), "\n")

results <- vector("list", length(CONFIGS))

for (ci in seq_along(CONFIGS)) {
  cfg <- CONFIGS[[ci]]

  cat(sprintf("\n%s\n[%s]\n%s\n", strrep("=", 60), cfg$label, strrep("=", 60)))

  ## ── Step 1: Generate TLPR instance ──────────────────────────────────────────
  tmpf <- tempfile(fileext = ".json")
  tryCatch(
    TLPR::generate_instance(
      nI = cfg$nI, nJ = cfg$nJ, tau = 4L,
      nB = 5L, nCS = cfg$nCS, nCO = 1L,
      rate = cfg$rate, regime = "balanced",
      seed = SEED, path = tmpf
    ),
    error   = function(e) { unlink(tmpf); stop(e) },
    finally = NULL
  )
  on.exit(unlink(tmpf), add = TRUE)

  env <- new.env()
  jsonlite::fromJSON(tmpf) |> list2env(envir = env)
  env$alpha  <- ALPHA_SCALE * env$alpha
  env$from_i <- as_int_list(env$from_i)
  env$to_j   <- as_int_list(env$to_j)
  env$CTb    <- env$CTb / 5.0

  cat(sprintf("  nSdx=%s  nCS=%d  nCO=%d  τ=%d  R=%d\n",
              format(env$nSdx, big.mark = ","),
              env$nCS, env$nCO, env$tau, env$R))

  ## Corrmat generated once; shared by x0 and x* evaluations for apples-to-apples
  nOD <- cfg$nI + cfg$nJ
  set.seed(SEED)
  corrmat <- MSTP::gen_corrmat(n_blocks   = 2L,
                                block_size = max(cfg$nI, cfg$nJ),
                                cross_corr = CROSS_CORR)
  if (nrow(corrmat) > nOD) corrmat <- corrmat[seq(nOD), seq(nOD)]

  ## ── Step 2: Reference scenario for LP ───────────────────────────────────────
  q_mode    <- env$Q$vals[which.max(env$Q$prob)]
  d_mode    <- env$D$vals[which.max(env$D$prob)]
  w_mode    <- env$W$vals[which.max(env$W$prob)]
  Q_ref     <- rep(q_mode, env$tau)
  D_ref     <- rep(d_mode, env$tau)
  env$CTo[] <- w_mode

  ## ── Step 3: Build multi-period LP and optimise capacity ─────────────────────
  ccx  <- carrier_capacity_padded(env)
  tlx  <- transition_logic(env, q = Q_ref[seq(env$nI)], d = D_ref[seq(env$nJ)])
  slx  <- storage_limits(env,  q = Q_ref[seq(env$nI)])
  obj_ <- c(env$alpha, env$CTb, env$CTo[1L, ], env$alpha)

  A   <- rbind(ccx$A, tlx$A, slx$A)
  rhs <- c(ccx$rhs, tlx$rhs, slx$rhs)
  sns <- c(ccx$sense, tlx$sense, slx$sense)

  mp <- multiperiod_expansion(env, Q_ref, D_ref, A, obj_, rhs, sns)
  mp$modelsense <- "min"
  mp$vtype      <- rep("C", ncol(mp$A))

  offset <- env$nI + 2L * env$nJ
  n_col  <- ncol(mp$A)
  A_s0   <- cbind(Matrix::Diagonal(offset),
                  Matrix::Matrix(0L, nrow = offset, ncol = n_col - offset))
  mp$A     <- rbind(mp$A, A_s0)
  mp$rhs   <- c(mp$rhs,   numeric(offset))
  mp$sense <- c(mp$sense, rep("=", offset))

  n_per_period <- nrow(A)
  capdx <- c(outer(seq(env$nCS + env$nCO),
                   seq(0L, n_per_period * (env$tau - 1L), n_per_period), "+"))

  x0 <- c(t(cbind(env$Cb, env$Co)))
  v  <- rep(mean(env$CTb) * 2.0, length(x0))

  dual_cache      <- new.env(parent = emptyenv())
  dual_cache$x    <- NULL
  dual_cache$dual <- NULL

  f_dual <- function(model, x, v, capdx) {
    model$rhs[capdx] <- x
    opt <- solve_lp(model)
    if (is.null(opt$objval)) return(Inf)
    dual_cache$x    <- x
    dual_cache$dual <- opt$dual
    opt$objval + sum(v * x)
  }

  gr_dual <- function(model, x, v, capdx) {
    if (!isTRUE(all.equal(dual_cache$x, x))) {
      model$rhs[capdx] <- x
      opt <- solve_lp(model)
      dual_cache$x    <- x
      dual_cache$dual <- opt$dual
    }
    dual_cache$dual[capdx] + v
  }

  v0_lp <- f_dual(mp, x0, v, capdx)
  cat(sprintf("  LP V(x0) = %.2f\n", v0_lp))

  ctrl  <- list(maxit = 1000L, factr = 1e3, pgtol = 1e-9)
  t0_lp <- proc.time()[["elapsed"]]
  res   <- optim(x0, f_dual, gr = gr_dual,
                 model = mp, v = v, capdx = capdx,
                 method = "L-BFGS-B", lower = 0.0, upper = 10.0, control = ctrl)
  t_lp  <- proc.time()[["elapsed"]] - t0_lp

  x_star   <- round(res$par)
  vstar_lp <- f_dual(mp, x_star, v, capdx)
  lp_reduction <- (1.0 - vstar_lp / v0_lp) * 100.0

  cat(sprintf("  LP V(x*) = %.2f  |  LP reduction = %.1f%%  |  LP time = %.1fs\n",
              vstar_lp, lp_reduction, t_lp))
  cat(sprintf("  x0    = %s\n", paste(x0,     collapse = " ")))
  cat(sprintf("  x*    = %s\n", paste(x_star, collapse = " ")))

  ## ── Step 4: SDDP at x0 ──────────────────────────────────────────────────────
  cat(sprintf("  SDDP at x0 (%d iters × %d trials)...\n", SDDP_ITER, N_TRIALS))
  sddp0      <- run_sddp(env, corrmat)
  ub0_with_v <- sddp0$ub + sum(v * x0)
  cat(sprintf("  SDDP(x0): LB=%.2f  UB=%.2f±%.2f  (UB+v·x0=%.2f)  Train=%.1fs\n",
              sddp0$lb, sddp0$ub, sddp0$ub_sd, ub0_with_v, sddp0$t_tr))

  ## ── Step 5: SDDP at x* ──────────────────────────────────────────────────────
  cat(sprintf("  SDDP at x* (%d iters × %d trials)...\n", SDDP_ITER, N_TRIALS))
  env_star      <- update_capacity(env, x_star)
  sddp_star     <- run_sddp(env_star, corrmat)
  ubstar_with_v <- sddp_star$ub + sum(v * x_star)
  cat(sprintf("  SDDP(x*): LB=%.2f  UB=%.2f±%.2f  (UB+v·x*=%.2f)  Train=%.1fs\n",
              sddp_star$lb, sddp_star$ub, sddp_star$ub_sd, ubstar_with_v, sddp_star$t_tr))

  sddp_reduction <- (1.0 - ubstar_with_v / ub0_with_v) * 100.0
  cat(sprintf("  SDDP reduction (incl. reservation): %.1f%%\n", sddp_reduction))

  results[[ci]] <- list(
    label          = cfg$label,
    nSdx           = env$nSdx,
    x0             = x0,
    x_star         = x_star,
    v0_lp          = v0_lp,
    vstar_lp       = vstar_lp,
    lp_reduction   = lp_reduction,
    t_lp           = t_lp,
    sddp0_ub       = sddp0$ub,
    sddp0_ub_v     = ub0_with_v,
    sddp_star_ub   = sddp_star$ub,
    sddp_star_ub_v = ubstar_with_v,
    sddp_reduction = sddp_reduction,
    t_sddp0        = sddp0$t_tr,
    t_sddp_star    = sddp_star$t_tr
  )

  rm(env_star, sddp0, sddp_star)
  gc(verbose = FALSE)
}

## ── Part A summary ────────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 88), "\n", sep = "")
cat("PART A — LP PROXY + SDDP VALIDATION\n")
cat(strrep("=", 88), "\n")

hdr <- sprintf("%-14s  %8s  %10s  %10s  %8s  %10s  %10s  %8s",
               "Config", "|S|",
               "LP V(x0)", "LP V(x*)", "LP red%",
               "SDDP(x0)+v", "SDDP(x*)+v", "SDDP red%")
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")

for (r in results) {
  cat(sprintf("%-14s  %8s  %10.1f  %10.1f  %7.1f%%  %10.1f  %10.1f  %7.1f%%\n",
              r$label,
              format(r$nSdx, big.mark = ","),
              r$v0_lp, r$vstar_lp, r$lp_reduction,
              r$sddp0_ub_v, r$sddp_star_ub_v, r$sddp_reduction))
}

cat(strrep("-", nchar(hdr)), "\n")
cat(sprintf("SDDP iters=%d  |  OOB trials=%d  |  α-scale=%.1f  |  λ=%.0f  |  seed=%d\n",
            SDDP_ITER, N_TRIALS, ALPHA_SCALE, LAMBDA_VAL, SEED))
cat("LP V(·)   = deterministic multi-period LP cost + reservation cost v·x\n")
cat("SDDP(·)+v = mean SDDP simulation cost + reservation cost v·x\n")

## ══════════════════════════════════════════════════════════════════════════════
## PART B — 6×6×20: first-stage optimisation via SDDP capacity duals
## ══════════════════════════════════════════════════════════════════════════════

cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("PART B — SDDP dual-based capacity optimisation (6×6×20, τ=12)\n")
cat(strrep("=", 70), "\n\n")

  inst <- MSTP::generate_instance(tau=12L, nOrigins=6L, nDestinations=6L,
                                   nCarriers=20L, seed=SEED)
  nI   <- inst$nOrigins
  nJ   <- inst$nDestinations
  nCS  <- inst$nCarriers
  nCSO <- inst$nSpotCarriers
  tau  <- inst$tau
  nOD  <- nI + nJ

  cat(sprintf("Instance: generated (6×6×20, τ=%d, seed=%d)\n", tau, SEED))
  cat(sprintf("  nI=%d  nJ=%d  nCS=%d  nCSO=%d  τ=%d\n", nI, nJ, nCS, nCSO, tau))

  lambda_vec <- rep(LAMBDA_VAL, nOD)
  set.seed(SEED)
  corrmat <- MSTP::gen_corrmat(n_blocks   = 2L,
                                block_size = max(nI, nJ),
                                cross_corr = CROSS_CORR)

  ## Reservation cost: mean transport rate × 2, one entry per (carrier, period)
  v_b <- rep(mean(inst$transport_coef) * 2.0, (nCS + nCSO) * tau)

  ## Capacity bounds for L-BFGS-B
  x0_b    <- as.numeric(inst$carrier_capacity)
  cap_max <- max(x0_b) * 2.0

  ## Cuts file persists across outer iterations for warm-starting
  cuts_file <- tempfile(fileext = ".json")
  on.exit(unlink(cuts_file), add = TRUE)

  ## Helper: train SDDP (cold or warm), write cuts, extract duals
  train_step <- function(x, warm = FALSE, n_iter = SDDP_ITER_WARM) {
    inst_x  <- MSTP::mstp_update_capacity(inst, x)
    config  <- MSTP::mstp_config(inst_x,
                                  lambda      = lambda_vec,
                                  corrmat     = corrmat,
                                  n_scenarios = N_SCENARIOS)
    model   <- if (warm) {
      MSTP::mstp_train_warm(config, n_iter, cuts_file)
    } else {
      MSTP::mstp_train(config, n_iter)
    }
    MSTP::mstp_write_cuts(model, cuts_file)
    lb    <- MSTP::mstp_bound(model)
    duals <- MSTP::mstp_capacity_duals(model, config, n_samples = N_CAP_SAMPLES)
    grad  <- duals + v_b
    list(model = model, config = config, lb = lb,
         obj   = lb + sum(v_b * x),
         grad  = grad, grad_norm = sqrt(sum(grad^2)))
  }

  conv_hist <- vector("list", N_CAP_ITER + 1L)

  ## ── Step 1: Cold start at x0 ────────────────────────────────────────────────
  cat(sprintf("Step 1: cold-start SDDP at x0 (%d iters)...\n", SDDP_ITER_COLD))
  t0_b  <- proc.time()[["elapsed"]]
  step  <- train_step(x0_b, warm = FALSE, n_iter = SDDP_ITER_COLD)
  t_s   <- proc.time()[["elapsed"]] - t0_b
  conv_hist[[1]] <- list(iter = 0L, obj = step$obj, lb = step$lb,
                          grad_norm = step$grad_norm, x_change = NA_real_,
                          t = t_s, warm = FALSE, n_iter = SDDP_ITER_COLD)
  cat(sprintf("  iter  0 | LB=%11.0f | obj=%13.0f | |∇|=%8.3f | Δt=%.1fs [cold]\n",
              step$lb, step$obj, step$grad_norm, t_s))
  x <- x0_b

  ## ── Step 2: Projected gradient with warm-start ──────────────────────────────
  cat(sprintf("\nStep 2: projected gradient (%d iters, αₖ = %.1f/√k, warm SDDP=%d iters)\n",
              N_CAP_ITER, STEP_SIZE, SDDP_ITER_WARM))

  t0_opt <- proc.time()[["elapsed"]]

  for (k in seq_len(N_CAP_ITER)) {
    alpha_k <- STEP_SIZE / sqrt(k)
    x_new   <- round(pmax(0.0, pmin(cap_max, x - alpha_k * step$grad)))
    x_change <- sqrt(sum((x_new - x)^2))

    t0_s  <- proc.time()[["elapsed"]]
    step  <- train_step(x_new, warm = TRUE, n_iter = SDDP_ITER_WARM)
    t_s   <- proc.time()[["elapsed"]] - t0_s

    conv_hist[[k + 1L]] <- list(iter = k, obj = step$obj, lb = step$lb,
                                 grad_norm = step$grad_norm, x_change = x_change,
                                 t = t_s, warm = TRUE, n_iter = SDDP_ITER_WARM)
    cat(sprintf("  iter %2d | LB=%11.0f | obj=%13.0f | |∇|=%8.3f | |Δx|=%7.1f | Δt=%.1fs [warm]\n",
                k, step$lb, step$obj, step$grad_norm, x_change, t_s))

    x <- x_new
    step$model <- NULL; gc(verbose = FALSE)
  }

  t_opt    <- proc.time()[["elapsed"]] - t0_opt
  x_star_b <- x

  ## ── Step 3: Baseline UB at x0 (full simulation) ─────────────────────────────
  cat(sprintf("\nStep 3: full simulation at x0 (%d iters × %d trials)...\n",
              SDDP_ITER, N_TRIALS))
  t0_b    <- proc.time()[["elapsed"]]
  config0 <- MSTP::mstp_config(inst, lambda=lambda_vec, corrmat=corrmat,
                                n_scenarios=N_SCENARIOS)
  model0  <- MSTP::mstp_train(config0, iterations = SDDP_ITER)
  lb0     <- MSTP::mstp_bound(model0)
  sims0   <- MSTP::mstp_simulate(model0, config0, trials = N_TRIALS)
  t_base  <- proc.time()[["elapsed"]] - t0_b
  ub0     <- mean(sims0$obj); ub0_sd <- sd(sims0$obj)
  ub0_with_v <- ub0 + sum(v_b * x0_b)
  cat(sprintf("  LB=%.0f  UB=%.0f±%.0f  UB+v·x0=%.0f  (%.1fs)\n",
              lb0, ub0, ub0_sd, ub0_with_v, t_base))
  rm(sims0, model0); gc(verbose = FALSE)

  ## ── Step 4: Validate at x* ──────────────────────────────────────────────────
  cat(sprintf("\nStep 4: SDDP at x* (%d iters × %d trials)...\n", SDDP_ITER, N_TRIALS))
  t0_b      <- proc.time()[["elapsed"]]
  inst_star <- MSTP::mstp_update_capacity(inst, x_star_b)
  config_s  <- MSTP::mstp_config(inst_star, lambda=lambda_vec, corrmat=corrmat,
                                  n_scenarios=N_SCENARIOS)
  model_s   <- MSTP::mstp_train_warm(config_s, SDDP_ITER, cuts_file)
  lb_star   <- MSTP::mstp_bound(model_s)
  sims_s    <- MSTP::mstp_simulate(model_s, config_s, trials = N_TRIALS)
  t_star    <- proc.time()[["elapsed"]] - t0_b
  ub_star   <- mean(sims_s$obj); ub_star_sd <- sd(sims_s$obj)
  ubstar_with_v <- ub_star + sum(v_b * x_star_b)
  cat(sprintf("  LB=%.0f  UB=%.0f±%.0f  UB+v·x*=%.0f  (%.1fs)\n",
              lb_star, ub_star, ub_star_sd, ubstar_with_v, t_star))
  rm(sims_s, model_s, inst_star); gc(verbose = FALSE)

  ## ── Part B summary ────────────────────────────────────────────────────────────
  sddp_b_reduction <- (1.0 - ubstar_with_v / ub0_with_v) * 100.0

  ## Convergence table
  cat("\n\n", strrep("=", 80), "\n", sep = "")
  cat("PART B — CONVERGENCE HISTORY\n")
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
  cat(sprintf("αₖ = %.1f/√k  |  cap ∈ [0, %.0f]  |  dual samples=%d\n",
              STEP_SIZE, cap_max, N_CAP_SAMPLES))

  ## Final comparison
  cat("\n\n", strrep("=", 70), "\n", sep = "")
  cat("PART B — SDDP DUAL-BASED CAPACITY OPTIMISATION (6×6×20)\n")
  cat(strrep("=", 70), "\n")
  cat(sprintf("%-20s  %14s  %14s  %8s\n",
              "Metric", "x0 (original)", "x* (optimised)", "Change"))
  cat(strrep("-", 62), "\n")
  cat(sprintf("%-20s  %14.0f  %14.0f  %7.1f%%\n",
              "SDDP LB",        lb0,        lb_star,
              100*(lb_star - lb0) / (abs(lb0) + 1e-8)))
  cat(sprintf("%-20s  %14.0f  %14.0f  %7.1f%%\n",
              "SDDP UB (mean)", ub0,        ub_star,
              100*(ub_star - ub0) / (abs(ub0) + 1e-8)))
  cat(sprintf("%-20s  %14.0f  %14.0f  %7.1f%%\n",
              "UB + v·x",       ub0_with_v, ubstar_with_v, sddp_b_reduction))
  cat(strrep("-", 62), "\n")
  cat(sprintf("Cold start: %d iters  |  Warm steps: %d × %d iters  |  Dual samples: %d\n",
              SDDP_ITER_COLD, N_CAP_ITER, SDDP_ITER_WARM, N_CAP_SAMPLES))
  cat(sprintf("Opt wall time: %.1fs  |  Validation: %.1fs\n", t_opt, t_star))
  cat("Gradient: ∂V/∂carrier_capacity from SDDP capacity constraint duals\n")
  cat("Step: projected gradient, αₖ = α₀/√k, projected onto [0, 2×max(x0)]\n")

  ## ── LaTeX snippet ─────────────────────────────────────────────────────────────
  cat("\n%% LaTeX snippet for Part B table (§6.7)\n")
  cat("\\begin{tabular}{lrr}\n")
  cat("Metric & $x^0$ & $x^*$ \\\\\n\\hline\n")
  cat(sprintf("SDDP LB       & $%.0f$ & $%.0f$ \\\\\n", lb0, lb_star))
  cat(sprintf("SDDP UB       & $%.0f$ & $%.0f$ \\\\\n", ub0, ub_star))
  cat(sprintf("UB + $v\\cdot x$ & $%.0f$ & $%.0f$ ($%.1f\\%%$) \\\\\n",
              ub0_with_v, ubstar_with_v, sddp_b_reduction))
  cat("\\end{tabular}\n")

## ── Part A LaTeX ─────────────────────────────────────────────────────────────

cat("\n%% LaTeX snippet for Part A table (§6.7)\n")
cat(paste0("Config & $|\\mathcal{S}|$ & $V_{\\rm LP}(x^0)$ & $V_{\\rm LP}(x^*)$",
           " & LP red.\\% & $V_{\\rm SDDP}(x^0)$ & $V_{\\rm SDDP}(x^*)$",
           " & SDDP red.\\% \\\\\n"))
cat("\\hline\n")
for (r in results)
  cat(sprintf(paste0("%s & %s & $%.1f$ & $%.1f$ & $%.1f\\%%$",
                     " & $%.1f$ & $%.1f$ & $%.1f\\%%$ \\\\\n"),
              r$label,
              format(r$nSdx, big.mark = ","),
              r$v0_lp, r$vstar_lp, r$lp_reduction,
              r$sddp0_ub_v, r$sddp_star_ub_v, r$sddp_reduction))

invisible(gc(verbose = FALSE))
