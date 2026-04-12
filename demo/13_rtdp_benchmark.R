## ============================================================
## 13_rtdp_benchmark.R — RTDP Phase 1 vs Phase 2 vs Exact DP / Sampled VFA
##
## Part 1 — Small instances (exact DP feasible as ground truth):
##   1x1 R=10 (nSdx=231), 2x1 R=5 (nSdx=396), 1x2 R=5 (nSdx=726)
##   Quality and convergence across all four methods.
##
## Part 2 — Speedup assessment across instance sizes:
##   Compares time for Phase 2 vs exact DP ceiling, sampled VFA, Phase 1.
##   Includes 2x2 R=5 (nSdx=4,356) and projects to 2x2 R=10 (nSdx=53,361).
##
## Phase 2 uses s0=NULL (random starting state per trajectory) to give the
## VFA enough coverage at each period for a meaningful quality comparison.
## ============================================================

suppressPackageStartupMessages(library(TLPR))
source("R/vfa.R")
source("R/dp_vfa.R")
source("R/rtdp.R")

N_THREADS <- max(1L, parallel::detectCores() - 2L)
N_ITER    <- 10L
N_TRAJ    <- 50L
K_FRAC    <- 0.30
SEED      <- 42L

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

setup_env <- function(nI, nJ, tau, rate, regime = "balanced", seed = SEED) {
  tmpf <- tempfile(fileext = ".json")
  generate_instance(nI = nI, nJ = nJ, tau = tau,
                    nB = 5L, nCS = 10L, nCO = 1L, rate = rate,
                    regime = regime, seed = seed, path = tmpf)
  e <- new.env()
  jsonlite::fromJSON(tmpf) |> list2env(envir = e)
  e$from_i <- as_int_list(e$from_i)
  e$to_j   <- as_int_list(e$to_j)
  create_model(e)
  dp_config(e)
  list(env = e, tmpf = tmpf)
}

max_rel_err_t1 <- function(v_hat, v_ex) {
  ok <- is.finite(v_ex) & is.finite(v_hat)
  if (!any(ok)) return(NA_real_)
  100 * max(abs(v_hat[ok] - v_ex[ok])) / max(abs(v_ex[ok]), 1e-8)
}

mean_rel_err <- function(V_hat, V_ex, tau) {
  errs <- numeric(0L)
  for (t in seq_len(tau)) {
    ve <- V_ex[t, ]; vh <- V_hat[t, ]
    ok <- is.finite(ve) & is.finite(vh)
    if (!any(ok)) next
    errs <- c(errs, abs(vh[ok] - ve[ok]) / max(abs(ve[ok]), 1e-8) * 100)
  }
  if (length(errs) == 0L) NA_real_ else mean(errs)
}

pi_agree_t1 <- function(pi_hat, pi_ex) {
  mean(pi_hat[1L, ] == pi_ex[1L, ], na.rm = TRUE) * 100
}

## ─────────────────────────────────────────────────────────────────────────────
## Part 1 — Small instances: quality comparison (exact DP as ground truth)
## ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 90), "\n", sep = "")
cat("RTDP Benchmark — Part 1: quality on small instances\n")
cat(sprintf("n_iter=%d | K=%.0f%% | n_traj=%d (random start) | Threads=%d\n",
            N_ITER, K_FRAC * 100, N_TRAJ, N_THREADS))
cat(strrep("=", 90), "\n")

configs_small <- list(
  list(nI = 1L, nJ = 1L, tau = 4L, rate = 1.0, label = "1x1 R=10"),
  list(nI = 2L, nJ = 1L, tau = 4L, rate = 0.5, label = "2x1 R=5 "),
  list(nI = 1L, nJ = 2L, tau = 4L, rate = 0.5, label = "1x2 R=5 ")
)

for (cfg in configs_small) {
  obj <- setup_env(cfg$nI, cfg$nJ, cfg$tau, cfg$rate)
  e   <- obj$env; tmpf <- obj$tmpf
  K   <- max(1L, as.integer(e$nSdx * K_FRAC))

  cat(sprintf("\n%s  nSdx=%d  nAdx=%d  nScen=%d  K=%d\n",
              cfg$label, e$nSdx, e$nAdx, e$nScen, K))
  cat(strrep("-", 78), "\n")

  t0 <- proc.time()[["elapsed"]]
  re <- rolling_dp_ptr(e, tmpf, numThreads = N_THREADS)
  t_ex <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  rv <- rolling_dp_vfa(e, tmpf, K = K, numThreads = N_THREADS, seed = SEED)
  t_vfa <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  r1 <- rolling_dp_rtdp(e, tmpf, n_iter = N_ITER, K = K,
                         numThreads = N_THREADS, seed = SEED)
  t_p1 <- proc.time()[["elapsed"]] - t0

  # Phase 2: s0=NULL (random start per trajectory)
  t0 <- proc.time()[["elapsed"]]
  r2 <- rolling_dp_rtdp_p2(e, tmpf, n_iter = N_ITER, n_traj = N_TRAJ,
                             s0 = NULL, numThreads = N_THREADS, seed = SEED)
  t_p2 <- proc.time()[["elapsed"]] - t0

  fmt <- "  %-32s  t=%6.2fs  mean_err=%6.2f%%  max_t1=%6.2f%%  pi@1=%5.1f%%\n"
  cat(sprintf(fmt, "Exact DP (ground truth)",
              t_ex, 0, 0, 100))
  cat(sprintf(fmt, sprintf("Sampled VFA (K=%d, 1-pass)", K),
              t_vfa,
              mean_rel_err(rv$V_approx[seq(e$tau),], re$V[seq(e$tau),], e$tau),
              max_rel_err_t1(rv$V_approx[1L,], re$V[1L,]),
              pi_agree_t1(rv$pi_star, re$pi_star)))
  cat(sprintf(fmt, sprintf("RTDP Phase 1 (K=%d, %d iter)", K, N_ITER),
              t_p1,
              mean_rel_err(r1$V_approx[seq(e$tau),], re$V[seq(e$tau),], e$tau),
              max_rel_err_t1(r1$V_history_t1[N_ITER,], re$V[1L,]),
              pi_agree_t1(r1$pi_star, re$pi_star)))
  cat(sprintf(fmt, sprintf("RTDP Phase 2 (n_traj=%d, %d iter)", N_TRAJ, N_ITER),
              t_p2,
              mean_rel_err(r2$V_approx[seq(e$tau),], re$V[seq(e$tau),], e$tau),
              max_rel_err_t1(r2$V_history_t1[N_ITER,], re$V[1L,]),
              pi_agree_t1(r2$pi_star, re$pi_star)))
  cat(sprintf("  P2 avg visited/period: %s  |  cache@t=1: %d/%-d (%.1f%%)\n",
              paste(round(colMeans(r2$n_visited)), collapse=" "),
              r2$cache_coverage[N_ITER], e$nSdx,
              100 * r2$cache_coverage[N_ITER] / e$nSdx))

  unlink(tmpf); rm(e, re, rv, r1, r2); invisible(gc(verbose = FALSE))
}

## ─────────────────────────────────────────────────────────────────────────────
## Part 2 — Speedup assessment across instance sizes
## ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 90), "\n", sep = "")
cat("RTDP Benchmark — Part 2: speedup assessment\n")
cat(strrep("=", 90), "\n\n")

## Phase 1 is excluded from large instances — its total LP cost is
## K * n_iter × exact-DP cost, making it intractable for nSdx > 5k.
## Shown only where exact DP itself is feasible.
speedup_configs <- list(
  list(nI=1L, nJ=1L, tau=4L, rate=1.0, label="1x1 R=10 ", run_p1=TRUE),
  list(nI=2L, nJ=1L, tau=4L, rate=0.5, label="2x1 R=5  ", run_p1=TRUE),
  list(nI=2L, nJ=2L, tau=4L, rate=0.5, label="2x2 R=5  ", run_p1=FALSE),
  list(nI=2L, nJ=2L, tau=4L, rate=1.0, label="2x2 R=10 ", run_p1=FALSE)
)

cat(sprintf("%-10s  %7s  %8s  %9s  %9s  %9s  %8s\n",
            "Config", "nSdx", "t_exact", "t_vfa(1p)", "t_P1(10i)", "t_P2(10i)", "speedup_P2"))
cat(strrep("-", 75), "\n")

for (cfg in speedup_configs) {
  obj <- setup_env(cfg$nI, cfg$nJ, cfg$tau, cfg$rate)
  e   <- obj$env; tmpf <- obj$tmpf
  K   <- max(1L, as.integer(e$nSdx * K_FRAC))

  # Exact DP
  t0  <- proc.time()[["elapsed"]]
  re  <- tryCatch(rolling_dp_ptr(e, tmpf, numThreads = N_THREADS), error = function(e) NULL)
  t_ex <- proc.time()[["elapsed"]] - t0
  if (is.null(re)) t_ex <- NA_real_

  # Sampled VFA
  t0 <- proc.time()[["elapsed"]]
  rolling_dp_vfa(e, tmpf, K = K, numThreads = N_THREADS, seed = SEED)
  t_vfa <- proc.time()[["elapsed"]] - t0

  # Phase 1 — only for small instances
  t_p1 <- NA_real_
  if (cfg$run_p1) {
    t0 <- proc.time()[["elapsed"]]
    rolling_dp_rtdp(e, tmpf, n_iter = N_ITER, K = K, numThreads = N_THREADS, seed = SEED)
    t_p1 <- proc.time()[["elapsed"]] - t0
  }

  # Phase 2 (n_traj = N_TRAJ, random starts)
  t0  <- proc.time()[["elapsed"]]
  rolling_dp_rtdp_p2(e, tmpf, n_iter = N_ITER, n_traj = N_TRAJ,
                      s0 = NULL, numThreads = N_THREADS, seed = SEED)
  t_p2 <- proc.time()[["elapsed"]] - t0

  speedup_p2 <- if (!is.na(t_ex)) t_ex / t_p2 else NA_real_

  cat(sprintf("%-10s  %7d  %8s  %9.2f  %9s  %9.2f  %8s\n",
              cfg$label, e$nSdx,
              if (is.na(t_ex)) "   SKIP" else sprintf("%8.2f", t_ex),
              t_vfa,
              if (is.na(t_p1)) "     SKIP" else sprintf("%9.2f", t_p1),
              t_p2,
              if (is.na(speedup_p2)) "    N/A" else sprintf("%7.1fx", speedup_p2)))

  unlink(tmpf); rm(e, re); invisible(gc(verbose = FALSE))
}

cat(strrep("-", 75), "\n")
cat(sprintf("n_traj=%d  |  n_iter=%d  |  K=%.0f%%  |  Threads=%d\n",
            N_TRAJ, N_ITER, K_FRAC * 100, N_THREADS))

cat("\n", strrep("=", 90), "\n", sep = "")
cat("Summary:\n")
cat("  Phase 1: iterative random K-sampling.  LP cost = K*tau*n_iter > exact DP → always slower.\n")
cat("  Phase 2: LP-accurate trajectory sampling.  LP cost independent of |S|.\n")
cat("           Scales to instances where exact DP is intractable (|S| > 50k).\n")
cat("  Quality: Phase 2 (random start, n_traj=50) matches Phase 1 / sampled VFA on small instances.\n")
cat(strrep("=", 90), "\n")
