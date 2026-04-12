## ============================================================
## 13_rtdp_benchmark.R — RTDP Phase 1 vs Phase 2/3 vs Exact DP / Sampled VFA
##
## Part 1 — Small instances (exact DP feasible as ground truth):
##   1x1 R=10 (nSdx=231), 2x1 R=5 (nSdx=396), 1x2 R=5 (nSdx=726)
##   Quality and convergence across all four methods.
##
## Part 2 — Speedup assessment across instance sizes:
##   Compares time for Phase 2 vs exact DP ceiling, sampled VFA, Phase 1.
##   Includes 2x2 R=5 (nSdx=4,356) and 2x2 R=10 (nSdx=53,361).
##
## Part 3 — Phase 3 convergence check (small instance):
##   Runs Phase 2 with tol=1e-3 / min_iter=3 and reports n_done vs n_iter,
##   theta_delta trajectory, and quality vs exact DP.
##
## Part 4 — Quality on large instance (2x2 R=10):
##   Compares Phase 2 quality against exact DP (run offline, ~2 min).
##
## Phase 2/3 uses s0=NULL (random starting state per trajectory) to give the
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

  t0 <- proc.time()[["elapsed"]]
  r2 <- rolling_dp_rtdp_p2(e, tmpf, n_iter = N_ITER, n_traj = N_TRAJ,
                             s0 = NULL, numThreads = N_THREADS, seed = SEED)
  t_p2 <- proc.time()[["elapsed"]] - t0

  fmt <- "  %-38s  t=%6.2fs  mean_err=%6.2f%%  max_t1=%6.2f%%  pi@1=%5.1f%%\n"
  cat(sprintf(fmt, "Exact DP (ground truth)",
              t_ex, 0, 0, 100))
  cat(sprintf(fmt, sprintf("Sampled VFA (K=%d, 1-pass)", K),
              t_vfa,
              mean_rel_err(rv$V_approx[seq(e$tau),], re$V[seq(e$tau),], e$tau),
              max_rel_err_t1(rv$V_approx[1L,], re$V[1L,]),
              pi_agree_t1(rv$pi_star, re$pi_star)))
  cat(sprintf(fmt, sprintf("RTDP Phase 1 (K=%d, %d iter → n_done=%d)", K, N_ITER, r1$n_done),
              t_p1,
              mean_rel_err(r1$V_approx[seq(e$tau),], re$V[seq(e$tau),], e$tau),
              max_rel_err_t1(r1$V_history_t1[r1$n_done,], re$V[1L,]),
              pi_agree_t1(r1$pi_star, re$pi_star)))
  cat(sprintf(fmt, sprintf("RTDP Phase 2/3 (n_traj=%d, %d iter → n_done=%d)", N_TRAJ, N_ITER, r2$n_done),
              t_p2,
              mean_rel_err(r2$V_approx[seq(e$tau),], re$V[seq(e$tau),], e$tau),
              max_rel_err_t1(r2$V_history_t1[r2$n_done,], re$V[1L,]),
              pi_agree_t1(r2$pi_star, re$pi_star)))
  cat(sprintf("  P2 avg visited/period: %s  |  cache@t=1: %d/%-d (%.1f%%)\n",
              paste(round(colMeans(r2$n_visited)), collapse=" "),
              r2$cache_coverage[r2$n_done], e$nSdx,
              100 * r2$cache_coverage[r2$n_done] / e$nSdx))

  unlink(tmpf); rm(e, re, rv, r1, r2); invisible(gc(verbose = FALSE))
}

## ─────────────────────────────────────────────────────────────────────────────
## Part 2 — Speedup assessment across instance sizes
## ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 90), "\n", sep = "")
cat("RTDP Benchmark — Part 2: speedup assessment\n")
cat(strrep("=", 90), "\n\n")

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

  t0  <- proc.time()[["elapsed"]]
  re  <- tryCatch(rolling_dp_ptr(e, tmpf, numThreads = N_THREADS), error = function(e) NULL)
  t_ex <- proc.time()[["elapsed"]] - t0
  if (is.null(re)) t_ex <- NA_real_

  t0 <- proc.time()[["elapsed"]]
  rolling_dp_vfa(e, tmpf, K = K, numThreads = N_THREADS, seed = SEED)
  t_vfa <- proc.time()[["elapsed"]] - t0

  t_p1 <- NA_real_
  if (cfg$run_p1) {
    t0 <- proc.time()[["elapsed"]]
    rolling_dp_rtdp(e, tmpf, n_iter = N_ITER, K = K, numThreads = N_THREADS, seed = SEED)
    t_p1 <- proc.time()[["elapsed"]] - t0
  }

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

## ─────────────────────────────────────────────────────────────────────────────
## Part 3 — Phase 3 convergence check
##   Convergence metric: mean relative change in V_approx[1,] (v_delta).
##   Stopping: max(tail(v_delta, patience)) < tol after min_iter iterations.
##   Tested on two instances: 1x1 R=10 (small) and 2x2 R=5 (medium).
## ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 90), "\n", sep = "")
cat("RTDP Benchmark — Part 3: Phase 3 convergence (tol=1e-2, patience=3, min_iter=3)\n")
cat(strrep("=", 90), "\n")

conv_configs <- list(
  list(nI=1L, nJ=1L, tau=4L, rate=1.0, label="1x1 R=10", max_iter=30L),
  list(nI=2L, nJ=2L, tau=4L, rate=0.5, label="2x2 R=5 ", max_iter=30L)
)

for (cfg3 in conv_configs) {
  obj3 <- setup_env(cfg3$nI, cfg3$nJ, cfg3$tau, cfg3$rate)
  e3   <- obj3$env; tmpf3 <- obj3$tmpf
  K3   <- max(1L, as.integer(e3$nSdx * K_FRAC))

  cat(sprintf("\n%s  nSdx=%d\n", cfg3$label, e3$nSdx))
  cat(strrep("-", 70), "\n")

  has_exact <- e3$nSdx < 10000L
  if (has_exact) re3 <- rolling_dp_ptr(e3, tmpf3, numThreads = N_THREADS)

  # Phase 1 with convergence stopping
  t0  <- proc.time()[["elapsed"]]
  r1c <- rolling_dp_rtdp(e3, tmpf3, n_iter=cfg3$max_iter, K=K3,
                          tol=1e-2, min_iter=3L, patience=3L,
                          numThreads=N_THREADS, seed=SEED)
  t_p1 <- proc.time()[["elapsed"]] - t0

  cat(sprintf("Phase 1  (max=%d): stopped at iter %d  t=%.2fs\n",
              cfg3$max_iter, r1c$n_done, t_p1))
  cat(sprintf("  v_delta: %s\n", paste(sprintf("%.3f", r1c$v_delta), collapse=" ")))
  if (has_exact)
    cat(sprintf("  mean_err=%.2f%%  pi@1=%.1f%%\n",
                mean_rel_err(r1c$V_approx[seq(e3$tau),], re3$V[seq(e3$tau),], e3$tau),
                pi_agree_t1(r1c$pi_star, re3$pi_star)))

  # Phase 2/3 with convergence + ε-greedy
  t0  <- proc.time()[["elapsed"]]
  r2c <- rolling_dp_rtdp_p2(e3, tmpf3, n_iter=cfg3$max_iter, n_traj=N_TRAJ,
                              s0=NULL, epsilon=0.2, tol=1e-2, min_iter=3L, patience=3L,
                              numThreads=N_THREADS, seed=SEED)
  t_p2 <- proc.time()[["elapsed"]] - t0

  cat(sprintf("Phase 2/3 (max=%d): stopped at iter %d  t=%.2fs\n",
              cfg3$max_iter, r2c$n_done, t_p2))
  cat(sprintf("  v_delta: %s\n", paste(sprintf("%.3f", r2c$v_delta), collapse=" ")))
  if (has_exact)
    cat(sprintf("  mean_err=%.2f%%  pi@1=%.1f%%\n",
                mean_rel_err(r2c$V_approx[seq(e3$tau),], re3$V[seq(e3$tau),], e3$tau),
                pi_agree_t1(r2c$pi_star, re3$pi_star)))

  if (has_exact) rm(re3)
  unlink(tmpf3); rm(e3, r1c, r2c); invisible(gc(verbose = FALSE))
}

## ─────────────────────────────────────────────────────────────────────────────
## Part 4 — Quality on large instance (2x2 R=10, nSdx=53,361)
##   Runs exact DP (~2 min) as offline ground truth, then Phase 2/3.
##   Warning: this part takes ~3 min total; skip for quick runs.
## ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 90), "\n", sep = "")
cat("RTDP Benchmark — Part 4: Phase 2/3 quality on 2x2 R=10 (nSdx≈53k)\n")
cat("  [Exact DP ground truth required — takes ~2 min]\n")
cat(strrep("=", 90), "\n")

obj4  <- setup_env(2L, 2L, 4L, 1.0)
e4    <- obj4$env; tmpf4 <- obj4$tmpf

cat(sprintf("nSdx=%d  nAdx=%d  nScen=%d\n", e4$nSdx, e4$nAdx, e4$nScen))

cat("  Running exact DP... "); flush.console()
t0 <- proc.time()[["elapsed"]]
re4 <- rolling_dp_ptr(e4, tmpf4, numThreads = N_THREADS)
t_ex4 <- proc.time()[["elapsed"]] - t0
cat(sprintf("done (%.1f s)\n", t_ex4))

cat("  Running Phase 2/3 (n_traj=50, max_iter=10, tol=1e-3)... "); flush.console()
t0 <- proc.time()[["elapsed"]]
r2_4 <- rolling_dp_rtdp_p2(e4, tmpf4,
                             n_iter   = N_ITER,
                             n_traj   = N_TRAJ,
                             s0       = NULL,
                             epsilon  = 0.2,
                             tol      = 1e-3,
                             min_iter = 3L,
                             numThreads = N_THREADS, seed = SEED)
t_p2_4 <- proc.time()[["elapsed"]] - t0
cat(sprintf("done (%.1f s, n_done=%d)\n", t_p2_4, r2_4$n_done))

cat(sprintf("\n  Exact DP:   %.1f s\n", t_ex4))
cat(sprintf("  Phase 2/3:  %.1f s  (speedup: %.1fx)\n", t_p2_4, t_ex4 / t_p2_4))
cat(sprintf("  mean_err:   %.2f%%\n",
            mean_rel_err(r2_4$V_approx[seq(e4$tau),], re4$V[seq(e4$tau),], e4$tau)))
cat(sprintf("  max_err@t1: %.2f%%\n",
            max_rel_err_t1(r2_4$V_history_t1[r2_4$n_done,], re4$V[1L,])))
cat(sprintf("  pi@t1:      %.1f%%\n",
            pi_agree_t1(r2_4$pi_star, re4$pi_star)))
cat(sprintf("  cache@t=1:  %d/%d (%.1f%%)\n",
            r2_4$cache_coverage[r2_4$n_done], e4$nSdx,
            100 * r2_4$cache_coverage[r2_4$n_done] / e4$nSdx))

unlink(tmpf4); rm(e4, re4, r2_4); invisible(gc(verbose = FALSE))

cat("\n", strrep("=", 90), "\n", sep = "")
cat("Summary:\n")
cat("  Phase 1: iterative random K-sampling.  LP cost = K*tau*n_iter > exact DP → always slower.\n")
cat("  Phase 2: LP-accurate trajectory sampling.  LP cost independent of |S|.\n")
cat("           Scales to instances where exact DP is intractable (|S| > 50k).\n")
cat("  Phase 3: convergence stopping (tol) + ε-greedy exploration (epsilon).\n")
cat("           Early stopping saves iterations when VFA has stabilised.\n")
cat("           ε-greedy prevents trajectory collapse around a fixed start.\n")
cat(strrep("=", 90), "\n")
