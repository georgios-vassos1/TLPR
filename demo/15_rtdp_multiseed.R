## ============================================================
## 15_rtdp_multiseed.R — Multi-seed RTDP Phase 2/3 confidence intervals
##
## Produces mean ± std for the scalability table (Table 3 of the paper):
##   - Wall-clock time (t_RTDP)
##   - Speedup over exact DP
##   - Mean relative error in V (all t)
##   - Max relative error in V at t=1
##   - Policy agreement at t=1
##
## Exact DP runs once per config (deterministic); RTDP reruns N_SEEDS times.
## Configs: 1x1 R=10, 2x1 R=5, 2x2 R=5, 2x2 R=10 (exact DP feasible for all).
## ============================================================

suppressPackageStartupMessages(library(TLPR))

N_THREADS <- max(1L, parallel::detectCores() - 2L)
N_ITER    <- 10L
N_TRAJ    <- 50L
K_FRAC    <- 0.30
N_SEEDS   <- 10L
SEEDS     <- seq_len(N_SEEDS) * 7L   # 7, 14, 21, ..., 70

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

setup_env <- function(nI, nJ, tau, rate, seed = 42L) {
  tmpf <- tempfile(fileext = ".json")
  generate_instance(nI = nI, nJ = nJ, tau = tau,
                    nB = 5L, nCS = 10L, nCO = 1L, rate = rate,
                    regime = "balanced", seed = seed, path = tmpf)
  e <- new.env()
  jsonlite::fromJSON(tmpf) |> list2env(envir = e)
  e$from_i <- as_int_list(e$from_i)
  e$to_j   <- as_int_list(e$to_j)
  create_model(e)
  dp_config(e)
  list(env = e, tmpf = tmpf)
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

max_rel_err_t1 <- function(v_hat, v_ex) {
  ok <- is.finite(v_ex) & is.finite(v_hat)
  if (!any(ok)) return(NA_real_)
  100 * max(abs(v_hat[ok] - v_ex[ok])) / max(abs(v_ex[ok]), 1e-8)
}

pi_agree_t1 <- function(pi_hat, pi_ex)
  mean(pi_hat[1L, ] == pi_ex[1L, ], na.rm = TRUE) * 100

fmt_msd <- function(m, s, fmt = "%.2f") {
  sprintf(paste0(fmt, " (", fmt, ")"), m, s)
}

## ─────────────────────────────────────────────────────────────────────────────
## Configs
## ─────────────────────────────────────────────────────────────────────────────

configs <- list(
  list(nI=1L, nJ=1L, tau=4L, rate=1.0, label="1x1 R=10"),
  list(nI=2L, nJ=1L, tau=4L, rate=0.5, label="2x1 R=5 "),
  list(nI=2L, nJ=2L, tau=4L, rate=0.5, label="2x2 R=5 "),
  list(nI=2L, nJ=2L, tau=4L, rate=1.0, label="2x2 R=10")
)

cat("\n", strrep("=", 90), "\n", sep = "")
cat(sprintf("RTDP Multi-seed Benchmark  |  N_SEEDS=%d  n_iter=%d  n_traj=%d  Threads=%d\n",
            N_SEEDS, N_ITER, N_TRAJ, N_THREADS))
cat(sprintf("Seeds: %s\n", paste(SEEDS, collapse = ", ")))
cat(strrep("=", 90), "\n\n")

results <- vector("list", length(configs))

for (ci in seq_along(configs)) {
  cfg <- configs[[ci]]
  cat(strrep("-", 90), "\n")
  cat(sprintf("Config: %s\n", cfg$label)); flush.console()

  ## Instance setup (fixed seed=42 so instance is the same across RTDP seeds)
  obj  <- setup_env(cfg$nI, cfg$nJ, cfg$tau, cfg$rate, seed = 42L)
  e    <- obj$env
  tmpf <- obj$tmpf
  K    <- max(1L, as.integer(e$nSdx * K_FRAC))

  cat(sprintf("  nSdx=%d  nAdx=%d  nScen=%d  K=%d\n",
              e$nSdx, e$nAdx, e$nScen, K))

  ## Exact DP — run once
  cat("  Running exact DP... "); flush.console()
  t0   <- proc.time()[["elapsed"]]
  re   <- rolling_dp_ptr(e, tmpf, numThreads = N_THREADS)
  t_ex <- proc.time()[["elapsed"]] - t0
  cat(sprintf("done (%.2f s)\n", t_ex)); flush.console()

  ## RTDP Phase 2 — loop over seeds
  metrics <- matrix(NA_real_, nrow = N_SEEDS,
                    ncol = 5L,
                    dimnames = list(NULL,
                      c("t_p2", "speedup", "mean_err", "max_err_t1", "pi_agree")))

  for (si in seq_len(N_SEEDS)) {
    cat(sprintf("  Seed %2d/%d ... ", SEEDS[si], max(SEEDS))); flush.console()
    t0  <- proc.time()[["elapsed"]]
    r2  <- rolling_dp_rtdp_p2(e, tmpf,
                               n_iter     = N_ITER,
                               n_traj     = N_TRAJ,
                               s0         = NULL,
                               epsilon    = 0.2,
                               tol        = 1e-3,
                               min_iter   = 3L,
                               numThreads = N_THREADS,
                               seed       = SEEDS[si])
    t_p2 <- proc.time()[["elapsed"]] - t0

    metrics[si, "t_p2"]       <- t_p2
    metrics[si, "speedup"]    <- t_ex / t_p2
    metrics[si, "mean_err"]   <- mean_rel_err(r2$V_approx[seq(e$tau), ],
                                               re$V[seq(e$tau), ], e$tau)
    metrics[si, "max_err_t1"] <- max_rel_err_t1(r2$V_history_t1[r2$n_done, ],
                                                  re$V[1L, ])
    metrics[si, "pi_agree"]   <- pi_agree_t1(r2$pi_star, re$pi_star)

    cat(sprintf("t=%.2fs  err=%.2f%%  pi=%.1f%%\n",
                t_p2, metrics[si, "mean_err"], metrics[si, "pi_agree"]))
    flush.console()
    rm(r2); invisible(gc(verbose = FALSE))
  }

  results[[ci]] <- list(
    label   = cfg$label,
    nSdx    = e$nSdx,
    t_ex    = t_ex,
    metrics = metrics
  )

  unlink(tmpf); rm(e, re); invisible(gc(verbose = FALSE))
}

## ─────────────────────────────────────────────────────────────────────────────
## Summary table
## ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 90), "\n", sep = "")
cat(sprintf("Summary: mean (std) across %d seeds\n", N_SEEDS))
cat(strrep("=", 90), "\n\n")

hdr <- sprintf("%-10s  %7s  %8s  %16s  %16s  %16s  %14s",
               "Config", "|S|", "t_exact",
               "t_RTDP [s]", "Speedup", "Mean err [%]", "pi@t=1 [%]")
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")

for (res in results) {
  m <- colMeans(res$metrics, na.rm = TRUE)
  s <- apply(res$metrics, 2L, sd, na.rm = TRUE)

  cat(sprintf("%-10s  %7d  %8.2f  %16s  %16s  %16s  %14s\n",
              res$label,
              res$nSdx,
              res$t_ex,
              fmt_msd(m["t_p2"],       s["t_p2"],       "%.2f"),
              fmt_msd(m["speedup"],    s["speedup"],    "%.1f"),
              fmt_msd(m["mean_err"],   s["mean_err"],   "%.2f"),
              fmt_msd(m["pi_agree"],   s["pi_agree"],   "%.1f")))
}

cat(strrep("-", nchar(hdr)), "\n")
cat(sprintf("n_iter=%d  |  n_traj=%d  |  K=%.0f%%  |  Threads=%d  |  N_seeds=%d\n",
            N_ITER, N_TRAJ, K_FRAC * 100, N_THREADS, N_SEEDS))

## ─────────────────────────────────────────────────────────────────────────────
## LaTeX table snippet (paste into paper)
## ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 90), "\n", sep = "")
cat("LaTeX table rows (mean (std) format):\n")
cat(strrep("=", 90), "\n")

for (res in results) {
  m <- colMeans(res$metrics, na.rm = TRUE)
  s <- apply(res$metrics, 2L, sd, na.rm = TRUE)

  cat(sprintf(
    "%s & %s & %.2f & %s & %s & %s & %s \\\\\n",
    gsub(" ", "", res$label),
    formatC(res$nSdx, format = "d", big.mark = ","),
    res$t_ex,
    sprintf("%.2f (%.2f)", m["t_p2"],     s["t_p2"]),
    sprintf("%.1f (%.1f)", m["speedup"],  s["speedup"]),
    sprintf("%.2f (%.2f)", m["mean_err"], s["mean_err"]),
    sprintf("%.1f (%.1f)", m["pi_agree"], s["pi_agree"])
  ))
}
