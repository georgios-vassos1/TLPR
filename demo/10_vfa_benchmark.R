## ============================================================
## 10_vfa_benchmark.R — Separable VFA vs Exact DP
##
## Part 1 — Approximation quality on 1x1 R=10 tau=4
##           (exact DP feasible; use as ground truth)
## Part 2 — Error and speedup across instance sizes
##           Compares exact DP, full-grid VFA, and sampled VFA (K states)
## ============================================================

suppressPackageStartupMessages(library(TLPR))

N_THREADS <- max(1L, parallel::detectCores() - 2L)

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

setup_env <- function(nI, nJ, tau, rate, seed = 42L) {
  tmpf <- tempfile(fileext = ".json")
  generate_instance(nI = nI, nJ = nJ, tau = tau,
                    nB = 5L, nCS = 10L, nCO = 1L, rate = rate,
                    seed = seed, path = tmpf)
  e <- new.env()
  jsonlite::fromJSON(tmpf) |> list2env(envir = e)
  e$from_i <- as_int_list(e$from_i)
  e$to_j   <- as_int_list(e$to_j)
  create_model(e)
  dp_config(e)
  list(env = e, tmpf = tmpf)
}

## ============================================================
## Part 1 — Approximation quality (1x1 R=10 tau=4)
## ============================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("Part 1 — VFA vs Exact DP (1x1, R=10, tau=4)\n")
cat(strrep("=", 60), "\n")

cfg1  <- setup_env(nI = 1L, nJ = 1L, tau = 4L, rate = 1.0)
e1    <- cfg1$env
tmpf1 <- cfg1$tmpf

t0        <- proc.time()[["elapsed"]]
res_exact <- rolling_dp_ptr(e1, tmpf1, numThreads = N_THREADS)
t_exact   <- proc.time()[["elapsed"]] - t0

# Full-grid VFA (approximate continuation, all states solved)
t0      <- proc.time()[["elapsed"]]
res_vfa <- rolling_dp_vfa(e1, tmpf1, K = NULL, numThreads = N_THREADS)
t_full  <- proc.time()[["elapsed"]] - t0

# Sampled VFA: K = 50% of states
K_half  <- max(1L, as.integer(e1$nSdx * 0.5))
t0      <- proc.time()[["elapsed"]]
res_smp <- rolling_dp_vfa(e1, tmpf1, K = K_half, numThreads = N_THREADS, seed = 1L)
t_samp  <- proc.time()[["elapsed"]] - t0

cat(sprintf("  nSdx = %d   K (50%%) = %d\n", e1$nSdx, K_half))
cat(sprintf("  %-12s  %9s  %9s  %9s\n", "Period", "Exact", "VFA-full", "VFA-samp"))
for (t in seq(e1$tau + 1L, 1L)) {
  v_ex <- res_exact$V[t, ]
  rel_full <- if (all(is.na(res_vfa$V_approx[t, ]))) NA_real_ else {
    ok <- is.finite(v_ex) & is.finite(res_vfa$V_approx[t, ])
    100 * max(abs(res_vfa$V_approx[t, ok] - v_ex[ok])) / max(abs(v_ex[ok]))
  }
  rel_samp <- if (all(is.na(res_smp$V_approx[t, ]))) NA_real_ else {
    ok <- is.finite(v_ex) & is.finite(res_smp$V_approx[t, ])
    100 * max(abs(res_smp$V_approx[t, ok] - v_ex[ok])) / max(abs(v_ex[ok]))
  }
  lbl <- if (t == e1$tau + 1L) "terminal   " else sprintf("t=%d        ", t)
  cat(sprintf("  %s  %8s  %8s  %8s\n", lbl,
              "0.000%",
              if (is.na(rel_full)) "N/A" else sprintf("%.3f%%", rel_full),
              if (is.na(rel_samp)) "N/A" else sprintf("%.3f%%", rel_samp)))
}

pi_full <- mean(res_exact$pi_star[1L,] == res_vfa$pi_star[1L,],  na.rm = TRUE) * 100
pi_samp <- mean(res_exact$pi_star[1L,] == res_smp$pi_star[1L,],  na.rm = TRUE) * 100
cat(sprintf("  pi_star agreement t=1: VFA-full = %.2f%%  VFA-samp = %.2f%%\n",
            pi_full, pi_samp))
cat(sprintf("  t_exact = %.2fs  t_full = %.2fs  t_samp = %.2fs  speedup = %.2fx\n",
            t_exact, t_full, t_samp, t_exact / t_samp))

unlink(tmpf1); rm(e1, res_exact, res_vfa, res_smp); invisible(gc(verbose = FALSE))

## ============================================================
## Part 2 — Error and speedup across instance sizes
## ============================================================

cat("\n", strrep("=", 80), "\n", sep = "")
cat("Part 2 — VFA error and speedup vs instance size  (K = 30% of states)\n")
cat(strrep("=", 80), "\n")

configs <- list(
  list(nI=1L, nJ=1L, tau=4L, rate=0.5, label="1x1 R=5  tau=4"),
  list(nI=1L, nJ=1L, tau=4L, rate=1.0, label="1x1 R=10 tau=4"),
  list(nI=2L, nJ=1L, tau=4L, rate=1.0, label="2x1 R=10 tau=4"),
  list(nI=1L, nJ=2L, tau=4L, rate=1.0, label="1x2 R=10 tau=4"),
  list(nI=2L, nJ=2L, tau=4L, rate=0.5, label="2x2 R=5  tau=4")
)

fmt_hdr <- "%-20s  %7s  %5s  %10s  %10s  %10s  %8s\n"
fmt_row <- "%-20s  %7d  %5d  %10.3f%%  %10.3f%%  %10.2fs  %8.2fx\n"
cat(sprintf(fmt_hdr, "Config", "nSdx", "K", "full rel%", "samp rel%",
            "t_exact", "speedup"))
cat(strrep("-", 82), "\n")

for (cfg in configs) {
  obj   <- setup_env(cfg$nI, cfg$nJ, cfg$tau, cfg$rate)
  e     <- obj$env
  tmpf  <- obj$tmpf
  K_cfg <- max(1L, as.integer(e$nSdx * 0.3))

  t0 <- proc.time()[["elapsed"]]
  re <- rolling_dp_ptr(e, tmpf, numThreads = N_THREADS)
  t_ex <- proc.time()[["elapsed"]] - t0

  rv <- rolling_dp_vfa(e, tmpf, K = NULL,  numThreads = N_THREADS)
  t0 <- proc.time()[["elapsed"]]
  rs <- rolling_dp_vfa(e, tmpf, K = K_cfg, numThreads = N_THREADS, seed = 1L)
  t_sp <- proc.time()[["elapsed"]] - t0

  v_ex <- re$V[1L, ]
  ok_f <- is.finite(v_ex) & is.finite(rv$V_approx[1L, ])
  ok_s <- is.finite(v_ex) & is.finite(rs$V_approx[1L, ])
  rel_f <- if (any(ok_f)) 100*max(abs(rv$V_approx[1L,ok_f]-v_ex[ok_f]))/max(abs(v_ex[ok_f])) else NA_real_
  rel_s <- if (any(ok_s)) 100*max(abs(rs$V_approx[1L,ok_s]-v_ex[ok_s]))/max(abs(v_ex[ok_s])) else NA_real_

  cat(sprintf(fmt_row, cfg$label, e$nSdx, K_cfg,
              rel_f, rel_s, t_ex, t_ex / t_sp))

  unlink(tmpf); rm(e, re, rv, rs); invisible(gc(verbose = FALSE))
}

cat(strrep("-", 82), "\n")
cat(sprintf("Threads: %d\n", N_THREADS))
cat("full rel%%: error with all states solved but approximate continuation.\n")
cat("samp rel%%: error with K sampled states — genuine speedup over exact DP.\n")
cat("speedup: t_exact / t_sampled_vfa.\n")
