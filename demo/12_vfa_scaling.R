## ============================================================
## 12_vfa_scaling.R — VFA execution time vs instance size
##
## Times rolling_dp_vfa (full-grid and K=30% sampled) across a
## range of instance sizes to find the 60-second frontier.
##
## Exact DP is timed where feasible (skipped once it exceeds
## EXACT_TIMEOUT seconds).  Traversal stops within each topology
## once full-grid VFA exceeds VFA_TIMEOUT seconds.
## ============================================================

suppressPackageStartupMessages(library(TLPR))
source("R/vfa.R")
source("R/dp_vfa.R")

N_THREADS   <- max(1L, parallel::detectCores() - 2L)
EXACT_TIMEOUT <- 60.0   # skip exact DP if prev run exceeded this
VFA_TIMEOUT   <- 90.0   # stop topology sweep if full-grid VFA exceeded this
K_FRAC        <- 0.30   # fraction of states for sampled VFA

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

run_config <- function(nI, nJ, tau, rate, skip_exact = FALSE) {
  tmpf <- tempfile(fileext = ".json")
  generate_instance(nI = nI, nJ = nJ, tau = tau,
                    nB = 5L, nCS = 10L, nCO = 1L, rate = rate,
                    seed = 42L, path = tmpf)
  e <- new.env()
  jsonlite::fromJSON(tmpf) |> list2env(envir = e)
  e$from_i <- as_int_list(e$from_i)
  e$to_j   <- as_int_list(e$to_j)
  create_model(e)
  dp_config(e)

  K <- max(1L, as.integer(e$nSdx * K_FRAC))

  t_exact <- NA_real_
  if (!skip_exact) {
    t0      <- proc.time()[["elapsed"]]
    rolling_dp_ptr(e, tmpf, numThreads = N_THREADS)
    t_exact <- proc.time()[["elapsed"]] - t0
  }

  t0     <- proc.time()[["elapsed"]]
  rolling_dp_vfa(e, tmpf, K = NULL, numThreads = N_THREADS)
  t_full <- proc.time()[["elapsed"]] - t0

  t0     <- proc.time()[["elapsed"]]
  rolling_dp_vfa(e, tmpf, K = K, numThreads = N_THREADS, seed = 1L)
  t_samp <- proc.time()[["elapsed"]] - t0

  unlink(tmpf)
  list(nSdx = e$nSdx, nAdx = e$nAdx, K = K,
       t_exact = t_exact, t_full = t_full, t_samp = t_samp)
}

## ── Configs: topology groups, increasing nSdx ─────────────────────────────────
topology_groups <- list(
  list(label = "1x1", nI = 1L, nJ = 1L, tau = 4L,
       rates = c(0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0)),

  list(label = "2x1", nI = 2L, nJ = 1L, tau = 4L,
       rates = c(0.5, 1.0, 1.5, 2.0, 2.5)),

  list(label = "1x2", nI = 1L, nJ = 2L, tau = 4L,
       rates = c(0.5, 1.0, 1.5, 2.0)),

  list(label = "2x2", nI = 2L, nJ = 2L, tau = 4L,
       rates = c(0.3, 0.5, 0.7, 1.0))
)

## ── Header ────────────────────────────────────────────────────────────────────
cat("\n", strrep("=", 95), "\n", sep = "")
cat(sprintf("VFA Scaling: execution time vs instance size  (K = %.0f%% of states, tau=4)\n",
            K_FRAC * 100))
cat(sprintf("Threads: %d  |  Exact DP skipped if prev > %.0fs  |  Sweep stops if VFA > %.0fs\n",
            N_THREADS, EXACT_TIMEOUT, VFA_TIMEOUT))
cat(strrep("=", 95), "\n\n")

fmt_hdr <- "%-6s  %5s  %5s  %6s  %7s  %10s  %10s  %10s  %9s  %9s\n"
fmt_row <- "%-6s  %5d  %5d  %6d  %7d  %10s  %10s  %10s  %9s  %9s\n"

cat(sprintf(fmt_hdr,
  "Topo", "R", "nSdx", "nAdx", "K (30%)",
  "t_exact", "t_vfa_full", "t_vfa_samp",
  "vs exact", "samp spdup"))
cat(strrep("-", 95), "\n")

for (grp in topology_groups) {
  cat(sprintf("\n  [%s]\n", grp$label))

  last_t_full  <- 0.0
  last_t_exact <- 0.0
  skip_exact   <- FALSE

  for (rate in grp$rates) {
    R_val <- as.integer(round(10 * rate))

    # Check if this group is already too slow
    if (last_t_full > VFA_TIMEOUT) {
      nSdx_est <- (R_val + 1L)^grp$nI * (2L * R_val + 1L)^grp$nJ
      cat(sprintf(fmt_row, grp$label, R_val, nSdx_est, R_val + 1L,
                  as.integer(nSdx_est * K_FRAC),
                  "—", "> timeout", "—", "—", "—"))
      next
    }

    if (last_t_exact > EXACT_TIMEOUT) skip_exact <- TRUE

    r <- run_config(grp$nI, grp$nJ, grp$tau, rate, skip_exact = skip_exact)

    exact_str <- if (is.na(r$t_exact)) "SKIP" else sprintf("%.1fs", r$t_exact)
    full_str  <- sprintf("%.1fs", r$t_full)
    samp_str  <- sprintf("%.1fs", r$t_samp)
    vs_exact  <- if (is.na(r$t_exact)) "—" else sprintf("%.2fx", r$t_exact / r$t_full)
    spdup     <- if (is.na(r$t_exact)) "—" else sprintf("%.2fx", r$t_exact / r$t_samp)

    # Mark if near or over 60s threshold
    mark_full <- if (r$t_full >= 60) " !" else ""
    mark_samp <- if (r$t_samp >= 60) " !" else ""

    cat(sprintf(fmt_row,
      grp$label, R_val, r$nSdx, r$nAdx, r$K,
      exact_str,
      paste0(full_str, mark_full),
      paste0(samp_str, mark_samp),
      vs_exact, spdup))

    last_t_full  <- r$t_full
    if (!is.na(r$t_exact)) last_t_exact <- r$t_exact

    invisible(gc(verbose = FALSE))
  }
}

cat("\n", strrep("-", 95), "\n", sep = "")
cat("t_vfa_full: VFA with all states solved per period (no state sampling).\n")
cat("t_vfa_samp: VFA with K=30% of states sampled per period.\n")
cat("vs exact: t_exact / t_vfa_full  (>1x means VFA overhead vs exact DP).\n")
cat("samp spdup: t_exact / t_vfa_samp  (speedup of sampled VFA over exact DP).\n")
cat("!: execution time >= 60s.\n")
