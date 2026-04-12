## ============================================================
## 07_rolling_dp.R — Standard vs rolling Bellman update
##
## Part 1 — Correctness cross-check on the reference 1x1 instance:
##   Compares two paths to the DP solution:
##
##   Path A (standard):
##     computeEnvironmentCx  — materialise full transit table
##                             (tau * nSdx * nAdx * nScen x 6 rows)
##     dynamic_programming   — R-level backward induction over table;
##                             sums raw scnpb[feasible] (not renormalised)
##
##   Path B (rolling, C++):
##     rolling_dp_cx         — one bellmanUpdateCx call per period;
##                             no transit table; O(nSdx * nAdx) memory;
##                             sums raw scnpb[feasible] identically to Path A
##
##   Both paths now aggregate identically: infeasible (s,a,omega) triples
##   are excluded; remaining probabilities are NOT renormalised.
##   V values and pi_star should agree exactly on the same instance.
##
## Part 2 — Scalability benchmark across the same config grid as 06:
##   For each config, times Path A (Cx + DP) and Path B (rolling).
##   Configs where the transit table exceeds MAX_TRANSIT_ROWS are
##   skipped for Path A but Path B continues (memory-safe).
## ============================================================

suppressPackageStartupMessages(library(TLPR))

## ============================================================
## Part 1 — Correctness cross-check (reference 1x1 instance)
## ============================================================

## -- Config / dependencies -----------------------------------
.demo_dir <- local({
  ofiles <- Filter(Negate(is.null), lapply(sys.frames(), `[[`, "ofile"))
  if (length(ofiles)) dirname(normalizePath(tail(ofiles, 1L)[[1L]])) else {
    farg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(farg)) dirname(normalizePath(sub("--file=", "", farg[1L]))) else getwd()
  }
})
if (!exists("JSON_PATH")) source(file.path(.demo_dir, "config/instance1x1_4.R"))

## 01_transit.R sets up env (alpha scaling, from_i/to_j conversion,
## dp_config) and builds / loads the full transit table.
if (!exists("transit")) source(file.path(.demo_dir, "01_transit.R"))

i0 <- which(
  env$stateSupport[env$Sdx[, 1L]] == S0[1L] &
  env$extendedStateSupport[env$Sdx[, 2L]] == S0[2L])

infeas_frac <- mean(is.na(transit[, 1L]))
cat(sprintf("\n[07] Instance: |S|=%d  |A|=%d  |Omega|=%d  tau=%d\n",
            env$nSdx, env$nAdx, env$nScen, env$tau))
cat(sprintf("[07] Transit table: %s rows x 6 cols  (%.1f MB)  infeasible: %.1f%%\n",
            format(nrow(transit), big.mark = ","),
            object.size(transit) / 1e6,
            infeas_frac * 100))

## ── Path A: dynamic_programming (R, full table) ─────────────
cat("\n[A] Running dynamic_programming (R, full transit table)...\n")
t0  <- proc.time()
dpA <- dynamic_programming(env, transit)
t_A <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("[A] Done in %.2f s\n", t_A))
cat(sprintf("[A] V[t=1, S0=(0,8)] = %.4f\n", dpA$V[1L, i0]))

## ── Path B: rolling_dp_cx (C++, no transit table) ───────────
cat("\n[B] Running rolling_dp_cx (C++, no transit table)...\n")
t0  <- proc.time()
dpB <- rolling_dp_cx(env, JSON_PATH, numThreads = 8L)
t_B <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("[B] Done in %.2f s\n", t_B))
cat(sprintf("[B] V[t=1, S0=(0,8)] = %.4f\n", dpB$V[1L, i0]))

## ── Comparison ───────────────────────────────────────────────
cat("\n[07] Comparison:\n")

# Policy agreement (where both have a defined greedy action)
both_def <- !is.na(dpA$pi_star[1L, ]) & !is.na(dpB$pi_star[1L, ])
pi_agree <- if (any(both_def)) {
  mean(dpA$pi_star[1L, both_def] == dpB$pi_star[1L, both_def]) * 100
} else NA_real_
cat(sprintf("  pi_star agreement at t=1: %.1f%% of %d defined states\n",
            pi_agree, sum(both_def)))

# V rank correlation across states (do both rank states similarly?)
both_v <- is.finite(dpA$V[1L, ]) & is.finite(dpB$V[1L, ])
rho <- cor(dpA$V[1L, both_v], dpB$V[1L, both_v], method = "spearman")
cat(sprintf("  V[1, ] Spearman rank corr: %.4f\n", rho))

v_diff_pct <- abs(dpB$V[1L, i0] / dpA$V[1L, i0] - 1) * 100
cat(sprintf("  V[S0] relative diff: %.6f%%  (0.000000%% => exact agreement)\n", v_diff_pct))

## ── Timing / memory summary ──────────────────────────────────
cat("\n[07] Timing and memory (reference instance):\n")
cat(sprintf("  Path A  dp=%.2fs  transit=%.1f MB\n",
            t_A, object.size(transit) / 1e6))
cat(sprintf("  Path B  rolling=%.2fs  working memory ~%.0f doubles (%.1f KB)\n",
            t_B, env$nSdx * env$nAdx * 4.0,
            env$nSdx * env$nAdx * 4.0 * 8 / 1024))
cat(sprintf("  Path B total (LP + DP): %.1fx %s\n",
            max(t_A, t_B) / min(t_A, t_B),
            if (t_B < t_A) "faster" else "slower (includes LP solves)"))

rm(dpA, dpB)
invisible(gc(verbose = FALSE))

## ============================================================
## Part 2 — Scalability benchmark across topology / state space
## ============================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Part 2 — Scalability: Path A (Cx+DP) vs Path B (rolling_dp_cx)\n")
cat(strrep("=", 70), "\n", sep = "")

MAX_TRANSIT_ROWS <- 25e6L
N_THREADS        <- max(1L, parallel::detectCores() - 2L)

configs <- list(
  ## A: State space — vary rate at fixed 1x1
  list(nI=1L, nJ=1L, tau=4L, rate=0.5, group="A: State space (1x1)"),
  list(nI=1L, nJ=1L, tau=4L, rate=1.0, group="A: State space (1x1)"),
  list(nI=1L, nJ=1L, tau=4L, rate=2.0, group="A: State space (1x1)"),
  list(nI=1L, nJ=1L, tau=4L, rate=4.0, group="A: State space (1x1)"),

  ## B: Time horizon — vary tau at fixed 1x1, rate=1
  list(nI=1L, nJ=1L, tau= 4L, rate=1.0, group="B: Time horizon (1x1, R=10)"),
  list(nI=1L, nJ=1L, tau= 8L, rate=1.0, group="B: Time horizon (1x1, R=10)"),
  list(nI=1L, nJ=1L, tau=12L, rate=1.0, group="B: Time horizon (1x1, R=10)"),

  ## C: Topology — curse of dimensionality
  list(nI=1L, nJ=1L, tau=4L, rate=1.0, group="C: Topology (tau=4, R=10)"),
  list(nI=2L, nJ=1L, tau=4L, rate=1.0, group="C: Topology (tau=4, R=10)"),
  list(nI=1L, nJ=2L, tau=4L, rate=1.0, group="C: Topology (tau=4, R=10)"),
  list(nI=2L, nJ=2L, tau=4L, rate=0.5, group="C: Topology (tau=4, R=5)"),
  list(nI=2L, nJ=2L, tau=4L, rate=1.0, group="C: Topology (tau=4, R=10)")
)

## Column widths:
##  label(22)  nSdx(7)  nAdx(5)  nScen(6)  transit_rows(12)
##  Cx_s(8)  DP_s(7)  Roll_s(8)  TransMB(9)  WorkKB(8)  Speedup(8)
fmt_hdr <- "%-22s  %7s  %5s  %6s  %12s  %8s  %7s  %8s  %9s  %8s  %8s\n"
fmt_row <- "%-22s  %7d  %5d  %6d  %12s  %8s  %7s  %8.2f  %9.1f  %8.1f  %8s\n"

cat(sprintf(fmt_hdr,
  "Config", "nSdx", "nAdx", "nScen", "Transit rows",
  "Cx (s)", "DP (s)", "Roll (s)", "Trans (MB)", "Work (KB)", "Speedup"))
cat(strrep("-", 110), "\n")

last_group <- ""

for (cfg in configs) {
  nI   <- cfg$nI
  nJ   <- cfg$nJ
  tau  <- cfg$tau
  rate <- cfg$rate
  grp  <- cfg$group

  if (grp != last_group) { cat(sprintf("\n  [%s]\n", grp)); last_group <- grp }

  R     <- as.integer(round(10 * rate))
  nSdx  <- (R + 1L)^nI * (2L * R + 1L)^nJ
  nAdx  <- R + 1L
  nScen <- 3L^nI * 3L^nJ * 3L   # nQ=3, nD=3, nW=3, nCO=1 (balanced regime)
  transit_rows <- tau * nSdx * nAdx * nScen

  label    <- sprintf("%dx%d tau=%d R=%d", nI, nJ, tau, R)
  trans_mb <- transit_rows * 6 * 8 / 1024^2
  work_kb  <- nSdx * nAdx * 4.0 * 8 / 1024   # 4 doubles per (s,a): sumCost,sumV,sumW,Q

  ## Generate instance and write to temp JSON
  tmpf <- tempfile(fileext = ".json")
  generate_instance(
    nI = nI, nJ = nJ, tau = tau,
    nB = 5L, nCS = 10L, nCO = 1L, rate = rate,
    seed = 42L, path = tmpf)

  ## Load and configure
  e <- new.env()
  jsonlite::fromJSON(tmpf) |> list2env(envir = e)
  as_int_list <- function(x)
    if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)
  e$from_i <- as_int_list(e$from_i)
  e$to_j   <- as_int_list(e$to_j)
  create_model(e)
  dp_config(e)

  ## -- Path A (Cx + DP) — skipped if transit table too large ----
  cx_str <- "SKIP"
  dp_str <- "SKIP"
  if (transit_rows <= MAX_TRANSIT_ROWS) {
    tbl <- matrix(NA_real_, nrow = transit_rows, ncol = 6L)
    t0  <- proc.time()[["elapsed"]]
    for (t in seq(tau)) {
      res <- computeEnvironmentCx(
        tmpf, t - 1L,
        seq(0L, e$R),
        e$Q$vals,
        numThreads = N_THREADS)
      offset <- (t - 1L) * e$nSdx * e$nAdx * e$nScen
      tbl[offset + seq(nrow(res)), ] <- res
    }
    cx_s  <- proc.time()[["elapsed"]] - t0
    t0    <- proc.time()[["elapsed"]]
    dpOut <- dynamic_programming(e, tbl)
    dp_s  <- proc.time()[["elapsed"]] - t0
    cx_str <- sprintf("%.2f", cx_s)
    dp_str <- sprintf("%.2f", dp_s)
    rm(tbl, dpOut)
    invisible(gc(verbose = FALSE))
  }

  ## -- Path B (rolling_dp_cx) -----------------------------------
  t0     <- proc.time()[["elapsed"]]
  dpRoll <- rolling_dp_cx(e, tmpf, numThreads = N_THREADS)
  roll_s <- proc.time()[["elapsed"]] - t0
  rm(dpRoll)
  invisible(gc(verbose = FALSE))

  ## Speedup = (Cx + DP) / rolling  (only when Path A ran)
  speedup_str <- if (cx_str != "SKIP") {
    ratio <- (cx_s + dp_s) / roll_s
    sprintf("%.1fx", ratio)
  } else {
    "N/A (A skip)"
  }

  trans_str <- if (transit_rows > MAX_TRANSIT_ROWS) {
    sprintf("%.1fM SKIP", transit_rows / 1e6)
  } else {
    format(transit_rows, big.mark = ",")
  }

  cat(sprintf(fmt_row,
    label, nSdx, nAdx, nScen, trans_str,
    cx_str, dp_str, roll_s,
    trans_mb, work_kb, speedup_str))

  unlink(tmpf)
  rm(e)
  invisible(gc(verbose = FALSE))
}

cat("\n", strrep("-", 110), "\n", sep = "")
cat(sprintf("Threads: %d  |  Transit row limit: %s\n",
  N_THREADS, format(MAX_TRANSIT_ROWS, big.mark = ",")))
cat("Speedup = (Cx time + DP time) / rolling time\n")
cat("Work KB = O(nSdx * nAdx) working memory used by rolling_dp_cx\n")
cat("Trans MB = hypothetical transit table size (6 doubles * rows)\n")
