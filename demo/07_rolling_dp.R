## ============================================================
## 07_rolling_dp.R — Standard vs rolling Bellman update
##
## Compares two paths to the DP solution:
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
##                             normalises by P(feasible | state, action)
##
## Note on infeasible scenarios:
##   When the LP has no solution for a (state, action, scenario) triple
##   (typically: action > S.I + Q), the two paths handle it differently:
##     Path A — excluded from the sum; remaining probs do NOT sum to 1
##     Path B — excluded; remaining probs renormalised to sum to 1
##   This causes V to differ whenever infeasible triples exist.
##   Both are valid approximations; path B gives E[cost | feasible].
## ============================================================

suppressPackageStartupMessages(library(TLPR))

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

cat(sprintf("  V difference: path B %.4f%% %s path A (renormalisation effect)\n",
            abs(dpB$V[1L, i0] / dpA$V[1L, i0] - 1) * 100,
            if (dpB$V[1L, i0] > dpA$V[1L, i0]) "higher than" else "lower than"))

## ── Timing / memory summary ──────────────────────────────────
cat("\n[07] Timing and memory:\n")
cat(sprintf("  Path A  dp=%.2fs  (transit %.1f MB already in memory)\n",
            t_A, object.size(transit) / 1e6))
cat(sprintf("  Path B  rolling=%.2fs  working memory ~%.0f doubles (%.1f KB)\n",
            t_B, env$nSdx * env$nAdx * 4.0,
            env$nSdx * env$nAdx * 4.0 * 8 / 1024))
cat(sprintf("  Path B is %.1fx %s for the DP step alone\n",
            max(t_A, t_B) / min(t_A, t_B),
            if (t_B < t_A) "faster" else "slower (includes LP solves)"))
cat(sprintf("  For large instances where transit cannot fit in RAM, path B is the only option.\n"))
