## ============================================================
## 09_hilbert_benchmark.R — Lexicographic vs Hilbert state traversal
##
## Part 1 — Correctness: both orders must produce identical V / Q.
## Part 2 — Timing A/B across several instance configs.
##
## Implementation note:
##   bellmanUpdateImpl calls Highs::setBasis() (no-arg) at the start of each
##   state to force a genuine cold start.  Without this, Highs::run() reuses
##   its internal optimal basis from the previous LP, causing traversal-order-
##   dependent tie-breaking on degenerate instances (multiple optimal routings).
##   The cold-start-per-state policy trades ~2× LP solve time for determinism:
##   both traversal orders produce bit-identical V and Q.
##   The residual speedup of Hilbert order is from CPU cache locality on the
##   constraint RHS data for geometrically adjacent states.
## ============================================================

suppressPackageStartupMessages(library(TLPR))

N_THREADS <- max(1L, parallel::detectCores() - 2L)
CHUNK     <- 32L

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

run_both <- function(tmpf, e, threads = N_THREADS, chunk = CHUNK) {
  prob <- loadProblemDataCx(tmpf)

  V_term <- -c(
    cbind(
      apply(e$Sdx[, e$I_,           drop = FALSE], 2L,
            function(sdx) e$stateSupport[sdx]),
      apply(e$Sdx[, e$nI + e$J_,   drop = FALSE], 2L,
            function(sdx) pmax(e$extendedStateSupport[sdx], 0L)),
     -apply(e$Sdx[, e$nI + e$J_,   drop = FALSE], 2L,
            function(sdx) pmin(e$extendedStateSupport[sdx], 0L))
    ) %*% c(e$alpha))

  t0   <- proc.time()[["elapsed"]]
  resL <- rolling_dp_ptr(e, tmpf, numThreads = threads,
                         traversalOrder = "lexicographic", chunkSize = chunk)
  t_lex <- proc.time()[["elapsed"]] - t0

  t0   <- proc.time()[["elapsed"]]
  resH <- rolling_dp_ptr(e, tmpf, numThreads = threads,
                         traversalOrder = "hilbert", chunkSize = chunk)
  t_hil <- proc.time()[["elapsed"]] - t0

  list(lex = resL, hil = resH, t_lex = t_lex, t_hil = t_hil)
}

## ============================================================
## Part 1 — Correctness (reference 1x1 R=10 instance)
## ============================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("Part 1 — Correctness check (1x1, R=10, tau=4)\n")
cat(strrep("=", 60), "\n")

tmpf1 <- tempfile(fileext = ".json")
generate_instance(nI = 1L, nJ = 1L, tau = 4L,
                  nB = 5L, nCS = 10L, nCO = 1L, rate = 1.0,
                  seed = 42L, path = tmpf1)

e1 <- new.env()
jsonlite::fromJSON(tmpf1) |> list2env(envir = e1)
e1$from_i <- as_int_list(e1$from_i)
e1$to_j   <- as_int_list(e1$to_j)
create_model(e1)
dp_config(e1)

res1 <- run_both(tmpf1, e1)

v_maxdiff <- max(abs(res1$lex$V - res1$hil$V), na.rm = TRUE)
q_maxdiff <- max(abs(res1$lex$Q - res1$hil$Q), na.rm = TRUE)
pi_agree  <- mean(res1$lex$pi_star == res1$hil$pi_star, na.rm = TRUE) * 100

cat(sprintf("  V  max |lex - hil| = %.2e  %s\n",
            v_maxdiff, if (v_maxdiff < 1e-8) "[PASS]" else "[FAIL]"))
cat(sprintf("  Q  max |lex - hil| = %.2e  %s\n",
            q_maxdiff, if (q_maxdiff < 1e-8) "[PASS]" else "[FAIL]"))
cat(sprintf("  pi_star agreement  = %.4f%%\n", pi_agree))
cat(sprintf("  t_lex = %.2fs  t_hil = %.2fs  ratio = %.3fx\n",
            res1$t_lex, res1$t_hil, res1$t_lex / res1$t_hil))

unlink(tmpf1); rm(e1, res1); invisible(gc(verbose = FALSE))

## ============================================================
## Part 2 — Timing benchmark
## ============================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Part 2 — Timing: Lexicographic vs Hilbert traversal order\n")
cat(strrep("=", 70), "\n")

MAX_NSDX <- 3000L   # skip configs whose state space exceeds this

configs <- list(
  list(nI=1L, nJ=1L, tau=4L,  rate=0.5,  label="1x1 R=5  tau=4"),
  list(nI=1L, nJ=1L, tau=4L,  rate=1.0,  label="1x1 R=10 tau=4"),
  list(nI=1L, nJ=1L, tau=4L,  rate=2.0,  label="1x1 R=20 tau=4"),
  list(nI=1L, nJ=1L, tau=8L,  rate=1.0,  label="1x1 R=10 tau=8"),
  list(nI=2L, nJ=1L, tau=4L,  rate=1.0,  label="2x1 R=10 tau=4"),
  list(nI=1L, nJ=2L, tau=4L,  rate=1.0,  label="1x2 R=10 tau=4"),
  list(nI=2L, nJ=2L, tau=4L,  rate=0.5,  label="2x2 R=5  tau=4"),
  list(nI=2L, nJ=2L, tau=4L,  rate=1.0,  label="2x2 R=10 tau=4"),
  list(nI=2L, nJ=2L, tau=1L,  rate=2.0,  label="2x2 R=20 tau=1")
)

fmt_hdr <- "%-22s  %8s  %5s  %6s  %9s  %9s  %8s  %8s\n"
fmt_row <- "%-22s  %8d  %5d  %6d  %9.2f  %9.2f  %8s  %8s\n"

cat(sprintf(fmt_hdr,
  "Config", "nSdx", "nAdx", "nScen",
  "Lex (s)", "Hil (s)", "Speedup", "V_diff"))
cat(strrep("-", 85), "\n")

for (cfg in configs) {
  R_cfg <- as.integer(round(10 * cfg$rate))
  nSdx_est <- (R_cfg + 1L)^cfg$nI * (2L * R_cfg + 1L)^cfg$nJ
  if (nSdx_est > MAX_NSDX) {
    cat(sprintf("%-22s  %8d  %5s  %6s  %9s  %9s  %8s  %8s\n",
                cfg$label, nSdx_est, "-", "-", "-", "-", "SKIP", "SKIP"))
    next
  }

  tmpf <- tempfile(fileext = ".json")
  generate_instance(
    nI = cfg$nI, nJ = cfg$nJ, tau = cfg$tau,
    nB = 5L, nCS = 10L, nCO = 1L, rate = cfg$rate,
    seed = 42L, path = tmpf)

  e <- new.env()
  jsonlite::fromJSON(tmpf) |> list2env(envir = e)
  e$from_i <- as_int_list(e$from_i)
  e$to_j   <- as_int_list(e$to_j)
  create_model(e)
  dp_config(e)

  nScen <- e$nQdx * e$nDdx * e$nWdx

  res <- run_both(tmpf, e)

  v_maxdiff <- max(abs(res$lex$V - res$hil$V), na.rm = TRUE)
  speedup <- res$t_lex / res$t_hil
  cat(sprintf(fmt_row,
    cfg$label, e$nSdx, e$nAdx, nScen,
    res$t_lex, res$t_hil,
    sprintf("%.3fx", speedup),
    if (v_maxdiff < 1e-8) "0 [OK]" else sprintf("%.1e!", v_maxdiff)))

  unlink(tmpf); rm(e, res); invisible(gc(verbose = FALSE))
}

cat(strrep("-", 85), "\n")
cat(sprintf("Threads: %d  |  Hilbert chunk size: %d\n", N_THREADS, CHUNK))
cat("Speedup > 1x: Hilbert is faster.  V_diff [OK]: identical DP solution.\n")
cat("Per-state Highs::setBasis() ensures cold start -> deterministic LP tie-breaking.\n")
