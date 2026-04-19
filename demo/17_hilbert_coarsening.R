## ============================================================
## 17_hilbert_coarsening.R — Hilbert coarsening Pareto curve
##
## For coarsening level b = 0..4:
##   - subset = every 2^b-th state in Hilbert order (sampling rate 2^(-b))
##   - run rolling_dp_ptr restricted to that subset via stateSubset
##   - fill unsampled V values by nearest-neighbour in Hilbert order
##   - report speedup vs exact (b=0) and relative V error
##
## Produces a Pareto curve: error vs speedup as b grows.
## ============================================================

suppressPackageStartupMessages(library(TLPR))

N_THREADS <- max(1L, parallel::detectCores() - 2L)
CHUNK     <- 32L

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

# Nearest-neighbour fill in Hilbert order:
# hperm[k] = lex index of k-th state in Hilbert order.
# V[hperm] gives V in Hilbert order; NA positions get the value of the
# nearest non-NA neighbour in that ordering (forward pass then backward pass).
fill_nn_hilbert <- function(V, hperm) {
  Vh <- V[hperm]
  n  <- length(Vh)
  # Forward pass: propagate last non-NA to the right
  last_val <- NA_real_
  for (k in seq_len(n)) {
    if (!is.na(Vh[k])) { last_val <- Vh[k] } else { Vh[k] <- last_val }
  }
  # Backward pass: fill any leading NAs from the right
  for (k in rev(seq_len(n))) {
    if (!is.na(Vh[k])) { last_val <- Vh[k] } else { Vh[k] <- last_val }
  }
  V_out <- numeric(n)
  V_out[hperm] <- Vh
  V_out
}

# Run coarsened DP: every step-th state in Hilbert order is computed exactly;
# the rest are filled by nearest-neighbour before being fed as V_next.
run_coarsened <- function(e, tmpf, hperm, step, threads = N_THREADS) {
  nSdx <- e$nSdx
  nHil <- length(hperm)
  subset_idx <- hperm[seq(1L, nHil, by = step)]   # 1-indexed lex positions

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

  V      <- matrix(NA_real_, nrow = e$tau + 1L, ncol = nSdx)
  V[e$tau + 1L, ] <- V_term

  t0 <- proc.time()[["elapsed"]]
  for (t in seq(e$tau, 1L)) {
    # Use filled V_next so unsampled successors are valid
    V_next_full <- fill_nn_hilbert(V[t + 1L, ], hperm)
    res <- bellmanUpdatePtr(
      problem_ptr  = prob,
      t            = t - 1L,          # 0-indexed in C++
      stateSupport = as.double(e$stateSupport),
      flowSupport  = as.double(e$Q$vals),
      scnpb        = e$scnpb,
      alpha        = e$alpha,
      V_next       = V_next_full,
      numThreads     = threads,
      traversalOrder = "hilbert",
      chunkSize      = CHUNK,
      stateSubset    = subset_idx
    )
    V_t <- res$V_t
    # Fill unsampled states by NN before storing
    V[t, ] <- fill_nn_hilbert(V_t, hperm)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  list(V = V, elapsed = elapsed, n_solved = length(subset_idx))
}

## ============================================================
## Setup: 1×1 R=10 instance (same as demo/09)
## ============================================================

cat("\n", strrep("=", 65), "\n", sep = "")
cat("Hilbert coarsening Pareto curve — 1x1, R=10, tau=4\n")
cat(strrep("=", 65), "\n")

tmpf <- tempfile(fileext = ".json")
generate_instance(nI = 1L, nJ = 1L, tau = 4L,
                  nB = 5L, nCS = 10L, nCO = 1L, rate = 1.0,
                  seed = 42L, path = tmpf)

e <- new.env()
jsonlite::fromJSON(tmpf) |> list2env(envir = e)
e$from_i <- as_int_list(e$from_i)
e$to_j   <- as_int_list(e$to_j)
create_model(e)
dp_config(e)

# Hilbert order permutation (once per instance)
hperm <- hilbert_order(e$nI, e$nJ, e$R)
nSdx  <- e$nSdx
stopifnot(length(hperm) == nSdx, all(sort(hperm) == seq(nSdx)))

# Exact DP (b = 0, step = 1 → all states)
cat(sprintf("\nComputing exact DP (b=0, all %d states)...\n", nSdx))
res_exact <- run_coarsened(e, tmpf, hperm, step = 1L)
V_exact   <- res_exact$V
t_exact   <- res_exact$elapsed
cat(sprintf("  Exact time: %.2fs\n", t_exact))

## ============================================================
## Pareto sweep: b = 1..4
## ============================================================

fmt_hdr <- "%-6s  %8s  %8s  %10s  %10s  %10s\n"
fmt_row <- "%-6s  %8d  %8d  %10.2f  %10.2fx  %10.2f%%\n"

cat("\n")
cat(sprintf(fmt_hdr, "b", "|S_exact|", "|S_coarse|", "Time (s)", "Speedup", "Rel V-err (%)"))
cat(strrep("-", 65), "\n")

cat(sprintf(fmt_row, "0 (ref)", nSdx, nSdx, t_exact, 1.0, 0.0))

for (b in 1:4) {
  step <- 2L^b
  n_subset <- ceiling(nSdx / step)
  res  <- run_coarsened(e, tmpf, hperm, step = step)

  rel_err <- 100 * mean(abs(res$V - V_exact) / (abs(V_exact) + 1e-8), na.rm = TRUE)
  speedup  <- t_exact / max(res$elapsed, 1e-4)

  cat(sprintf(fmt_row, sprintf("b=%d", b), nSdx, n_subset,
              res$elapsed, speedup, rel_err))
}

cat(strrep("-", 65), "\n")
cat(sprintf("Threads: %d  |  Instance: 1x1 R=10 tau=4 seed=42\n", N_THREADS))
cat("Speedup > 1x: coarsening is faster. Rel V-err: mean |V_coarse-V_exact|/(|V_exact|+eps).\n")
cat("Nearest-neighbour fill in Hilbert order; no LP solved for unsampled states.\n")

unlink(tmpf); rm(e, V_exact); invisible(gc(verbose = FALSE))

## ============================================================
## Second instance: 2×2 R=5 (nSdx=4356 — demonstrates larger speedup)
## ============================================================

cat("\n", strrep("=", 65), "\n", sep = "")
cat("Hilbert coarsening Pareto curve — 2x2, R=5, tau=4\n")
cat(strrep("=", 65), "\n")

tmpf2 <- tempfile(fileext = ".json")
generate_instance(nI = 2L, nJ = 2L, tau = 4L,
                  nB = 5L, nCS = 10L, nCO = 1L, rate = 0.5,
                  seed = 42L, path = tmpf2)

e2 <- new.env()
jsonlite::fromJSON(tmpf2) |> list2env(envir = e2)
e2$from_i <- as_int_list(e2$from_i)
e2$to_j   <- as_int_list(e2$to_j)
create_model(e2)
dp_config(e2)

hperm2 <- hilbert_order(e2$nI, e2$nJ, e2$R)
nSdx2  <- e2$nSdx
stopifnot(length(hperm2) == nSdx2, all(sort(hperm2) == seq(nSdx2)))

cat(sprintf("\nComputing exact DP (b=0, all %d states)...\n", nSdx2))
res_exact2 <- run_coarsened(e2, tmpf2, hperm2, step = 1L)
V_exact2   <- res_exact2$V
t_exact2   <- res_exact2$elapsed
cat(sprintf("  Exact time: %.2fs\n", t_exact2))

cat("\n")
cat(sprintf(fmt_hdr, "b", "|S_exact|", "|S_coarse|", "Time (s)", "Speedup", "Rel V-err (%)"))
cat(strrep("-", 65), "\n")
cat(sprintf(fmt_row, "0 (ref)", nSdx2, nSdx2, t_exact2, 1.0, 0.0))

for (b in 1:4) {
  step <- 2L^b
  n_subset <- ceiling(nSdx2 / step)
  res <- run_coarsened(e2, tmpf2, hperm2, step = step)
  rel_err <- 100 * mean(abs(res$V - V_exact2) / (abs(V_exact2) + 1e-8), na.rm = TRUE)
  speedup  <- t_exact2 / max(res$elapsed, 1e-4)
  cat(sprintf(fmt_row, sprintf("b=%d", b), nSdx2, n_subset, res$elapsed, speedup, rel_err))
}

cat(strrep("-", 65), "\n")
cat(sprintf("Threads: %d  |  Instance: 2x2 R=5 tau=4 seed=42\n", N_THREADS))

unlink(tmpf2)
invisible(NULL)
