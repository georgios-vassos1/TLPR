## ============================================================
## 16_head_to_head.R — Validate mstp_to_tlpr_json() adapter
##
## Tests MSTP::mstp_to_tlpr_json() via a round-trip:
##   1. Generate a TLPR instance (source of truth)
##   2. Convert the TLPR env to an MSTP-style inst list
##   3. Call mstp_to_tlpr_json() → TLPR JSON
##   4. Run exact DP and RTDP on the adapter-produced JSON
##   5. Compare against exact DP on the original JSON
##
## Instance: 1×1, R=5, tau=4 (nSdx=66; exact DP < 1 s)
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(MSTP)
})

N_THREADS <- max(1L, parallel::detectCores() - 2L)

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

# ── Helper: build an MSTP-style inst from a loaded TLPR env ───────────────────
# This mirrors the fields expected by MSTP::mstp_to_tlpr_json().
tlpr_env_to_mstp_inst <- function(e) {
  list(
    tau             = e$tau,
    nOrigins        = e$nI,
    nDestinations   = e$nJ,
    nCarriers       = e$nCS,
    nSpotCarriers   = e$nCO,
    # Bids and Winners: use the TLPR auction structure directly
    Bids            = e$B,
    Winners         = setNames(e$winner, NULL),  # unnamed list of bid-index vecs
    Ldx             = unlist(e$L_, use.names = FALSE),
    nLc             = e$nLc,
    # Capacities (same matrices as in the TLPR JSON)
    entry_capacity  = rep(e$R, e$nI),
    exit_capacity   = rep(e$R, e$nJ),
    # R stores matrices column-major; MSTP carrier_capacity is a flat column-major
    # vector: position k = Cb[(k-1)%%tau + 1, (k-1)%/%tau + 1]
    carrier_capacity = c(as.vector(e$Cb),   # nCS*tau entries (strategic)
                         as.vector(e$Co)),  # nCO*tau entries (spot)
    # Cost coefficients
    entry_store_coef = e$alpha[seq(e$nI)],
    exit_store_coef  = e$alpha[e$nI + seq(e$nJ)],
    exit_short_coef  = e$alpha[e$nI + e$nJ + seq(e$nJ)],
    transport_coef   = unlist(e$CTb_list, use.names = FALSE),
    # Spot rates: CTo is tau × (nCO*nL), stored column-major → flat vector
    spot_coef        = as.vector(e$CTo),
    # Initial stocks (zero for DP comparison)
    entry_stock_0    = rep(0L, e$nI),
    exit_stock_0     = rep(0L, e$nJ),
    exit_short_0     = rep(0L, e$nJ)
  )
}

cat("\n", strrep("=", 72), "\n", sep = "")
cat("16_head_to_head.R — mstp_to_tlpr_json() round-trip validation\n")
cat(strrep("=", 72), "\n\n")

# ── 1. Generate TLPR instance (source of truth) ────────────────────────────────
cat("--- Step 1: TLPR generate_instance (1×1, R=5, tau=4) ---\n")
tmpf_orig <- tempfile(fileext = ".json")
TLPR::generate_instance(
  nI = 1L, nJ = 1L, tau = 4L,
  nB = 5L, nCS = 5L, nCO = 1L, rate = 0.5,
  seed = 7L, path = tmpf_orig
)

e_orig <- new.env()
jsonlite::fromJSON(tmpf_orig) |> list2env(envir = e_orig)
e_orig$from_i <- as_int_list(e_orig$from_i)
e_orig$to_j   <- as_int_list(e_orig$to_j)
dp_config(e_orig)
cat(sprintf("  nI=%d  nJ=%d  R=%d  tau=%d  nSdx=%d  nCS=%d  nCO=%d  nBids=%d\n",
            e_orig$nI, e_orig$nJ, e_orig$R, e_orig$tau,
            e_orig$nSdx, e_orig$nCS, e_orig$nCO, length(e_orig$B)))

# ── 2. Convert TLPR env → MSTP inst → TLPR JSON ───────────────────────────────
cat("\n--- Step 2: convert env → MSTP inst → TLPR JSON (via adapter) ---\n")
inst    <- tlpr_env_to_mstp_inst(e_orig)
tmpf_rt <- tempfile(fileext = ".json")

MSTP::mstp_to_tlpr_json(
  inst = inst,
  R    = e_orig$R,
  Q    = e_orig$Q,
  D    = e_orig$D,   # pass demand distribution explicitly (may differ from Q)
  W    = e_orig$W,
  path = tmpf_rt
)

e_rt <- new.env()
jsonlite::fromJSON(tmpf_rt) |> list2env(envir = e_rt)
e_rt$from_i <- as_int_list(e_rt$from_i)
e_rt$to_j   <- as_int_list(e_rt$to_j)
dp_config(e_rt)
cat(sprintf("  Adapter env: nSdx=%d  nAdx=%d  nScen=%d\n",
            e_rt$nSdx, e_rt$nAdx, e_rt$nScen))

# ── 3. Sanity checks on reconstructed env ─────────────────────────────────────
cat("\n--- Step 3: field-level sanity checks ---\n")
checks <- list(
  R         = all(e_orig$R       == e_rt$R),
  tau       = all(e_orig$tau     == e_rt$tau),
  nSdx      = all(e_orig$nSdx    == e_rt$nSdx),
  nScen     = all(e_orig$nScen   == e_rt$nScen),
  stateKeys = all(e_orig$stateKeys == e_rt$stateKeys),
  flowKeys  = all(e_orig$flowKeys  == e_rt$flowKeys),
  scnpb_sum  = isTRUE(all.equal(sum(e_orig$scnpb), sum(e_rt$scnpb))),
  scnpb_vals = isTRUE(all.equal(e_orig$scnpb,      e_rt$scnpb))
)
for (nm in names(checks)) {
  cat(sprintf("  %-14s  %s\n", nm, if (checks[[nm]]) "OK" else "MISMATCH"))
}

# ── 4. Exact DP on original vs adapter JSON ────────────────────────────────────
cat("\n--- Step 4: exact DP comparison ---\n")

t0        <- proc.time()[["elapsed"]]
res_orig  <- rolling_dp_ptr(e_orig, tmpf_orig, numThreads = N_THREADS)
t_orig    <- proc.time()[["elapsed"]] - t0

t0        <- proc.time()[["elapsed"]]
res_rt    <- rolling_dp_ptr(e_rt,   tmpf_rt,   numThreads = N_THREADS)
t_rt      <- proc.time()[["elapsed"]] - t0

# V comparison at t=1
ok <- is.finite(res_orig$V[1L, ]) & is.finite(res_rt$V[1L, ])
max_ref   <- max(abs(res_orig$V[1L, ok]))
rel_err_V <- if (any(ok) && max_ref > 0) {
  100 * max(abs(res_rt$V[1L, ok] - res_orig$V[1L, ok])) / max_ref
} else NA_real_

pi_agr <- mean(res_orig$pi_star[1L, ] == res_rt$pi_star[1L, ], na.rm = TRUE) * 100

cat(sprintf("  Original JSON:   t=%.3fs  V[t=1] in [%.2f, %.2f]\n",
            t_orig,
            min(res_orig$V[1L,], na.rm=TRUE), max(res_orig$V[1L,], na.rm=TRUE)))
cat(sprintf("  Adapter JSON:    t=%.3fs  V[t=1] in [%.2f, %.2f]\n",
            t_rt,
            min(res_rt$V[1L,], na.rm=TRUE), max(res_rt$V[1L,], na.rm=TRUE)))
cat(sprintf("  max rel |V-err| at t=1: %.4f%%\n", rel_err_V))
cat(sprintf("  pi_star agreement t=1:  %.1f%%\n", pi_agr))

# ── 5. RTDP Phase 2 on adapter JSON ───────────────────────────────────────────
cat("\n--- Step 5: RTDP Phase 2 on adapter JSON (n_iter=10, n_traj=50) ---\n")
t0 <- proc.time()[["elapsed"]]
res_rtdp <- rolling_dp_rtdp_p2(
  e_rt, tmpf_rt,
  n_iter     = 10L,
  n_traj     = 50L,
  epsilon    = 0.2,
  tol        = 1e-3,
  numThreads = N_THREADS,
  seed       = 7L
)
t_rtdp <- proc.time()[["elapsed"]] - t0

ok_r <- is.finite(res_orig$V[1L,]) & is.finite(res_rtdp$V_approx[1L,])
rel_rtdp <- if (any(ok_r) && max_ref > 0) {
  100 * mean(abs(res_rtdp$V_approx[1L, ok_r] - res_orig$V[1L, ok_r])) / max_ref
} else NA_real_
pi_rtdp <- mean(res_orig$pi_star[1L,] == res_rtdp$pi_star[1L,], na.rm=TRUE) * 100

cat(sprintf("  t_rtdp = %.3fs  |  speedup = %.2fx\n", t_rtdp, t_orig / t_rtdp))
cat(sprintf("  mean rel |V-err| at t=1: %.2f%%  |  pi_agr: %.1f%%\n", rel_rtdp, pi_rtdp))

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n", strrep("-", 72), "\n", sep = "")
cat(sprintf("%-26s  %7s  %7s  %9s  %8s\n",
            "Method", "t (s)", "speedup", "rel err%", "pi agr%"))
cat(strrep("-", 72), "\n")
cat(sprintf("%-26s  %7.3f  %7s  %9s  %8s\n",
            "Exact DP (original JSON)", t_orig, "1.00x", "0.000%", "100.0%"))
cat(sprintf("%-26s  %7.3f  %7.2fx  %9.4f%%  %8.1f%%\n",
            "Exact DP (adapter JSON)", t_rt, t_orig / t_rt, rel_err_V, pi_agr))
cat(sprintf("%-26s  %7.3f  %7.2fx  %9.2f%%  %8.1f%%\n",
            "RTDP-P2 (adapter JSON)", t_rtdp, t_orig / t_rtdp, rel_rtdp, pi_rtdp))
cat(strrep("-", 72), "\n")
cat("Round-trip: TLPR env → MSTP inst → mstp_to_tlpr_json() → TLPR env.\n")
cat("Exact DP adapter error should be ~0 (same instance, same JSON structure).\n")

unlink(c(tmpf_orig, tmpf_rt))
invisible(gc(verbose = FALSE))
