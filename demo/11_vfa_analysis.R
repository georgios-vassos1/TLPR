## ============================================================
## 11_vfa_analysis.R — VFA approximation quality analysis
##
## Investigates when and why the separable VFA succeeds or fails.
## Uses multiple seeds and reports mean/median/max relative error,
## R², and policy agreement — across topologies and regimes.
##
## Configs (fast exact DP, < 2s each):
##   1x1 R=10  (nSdx=231)
##   2x1 R=5   (nSdx=396)
##   1x2 R=5   (nSdx=726)
## Regimes: balanced, tight, loose  (3 seeds each)
## Total: 3 × 3 × 3 = 27 exact+VFA pairs  (~75s)
## ============================================================

suppressPackageStartupMessages(library(TLPR))
source("R/vfa.R")
source("R/dp_vfa.R")

N_THREADS <- max(1L, parallel::detectCores() - 2L)
SEEDS     <- c(1L, 42L, 99L)
REGIMES   <- c("balanced", "tight", "loose")

as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

setup_env <- function(nI, nJ, tau, rate, regime, seed) {
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

# Aggregate relative errors and R² across all (period, state) pairs.
# Denominator per period: max|V_exact| (consistent across periods).
vfa_metrics <- function(re, rv, e) {
  all_err <- numeric(0L)
  r2_vals <- numeric(e$tau)

  for (t in seq(e$tau)) {
    v_ex  <- re$V[t, ]
    v_hat <- rv$V_approx[t, ]
    ok    <- is.finite(v_ex) & is.finite(v_hat)
    if (!any(ok)) { r2_vals[t] <- NA_real_; next }

    denom      <- max(abs(v_ex[ok]))
    all_err    <- c(all_err, abs(v_hat[ok] - v_ex[ok]) / max(denom, 1e-8) * 100)
    ss_res     <- sum((v_hat[ok] - v_ex[ok])^2)
    ss_tot     <- sum((v_ex[ok] - mean(v_ex[ok]))^2)
    r2_vals[t] <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  }

  list(
    err_mean    = mean(all_err,   na.rm = TRUE),
    err_median  = median(all_err, na.rm = TRUE),
    err_max     = max(all_err,    na.rm = TRUE),
    r2          = mean(r2_vals,   na.rm = TRUE),
    pi_all      = mean(re$pi_star == rv$pi_star, na.rm = TRUE) * 100,
    pi_t1       = mean(re$pi_star[1L, ] == rv$pi_star[1L, ], na.rm = TRUE) * 100
  )
}

configs <- list(
  list(nI=1L, nJ=1L, tau=4L, rate=1.0, label="1x1 R=10"),
  list(nI=2L, nJ=1L, tau=4L, rate=0.5, label="2x1 R=5 "),
  list(nI=1L, nJ=2L, tau=4L, rate=0.5, label="1x2 R=5 ")
  ## Optional — add ~72s: uncomment if a larger 1x1 case is needed
  # list(nI=1L, nJ=1L, tau=4L, rate=2.0, label="1x1 R=20")
)

## ── Header ────────────────────────────────────────────────────────────────────
cat("\n", strrep("=", 90), "\n", sep = "")
cat("VFA Approximation Quality Analysis  (full-grid VFA vs exact DP)\n")
cat(sprintf("Topologies: %d  |  Regimes: %d  |  Seeds: %d  |  Threads: %d\n",
            length(configs), length(REGIMES), length(SEEDS), N_THREADS))
cat(strrep("=", 90), "\n\n")

fmt_hdr <- "%-10s  %-9s  %5s  %9s  %11s  %9s  %5s  %9s  %9s\n"
fmt_row <- "%-10s  %-9s  %5d  %9.2f%%  %11.2f%%  %9.2f%%  %5.2f  %9.1f%%  %9.1f%%\n"
fmt_sep <- "%-10s  %s\n"

cat(sprintf(fmt_hdr,
            "Config", "Regime", "nSdx",
            "mean err%", "median err%", "max err%", "R²",
            "π agr all", "π agr t=1"))
cat(strrep("-", 90), "\n")

## ── Main sweep ────────────────────────────────────────────────────────────────
for (cfg in configs) {
  R_val <- as.integer(round(10 * cfg$rate))
  nSdx  <- (R_val + 1L)^cfg$nI * (2L * R_val + 1L)^cfg$nJ

  cat(sprintf(fmt_sep, cfg$label, sprintf("[nSdx = %d]", nSdx)))

  for (regime in REGIMES) {
    # Collect metrics across seeds
    m <- lapply(SEEDS, function(seed) {
      obj <- setup_env(cfg$nI, cfg$nJ, cfg$tau, cfg$rate, regime, seed)
      e   <- obj$env; tmpf <- obj$tmpf
      re  <- rolling_dp_ptr(e, tmpf, numThreads = N_THREADS)
      rv  <- rolling_dp_vfa(e, tmpf, K = NULL, numThreads = N_THREADS)
      out <- vfa_metrics(re, rv, e)
      unlink(tmpf); rm(e, re, rv); invisible(gc(verbose = FALSE))
      out
    })

    # Average across seeds
    avg <- function(key) mean(sapply(m, `[[`, key), na.rm = TRUE)

    cat(sprintf(fmt_row,
                "", regime, nSdx,
                avg("err_mean"), avg("err_median"), avg("err_max"),
                avg("r2"), avg("pi_all"), avg("pi_t1")))
  }
  cat(strrep("-", 90), "\n")
}

cat(sprintf("Seeds: %s\n", paste(SEEDS, collapse = ", ")))
cat("mean/median/max: |V_approx - V_exact| / max|V_exact| per period, pooled across t=1..tau.\n")
cat("R²: goodness of separable fit averaged across periods. 1.0 = perfect additive structure.\n")
cat("π agr all: policy agreement across all (t, s). π agr t=1: agreement at the first decision.\n")
