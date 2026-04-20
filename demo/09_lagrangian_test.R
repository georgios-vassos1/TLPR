## ============================================================
## 09_lagrangian_test.R — Lagrangian heuristic vs HiGHS LP
##
## Section A: regression on the 5 canonical instances (JSON files)
## Section B: generated instances covering a wider range of
##   topologies, regimes, and LP sizes
##
## For each instance the script randomly samples LP calls across
## (t, state, inflow, spot-rate, action) and compares:
##   optimizeModelFromJSON  — exact LP (HiGHS)
##   solveLPHeuristicCx     — Lagrangian dual ascent + greedy
##   solveLPGreedyCx        — pure greedy (no dual ascent)
##
## Outputs per-instance and global summary statistics:
##   - LP variable count
##   - feasibility rate
##   - relative objective error vs HiGHS
##   - mean dual gap and mean iteration count
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(jsonlite)
})

## ── locate demo dir ────────────────────────────────────────────────────────
.demo_dir <- local({
  ofiles <- Filter(Negate(is.null), lapply(sys.frames(), `[[`, "ofile"))
  if (length(ofiles)) dirname(normalizePath(tail(ofiles, 1L)[[1L]])) else {
    farg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(farg)) dirname(normalizePath(sub("--file=", "", farg[1L]))) else getwd()
  }
})

TLPR_ROOT <- path.expand("~/drayage/TLPR")
INST_DIR  <- file.path(TLPR_ROOT, "src/instances")

## Tuning knobs
N_CASES   <- 500L   # random LP calls per instance
LAGR_ITER <- 50L    # Lagrangian dual-ascent iterations
set.seed(42L)

## ── helper: run one instance from its JSON path ───────────────────────────
run_instance <- function(json_path, label = basename(json_path)) {
  env_raw <- jsonlite::fromJSON(json_path)

  R_     <- env_raw$R
  nI     <- env_raw$nI
  nJ     <- env_raw$nJ
  nQ     <- env_raw$nQ
  nW     <- env_raw$nW
  nCO    <- env_raw$nCO
  nL     <- env_raw$nL
  tau    <- env_raw$tau
  n_vars <- env_raw$nL_ + nCO * nL

  ss     <- seq(0L, R_)
  es     <- c(-rev(ss[-1L]), ss)
  Q_vals <- env_raw$Q$vals
  W_vals <- env_raw$W$vals

  t_vec      <- sample(0L:(tau - 1L),   N_CASES, replace = TRUE)
  action_vec <- sample(0L:R_,           N_CASES, replace = TRUE)
  k1_vec     <- sample(0L:(nQ - 1L),   N_CASES, replace = TRUE)
  k3_vec     <- sample(0L:(nW - 1L),   N_CASES, replace = TRUE)
  orig_idx_mat <- matrix(sample(0L:R_,         N_CASES * nI, replace = TRUE), N_CASES, nI)
  dest_idx_mat <- matrix(sample(0L:(2L * R_),  N_CASES * nJ, replace = TRUE), N_CASES, nJ)

  rows    <- vector("list", N_CASES)
  n_valid <- 0L

  for (rep in seq_len(N_CASES)) {
    t      <- t_vec[rep]
    action <- action_vec[rep]
    k1     <- k1_vec[rep]
    k3     <- k3_vec[rep]
    s_orig <- ss[orig_idx_mat[rep, ] + 1L]
    s_dest <- es[dest_idx_mat[rep, ] + 1L]

    limits <- as.integer(c(s_orig + Q_vals[k1 + 1L], R_ - s_dest))
    sp     <- rep(W_vals[k3 + 1L], nCO * nL)

    ref <- tryCatch(
      optimizeModelFromJSON(json_path, t, sp, limits, action),
      error = function(e) list(status = "ERROR")
    )
    if (ref$status != "OPTIMAL") next

    heur <- tryCatch(
      solveLPHeuristicCx(json_path, t, sp, limits, action, nIter = LAGR_ITER),
      error = function(e) list(objval = NA_real_, feasible = FALSE,
                               dual_gap = NA_real_, n_iter = NA_integer_)
    )
    grdy <- tryCatch(
      solveLPGreedyCx(json_path, t, sp, limits, action),
      error = function(e) list(objval = NA_real_, feasible = FALSE)
    )

    ref_val  <- ref$objval
    abs_base <- abs(ref_val) + 1e-8
    n_valid  <- n_valid + 1L
    rows[[n_valid]] <- data.frame(
      instance  = label,
      t         = t,
      action    = action,
      n_vars    = n_vars,
      ref_val   = ref_val,
      heur_val  = ifelse(heur$feasible, heur$objval, NA_real_),
      grdy_val  = ifelse(grdy$feasible, grdy$objval, NA_real_),
      heur_feas = heur$feasible,
      grdy_feas = grdy$feasible,
      heur_gap  = ifelse(heur$feasible, heur$dual_gap, NA_real_),
      heur_iter = heur$n_iter,
      heur_rerr = ifelse(heur$feasible, abs(heur$objval - ref_val) / abs_base, NA_real_),
      grdy_rerr = ifelse(grdy$feasible, abs(grdy$objval - ref_val) / abs_base, NA_real_),
      stringsAsFactors = FALSE
    )
  }

  if (n_valid == 0L) {
    cat(sprintf("\n[%s]  No feasible LPs sampled — skipping.\n", label))
    return(NULL)
  }

  df <- do.call(rbind, rows[seq_len(n_valid)])

  cat(sprintf("\n[%s]\n", label))
  cat(sprintf("  Sampled / feasible   : %d / %d\n", N_CASES, n_valid))
  cat(sprintf("  n=vars               : %d\n", n_vars))

  for (method in c("Lagrangian", "Greedy")) {
    feas_col  <- if (method == "Lagrangian") "heur_feas" else "grdy_feas"
    rerr_col  <- if (method == "Lagrangian") "heur_rerr" else "grdy_rerr"
    feas_rate <- mean(df[[feas_col]], na.rm = TRUE)
    df_feas   <- df[!is.na(df[[rerr_col]]), ]
    if (nrow(df_feas) > 0L) {
      cat(sprintf("  %-10s feasible : %5.1f%%  mean rerr : %7.4f%%  max rerr : %7.4f%%\n",
                  method, 100 * feas_rate,
                  100 * mean(df_feas[[rerr_col]]),
                  100 * max(df_feas[[rerr_col]])))
    } else {
      cat(sprintf("  %-10s feasible : %5.1f%%  (no feasible rows)\n",
                  method, 100 * feas_rate))
    }
  }
  df_hf <- df[df$heur_feas, ]
  if (nrow(df_hf) > 0L) {
    cat(sprintf("  Mean dual gap        : %.4f\n", mean(df_hf$heur_gap, na.rm = TRUE)))
    cat(sprintf("  Mean Lagr iterations : %.1f / %d\n",
                mean(df_hf$heur_iter, na.rm = TRUE), LAGR_ITER))
  }
  df
}

## ============================================================
## Section A: canonical instances
## ============================================================
cat("\n=== SECTION A: CANONICAL INSTANCES ===\n")

canonical_files <- c(
  "instance1x1_4_001.json",
  "instance1x1_12_001.json",
  "instance1x2_4_001.json",
  "instance2x1_4_001.json",
  "instance2x2_4_001.json"
)

results_A <- lapply(canonical_files, function(f) {
  run_instance(file.path(INST_DIR, f), label = f)
})

## ============================================================
## Section B: generated instances
## ============================================================
cat("\n=== SECTION B: GENERATED INSTANCES ===\n")

## Descriptor list: each entry passed to generate_instance() + a label.
## nCO == nCS throughout (spot and strategic carrier counts are equal by design).
gen_specs <- list(
  ## ── small-medium (n ~ 20-60) ──────────────────────────────────────────────
  list(label = "2x2_bal_nCS3",   nI=2L, nJ=2L, tau=6L,  nB=6L,  nCS=3L, nCO=3L, rate=1.0, regime="balanced",  seed=101L),
  list(label = "2x3_tight_nCS3", nI=2L, nJ=3L, tau=4L,  nB=8L,  nCS=3L, nCO=3L, rate=1.0, regime="tight",     seed=102L),
  list(label = "3x2_vol_nCS3",   nI=3L, nJ=2L, tau=6L,  nB=8L,  nCS=3L, nCO=3L, rate=1.0, regime="volatile",  seed=103L),
  list(label = "3x3_bal_nCS4",   nI=3L, nJ=3L, tau=4L,  nB=10L, nCS=4L, nCO=4L, rate=1.0, regime="balanced",  seed=201L),
  list(label = "3x3_tight_nCS4", nI=3L, nJ=3L, tau=6L,  nB=10L, nCS=4L, nCO=4L, rate=1.0, regime="tight",     seed=202L),
  list(label = "3x3_vol_nCS4",   nI=3L, nJ=3L, tau=6L,  nB=10L, nCS=4L, nCO=4L, rate=1.0, regime="volatile",  seed=203L),
  ## ── medium-large (n ~ 100-250) ────────────────────────────────────────────
  list(label = "4x3_tight_nCS5", nI=4L, nJ=3L, tau=6L,  nB=15L, nCS=5L, nCO=5L, rate=1.0, regime="tight",     seed=301L),
  list(label = "3x4_vol_nCS5",   nI=3L, nJ=4L, tau=6L,  nB=15L, nCS=5L, nCO=5L, rate=1.0, regime="volatile",  seed=302L),
  list(label = "4x4_bal_nCS6",   nI=4L, nJ=4L, tau=6L,  nB=20L, nCS=6L, nCO=6L, rate=1.0, regime="balanced",  seed=401L),
  list(label = "4x4_tight_nCS6", nI=4L, nJ=4L, tau=4L,  nB=20L, nCS=6L, nCO=6L, rate=1.0, regime="tight",     seed=402L),
  list(label = "4x4_vol_nCS6",   nI=4L, nJ=4L, tau=8L,  nB=20L, nCS=6L, nCO=6L, rate=1.0, regime="volatile",  seed=403L),
  ## ── large (n ~ 300-600) ───────────────────────────────────────────────────
  list(label = "5x4_tight_nCS7", nI=5L, nJ=4L, tau=6L,  nB=25L, nCS=7L, nCO=7L, rate=1.0, regime="tight",     seed=501L),
  list(label = "4x5_bal_nCS7",   nI=4L, nJ=5L, tau=6L,  nB=25L, nCS=7L, nCO=7L, rate=1.0, regime="balanced",  seed=502L),
  list(label = "5x5_bal_nCS8",   nI=5L, nJ=5L, tau=4L,  nB=30L, nCS=8L, nCO=8L, rate=1.0, regime="balanced",  seed=601L),
  list(label = "5x5_tight_nCS8", nI=5L, nJ=5L, tau=6L,  nB=30L, nCS=8L, nCO=8L, rate=1.0, regime="tight",     seed=602L),
  list(label = "5x5_vol_nCS8",   nI=5L, nJ=5L, tau=6L,  nB=35L, nCS=8L, nCO=8L, rate=1.0, regime="volatile",  seed=603L)
)

results_B <- lapply(gen_specs, function(spec) {
  ## Write instance to a temp file so run_instance can use the standard JSON path
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  do.call(generate_instance, c(spec[setdiff(names(spec), "label")], list(path = tmp)))
  run_instance(tmp, label = spec$label)
})

## ============================================================
## Global summary (both sections)
## ============================================================
cat("\n=== GLOBAL SUMMARY ===\n")
df_all <- do.call(rbind, Filter(Negate(is.null), c(results_A, results_B)))

if (is.null(df_all) || nrow(df_all) == 0L) {
  cat("No results to summarise.\n")
} else {
  cat(sprintf("Total feasible LPs       : %d\n", nrow(df_all)))
  cat(sprintf("Lagrangian feasibility   : %.1f%%\n", 100 * mean(df_all$heur_feas, na.rm = TRUE)))
  cat(sprintf("Greedy    feasibility    : %.1f%%\n", 100 * mean(df_all$grdy_feas, na.rm = TRUE)))

  df_hf <- df_all[df_all$heur_feas & !is.na(df_all$heur_rerr), ]
  df_gf <- df_all[df_all$grdy_feas & !is.na(df_all$grdy_rerr), ]
  if (nrow(df_hf) > 0L) {
    cat(sprintf("Lagrangian mean rel err  : %.4f%%\n", 100 * mean(df_hf$heur_rerr)))
    cat(sprintf("Lagrangian max rel err   : %.4f%%\n", 100 * max(df_hf$heur_rerr)))
  }
  if (nrow(df_gf) > 0L) {
    cat(sprintf("Greedy    mean rel err   : %.4f%%\n", 100 * mean(df_gf$grdy_rerr)))
    cat(sprintf("Greedy    max rel err    : %.4f%%\n", 100 * max(df_gf$grdy_rerr)))
  }

  ## Per-n_vars bucket summary
  cat("\nSummary by LP size (n_vars):\n")
  breaks <- c(0L, 5L, 10L, 20L, 50L, 100L, Inf)
  labels <- c("≤5", "6-10", "11-20", "21-50", "51-100", ">100")
  df_hf$bucket <- cut(df_hf$n_vars, breaks = breaks, labels = labels, right = TRUE)
  df_gf$bucket <- cut(df_gf$n_vars, breaks = breaks, labels = labels, right = TRUE)
  for (b in labels) {
    h <- df_hf[df_hf$bucket == b, ]
    g <- df_gf[df_gf$bucket == b, ]
    if (nrow(h) == 0L) next
    cat(sprintf("  n=%s  n_lp=%d  lagr %.4f%%  greedy %.4f%%\n",
                b, nrow(h),
                100 * mean(h$heur_rerr),
                if (nrow(g) > 0L) 100 * mean(g$grdy_rerr) else NA_real_))
  }

  cat("\nWorst Lagrangian cases (by relative error):\n")
  df_hf_sorted <- df_hf[order(df_hf$heur_rerr, decreasing = TRUE), ]
  top_n <- min(10L, nrow(df_hf_sorted))
  print(df_hf_sorted[seq_len(top_n),
                      c("instance", "n_vars", "t", "action", "ref_val",
                        "heur_val", "heur_rerr", "heur_gap", "heur_iter")],
        row.names = FALSE, digits = 4L)
}

cat("\n[09] Lagrangian test complete.\n")
