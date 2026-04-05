## ============================================================
## reproduce_paper.R — End-to-end numerical reproduction
##
## Sources scripts 01–03 in order and asserts all paper
## reference values from config/instance1x1_4.R.
##
## Usage:
##   Rscript demo/reproduce_paper.R
##   source("demo/reproduce_paper.R")   # from package root
##
## Optional environment variable:
##   TLPR_FULL=1   Run the full 1M Monte Carlo and EV/regret
##                 analysis (03_capacity_opt.R sections 5.2.3).
##                 Defaults to skipping those sections.
##
## Outputs checked against paper (arxiv 2505.01808):
##   v0      ≈ 557.2  (Table 1 / Section 5.2.2)
##   v_opt   ≈ 439.2  (Table 2)
##   reduction ≈ 21.2%
## ============================================================

suppressPackageStartupMessages(library(TLPR))

## -- Locate demo directory ------------------------------------
DEMO_DIR <- tryCatch(
  dirname(sys.frame(1L)$ofile),
  error = function(e) file.path(path.expand("~/drayage/TLPR"), "demo")
)

## -- Config --------------------------------------------------
source(file.path(DEMO_DIR, "config/instance1x1_4.R"))

## -- Step 1: Transit matrix ----------------------------------
cat("\n========  Step 1: Transit Matrix  ========\n")
source(file.path(DEMO_DIR, "01_transit.R"))

## -- Step 2: Dynamic programming ----------------------------
cat("\n========  Step 2: Dynamic Programming  ========\n")
source(file.path(DEMO_DIR, "02_dp.R"))

## -- Step 3: Capacity optimisation --------------------------
##    Only the fixed-scenario LP and L-BFGS-B optimisation
##    sections are required for the paper assertions.
##    The 1M Monte Carlo and EV/regret sections are optional.
cat("\n========  Step 3: Capacity Optimisation  ========\n")
source(file.path(DEMO_DIR, "03_capacity_opt.R"))

## -- Numerical assertions ------------------------------------
cat("\n========  Paper Assertions  ========\n")

## 1. Initial cost V1*(x0)
stopifnot(
  "V1*(x0) must be within 1.0 of paper value 557.2" =
    abs(v0 - PAPER_COST_INITIAL) < PAPER_ASSERT_TOL
)
cat(sprintf("  [PASS] V1*(x0)     = %7.2f   paper = %.1f\n", v0, PAPER_COST_INITIAL))

## 2. Optimal cost V1*(x*)
stopifnot(
  "V1*(x*) must be within 1.0 of paper value 439.2" =
    abs(v_opt - PAPER_COST_OPTIMAL) < PAPER_ASSERT_TOL
)
cat(sprintf("  [PASS] V1*(x*)     = %7.2f   paper = %.1f\n", v_opt, PAPER_COST_OPTIMAL))

## 3. Percentage cost reduction
reduction_pct <- (1.0 - v_opt / v0) * 100.0
stopifnot(
  "Reduction % must be within 0.5pp of paper value 21.2%" =
    abs(reduction_pct - PAPER_REDUCTION_PCT) < 0.5
)
cat(sprintf("  [PASS] Reduction   = %7.1f%%   paper = %.1f%%\n",
            reduction_pct, PAPER_REDUCTION_PCT))

cat("\nAll paper assertions PASSED.\n")
