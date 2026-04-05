## ============================================================
## Shared configuration for the 1x1, tau=4 instance
## All demo scripts source this file first.
## ============================================================

## Package root — only path to update if the repo moves
TLPR_ROOT <- path.expand("~/drayage/TLPR")

## Instance
JSON_PATH <- file.path(TLPR_ROOT, "src/instances/instance1x1_4_001.json")

## Fixed scenario used in the paper (Table 1):
##   t=1: q=8, d=8, w=7
##   t=2: q=8, d=8, w=22
##   t=3: q=0, d=8, w=7
##   t=4: q=0, d=0, w=22
FIXED_VARPHIDX <- c(9L, 18L, 7L, 10L)

## Alpha scaling applied in all experiments
ALPHA_SCALE <- 3.0

## Reservation costs (Table 2, strategic carrier, per-period)
RESERVATION_COSTS <- c(8.52, 4.46, 4.25, 9.70)

## Initial state for simulation
S0 <- c(0L, 8L)

## Cache files (written to tempdir; delete to force recomputation)
TRANSIT_CACHE <- file.path(tempdir(), "tlpr_transit_1x1_4.rds")
DP_CACHE      <- file.path(tempdir(), "tlpr_dp_1x1_4.rds")

## Paper reference values for numerical assertions
PAPER_COST_INITIAL  <- 557.2   # V1*(x0) — Table 1 / Section 5.2.2
PAPER_COST_OPTIMAL  <- 439.2   # V1*(x*) — Table 2
PAPER_REDUCTION_PCT <- 21.2    # percentage cost reduction
PAPER_ASSERT_TOL    <- 1.0     # acceptable absolute tolerance

cat(sprintf("[config] JSON: %s\n", JSON_PATH))
cat(sprintf("[config] Fixed scenario varphidx: %s\n", paste(FIXED_VARPHIDX, collapse = ", ")))
cat(sprintf("[config] Cache: %s\n", TRANSIT_CACHE))
