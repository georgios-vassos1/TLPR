## ============================================================
## 08_xptr_equivalence.R — XPtr API equivalence tests
##
## Verifies that the three new XPtr-based functions produce output
## numerically identical to their JSON-file counterparts:
##
##   loadProblemDataCx + computeEnvironmentPtr  ==  computeEnvironmentCx
##   loadProblemDataCx + bellmanUpdatePtr        ==  bellmanUpdateCx
##   createHIGHSmodel + createTransportVarsCx +
##     addCapacity/StorageLimit/VolumeConstraintCx
##     + solveModelCx                            ==  optimizeModelFromJSON
##
## Also cross-checks rolling_dp_ptr against rolling_dp_cx (both should
## produce the same V / Q / pi_star / pi_rand matrices).
## ============================================================

suppressPackageStartupMessages(library(TLPR))

## -- Config --------------------------------------------------
.demo_dir <- local({
  ofiles <- Filter(Negate(is.null), lapply(sys.frames(), `[[`, "ofile"))
  if (length(ofiles)) dirname(normalizePath(tail(ofiles, 1L)[[1L]])) else {
    farg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(farg)) dirname(normalizePath(sub("--file=", "", farg[1L]))) else getwd()
  }
})
if (!exists("JSON_PATH")) source(file.path(.demo_dir, "config/instance1x1_4.R"))

env_raw <- jsonlite::fromJSON(JSON_PATH)

## dp_config initialises env with all index arrays needed by rolling_dp_*
if (!exists("env") || !exists("JSON_PATH")) source(file.path(.demo_dir, "01_transit.R"))

## ============================================================
## 1. computeEnvironmentPtr  vs  computeEnvironmentCx
## ============================================================
cat("\n[08] Test 1: computeEnvironmentPtr vs computeEnvironmentCx\n")

t_test  <- 0L
ss      <- seq(0L, env_raw$R)
fs      <- env_raw$Q$vals
nT      <- 2L  # check first 2 periods only (it's slow)

prob <- loadProblemDataCx(JSON_PATH)

for (t in seq(0L, nT - 1L)) {
  ref <- computeEnvironmentCx(JSON_PATH, t, ss, fs, numThreads = 4L)
  ptr <- computeEnvironmentPtr(prob,      t, ss, fs, numThreads = 4L)

  stopifnot(
    "dim mismatch"          = identical(dim(ref), dim(ptr)),
    "values not identical"  = isTRUE(all.equal(ref, ptr, tolerance = 1e-12))
  )
  cat(sprintf("  t=%d  OK  (dim %d x %d, max |diff| = %.2e)\n",
              t, nrow(ref), ncol(ref), max(abs(ref - ptr), na.rm = TRUE)))
}

cat("[08] Test 1 PASSED\n")

## ============================================================
## 2. bellmanUpdatePtr  vs  bellmanUpdateCx  (one period)
## ============================================================
cat("\n[08] Test 2: bellmanUpdatePtr vs bellmanUpdateCx\n")

## Terminal value function (same as rolling_dp_cx)
V_terminal <- -c(
  cbind(
    apply(env$Sdx[, env$I_,          drop = FALSE], 2L,
          function(sdx) env$stateSupport[sdx]),
    apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
          function(sdx) pmax(env$extendedStateSupport[sdx], 0L)),
   -apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
          function(sdx) pmin(env$extendedStateSupport[sdx], 0L))
  ) %*% c(env$alpha))

for (t in seq(env$tau, max(env$tau - 1L, 1L))) {
  ref <- bellmanUpdateCx(
    jsonFile     = JSON_PATH,
    t            = t - 1L,
    stateSupport = as.double(env$stateSupport),
    flowSupport  = as.double(env$Q$vals),
    scnpb        = env$scnpb,
    alpha        = env$alpha,
    V_next       = V_terminal,
    numThreads   = 4L
  )
  ptr <- bellmanUpdatePtr(
    problem_ptr  = prob,
    t            = t - 1L,
    stateSupport = as.double(env$stateSupport),
    flowSupport  = as.double(env$Q$vals),
    scnpb        = env$scnpb,
    alpha        = env$alpha,
    V_next       = V_terminal,
    numThreads   = 4L
  )

  for (field in c("V_t", "Q_t", "pi_star_t", "pi_rand_t")) {
    stopifnot(isTRUE(all.equal(ref[[field]], ptr[[field]], tolerance = 1e-12)))
  }
  cat(sprintf("  t=%d  OK\n", t))
}

cat("[08] Test 2 PASSED\n")

## ============================================================
## 3. solveModelCx  vs  optimizeModelFromJSON
## ============================================================
cat("\n[08] Test 3: solveModelCx vs optimizeModelFromJSON\n")

t_lp    <- 0L
limits  <- rep(env_raw$R, env_raw$nI + env_raw$nJ)
sp      <- env_raw$CTo[1L, ]

## jsonlite returns B and L as matrices; the XPtr chain expects lists of
## integer vectors (same representation the C++ code sees from JSON).
bids_list  <- lapply(seq_len(nrow(env_raw$B)), function(i) as.integer(env_raw$B[i, ]))
lanes_list <- lapply(seq_len(nrow(env_raw$L)), function(i) as.integer(env_raw$L[i, ]))

ref_lp <- optimizeModelFromJSON(JSON_PATH, t_lp, sp, limits, env_raw$R)

## Build the same model via the XPtr chain
model_ptr     <- createHIGHSmodel()
n_vars        <- env_raw$nL_ + env_raw$nCO * env_raw$nL
transport_ptr <- createTransportVarsCx(
  model_ptr     = model_ptr,
  n             = n_vars,
  winnerKeys    = env_raw$winnerKey,
  winners       = env_raw$winner,
  bids          = bids_list,
  lanes         = lanes_list,
  contractRates = env_raw$CTb_list,
  spotRates     = sp,
  nSpotCarriers = env_raw$nCO,
  nLanes        = env_raw$nL
)
addCapacityConstraintsCx(
  model_ptr             = model_ptr,
  transport_ptr         = transport_ptr,
  winnerKeys            = env_raw$winnerKey,
  winners               = env_raw$winner,
  bids                  = bids_list,
  carrierIdx            = env_raw$carrierIdx,
  strategicCapacities   = env_raw$Cb[t_lp + 1L, ],
  spotCapacities        = env_raw$Co[t_lp + 1L, ],
  nSpotSources          = env_raw$nCO,
  nLanes                = env_raw$nL
)
addStorageLimitConstraintsCx(
  model_ptr     = model_ptr,
  transport_ptr = transport_ptr,
  winnerKeys    = env_raw$winnerKey,
  winners       = env_raw$winner,
  bids          = bids_list,
  lanes         = lanes_list,
  nSpotCarriers = env_raw$nCO,
  nLanes        = env_raw$nL,
  nWarehouses   = c(env_raw$nI, env_raw$nJ),
  limits        = limits
)
addVolumeConstraintCx(
  model_ptr     = model_ptr,
  transport_ptr = transport_ptr,
  n             = n_vars,
  At            = env_raw$R
)

ptr_lp <- solveModelCx(model_ptr, transport_ptr)

stopifnot(
  "status mismatch"  = identical(ref_lp$status, ptr_lp$status),
  "objval mismatch"  = isTRUE(all.equal(ref_lp$objval, ptr_lp$objval, tolerance = 1e-10)),
  "x mismatch"       = isTRUE(all.equal(ref_lp$x,      ptr_lp$x,      tolerance = 1e-10))
)
cat(sprintf("  status: %s  objval (ref): %.4f  objval (ptr): %.4f  max |x diff|: %.2e\n",
            ptr_lp$status, ref_lp$objval, ptr_lp$objval,
            max(abs(ref_lp$x - ptr_lp$x))))

cat("[08] Test 3 PASSED\n")

## ============================================================
## 4. rolling_dp_ptr  vs  rolling_dp_cx  (full DP)
## ============================================================
cat("\n[08] Test 4: rolling_dp_ptr vs rolling_dp_cx\n")

dp_cx  <- rolling_dp_cx (env, JSON_PATH, numThreads = 4L)
dp_ptr <- rolling_dp_ptr(env, JSON_PATH, numThreads = 4L)

for (field in c("V", "Q", "pi_star", "pi_rand")) {
  stopifnot(isTRUE(all.equal(dp_cx[[field]], dp_ptr[[field]], tolerance = 1e-12)))
  cat(sprintf("  %-8s OK  max |diff| = %.2e\n", field,
              max(abs(dp_cx[[field]] - dp_ptr[[field]]), na.rm = TRUE)))
}

cat("[08] Test 4 PASSED\n")
cat("\n[08] All equivalence tests PASSED.\n")
