## ============================================================
## 01_transit.R — Compute and cache the transit matrix
##
## Tries to load from TRANSIT_CACHE; if absent, computes via
## computeEnvironmentCx (C++/OpenMP, ~0.5 s) with an Rx
## cross-check on a random subset, then saves to disk.
##
## Outputs (in parent environment or .GlobalEnv):
##   env     — configured environment
##   model   — single-period HiGHS model
##   transit — tau * nSdx * nAdx * nScen x 6 matrix
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(foreach)
  library(doParallel)
})

## -- Config --------------------------------------------------
.demo_dir <- local({
  ofiles <- Filter(Negate(is.null), lapply(sys.frames(), `[[`, "ofile"))
  if (length(ofiles)) dirname(normalizePath(tail(ofiles, 1L)[[1L]])) else {
    farg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(farg)) dirname(normalizePath(sub("--file=", "", farg[1L]))) else getwd()
  }
})
if (!exists("JSON_PATH")) source(file.path(.demo_dir, "config/instance1x1_4.R"))

## -- Instance ------------------------------------------------
env <- new.env()
jsonlite::fromJSON(JSON_PATH) |> list2env(envir = env)

env$alpha  <- ALPHA_SCALE * env$alpha
env$from_i <- apply(t(env$from_i), 2L, as.integer, simplify = FALSE)
env$to_j   <- apply(t(env$to_j),   2L, as.integer, simplify = FALSE)

model <- create_model(env)
dp_config(env)

cat(sprintf("[01] |S|=%d  |A|=%d  |Omega|=%d  tau=%d\n",
            env$nSdx, env$nAdx, env$nScen, env$tau))
cat(sprintf("[01] Transit rows: %d\n", env$tau * env$nSdx * env$nAdx * env$nScen))

## -- Load or compute -----------------------------------------
if (file.exists(TRANSIT_CACHE)) {
  cat("[01] Loading transit from cache...\n")
  transit <- readRDS(TRANSIT_CACHE)
  cat(sprintf("[01] Loaded: %d rows, %d NA\n", nrow(transit), sum(is.na(transit[, 1L]))))
} else {
  cat("[01] Computing transit via computeEnvironmentCx (C++/OpenMP)...\n")
  timex   <- Sys.time()
  transit <- matrix(NA, nrow = env$tau * env$nSdx * env$nAdx * env$nScen, ncol = 6L)

  for (t in seq(env$tau)) {
    cat(sprintf("  t = %d\n", t))
    resultCx <- computeEnvironmentCx(
      JSON_PATH, t - 1L,
      seq(0L, env$R),
      env$Q$vals,
      numThreads = max(1L, parallel::detectCores() - 2L))
    transit[(t - 1L) * env$nSdx * env$nAdx * env$nScen + seq(nrow(resultCx)), ] <- resultCx
  }

  elapsed <- as.numeric(Sys.time() - timex, units = "secs")
  cat(sprintf("[01] Cx done in %.1f s | NA rows: %d / %d\n",
              elapsed, sum(is.na(transit[, 1L])), nrow(transit)))

  ## -- Cross-check: compare Cx vs Rx on t=1, first 20 states --
  cat("[01] Cross-checking Cx vs Rx on t=1, states 1:20...\n")
  n_chunks <- min(4L, parallel::detectCores() - 1L)
  chunks   <- chunkup(20L, n_chunks)
  cl       <- makeCluster(n_chunks)
  clusterExport(cl, c("env", "model"))
  registerDoParallel(cl)
  resultRx <- foreach(idx = seq(n_chunks), .packages = "TLPR", .combine = "rbind") %dopar% {
    computeEnvironmentRx(env, model, 1L, chunks[idx], chunks[idx + 1L] - 1L)
  }
  stopCluster(cl)

  n_cx      <- nrow(resultRx)
  cx_sub    <- transit[seq(n_cx), ]
  feas      <- !is.na(cx_sub[, 1L]) & !is.na(resultRx[, 1L])
  ns_agree  <- all(cx_sub[feas, 1L] == resultRx[feas, 1L])
  cost_agree <- isTRUE(all.equal(cx_sub[feas, 2L], resultRx[feas, 2L], tolerance = 1e-6))
  cat(sprintf("[01] Cross-check — next_state match: %s | cost match: %s\n",
              ns_agree, cost_agree))
  if (!ns_agree || !cost_agree) warning("[01] Cx/Rx discrepancy detected — inspect before proceeding.")

  saveRDS(transit, TRANSIT_CACHE)
  cat(sprintf("[01] Transit saved to: %s\n", TRANSIT_CACHE))
}
