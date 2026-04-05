## ============================================================
## 06_dp_scalability.R — DP scalability benchmark
##
## Measures how the backward-induction pipeline scales with
## instance size. Three experiments:
##
##   A) State space (R): vary rate at fixed 1x1 topology
##   B) Time horizon (tau): vary tau at fixed 1x1, rate=1
##   C) Topology (nI, nJ): curse of dimensionality at rate=1
##
## Key scaling identities:
##   R     = 10 * rate
##   nSdx  = (R+1)^nI * (2R+1)^nJ   [state space]
##   nAdx  = R+1                      [action space]
##   nScen = nQ^nI * nD^nJ * nW^nCO  [scenario space]
##   Transit rows = tau * nSdx * nAdx * nScen
##
## Instances are generated programmatically via generate_instance().
## Configs with transit rows > MAX_TRANSIT_ROWS are skipped.
## ============================================================

suppressPackageStartupMessages(library(TLPR))

MAX_TRANSIT_ROWS <- 25e6L
N_THREADS        <- max(1L, parallel::detectCores() - 2L)

## -- Configurations ------------------------------------------
configs <- list(
  ## A: State space — vary rate at fixed 1x1
  list(nI=1L, nJ=1L, tau=4L, rate=0.5, group="A: State space (1x1)"),
  list(nI=1L, nJ=1L, tau=4L, rate=1.0, group="A: State space (1x1)"),
  list(nI=1L, nJ=1L, tau=4L, rate=2.0, group="A: State space (1x1)"),
  list(nI=1L, nJ=1L, tau=4L, rate=4.0, group="A: State space (1x1)"),

  ## B: Time horizon — vary tau at fixed 1x1, rate=1
  list(nI=1L, nJ=1L, tau= 4L, rate=1.0, group="B: Time horizon (1x1, R=10)"),
  list(nI=1L, nJ=1L, tau= 8L, rate=1.0, group="B: Time horizon (1x1, R=10)"),
  list(nI=1L, nJ=1L, tau=12L, rate=1.0, group="B: Time horizon (1x1, R=10)"),

  ## C: Topology — curse of dimensionality
  list(nI=1L, nJ=1L, tau=4L, rate=1.0, group="C: Topology (tau=4, rate=1)"),
  list(nI=2L, nJ=1L, tau=4L, rate=1.0, group="C: Topology (tau=4, rate=1)"),
  list(nI=1L, nJ=2L, tau=4L, rate=1.0, group="C: Topology (tau=4, rate=1)"),
  list(nI=2L, nJ=2L, tau=4L, rate=0.5, group="C: Topology (tau=4, rate=0.5 max)"),
  list(nI=2L, nJ=2L, tau=4L, rate=1.0, group="C: Topology (tau=4, rate=1)")
)

## -- Header --------------------------------------------------
fmt_hdr <- "%-30s  %4s  %4s  %4s  %8s  %5s  %6s  %14s  %8s  %8s  %8s\n"
fmt_row <- "%-30s  %4d  %4.1f  %4d  %8d  %5d  %6d  %14s  %8.2f  %8.2f  %8.1f\n"
fmt_skp <- "%-30s  %4d  %4.1f  %4d  %8d  %5d  %6d  %14s  %8s  %8s  %8s\n"
cat(sprintf(fmt_hdr, "Config", "nI*nJ", "rate", "tau",
            "nSdx", "nAdx", "nScen", "Transit rows",
            "Cx (s)", "DP (s)", "Mem (MB)"))
cat(strrep("-", 105), "\n")

last_group <- ""

for (cfg in configs) {
  nI   <- cfg$nI
  nJ   <- cfg$nJ
  tau  <- cfg$tau
  rate <- cfg$rate
  grp  <- cfg$group

  if (grp != last_group) { cat(sprintf("\n  [%s]\n", grp)); last_group <- grp }

  R     <- as.integer(round(10 * rate))
  nSI   <- R + 1L
  nSJ   <- 2L * R + 1L
  nSdx  <- nSI^nI * nSJ^nJ
  nAdx  <- R + 1L
  nScen <- 3L^nI * 3L^nJ * 2L    # nQ=3, nD=3, nW=2, nCO=1
  transit_rows <- tau * nSdx * nAdx * nScen

  label <- sprintf("%dx%d tau=%d R=%d", nI, nJ, tau, R)

  if (transit_rows > MAX_TRANSIT_ROWS) {
    cat(sprintf(fmt_skp, label, nI*nJ, rate, tau, nSdx, nAdx, nScen,
      sprintf("%.1fM — SKIP", transit_rows/1e6), "", "", ""))
    next
  }

  ## Generate instance and write to temp JSON
  tmpf <- tempfile(fileext = ".json")
  generate_instance(
    nI = nI, nJ = nJ, tau = tau,
    nB = 5L, nCS = 10L, nCO = 1L, rate = rate,
    seed = 42L, path = tmpf)

  ## Load and configure
  env <- new.env()
  jsonlite::fromJSON(tmpf) |> list2env(envir = env)
  as_int_list <- function(x)
    if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)
  env$from_i <- as_int_list(env$from_i)
  env$to_j   <- as_int_list(env$to_j)

  model <- create_model(env)
  dp_config(env)

  ## -- Transit (Cx) ------------------------------------------
  transit <- matrix(NA_real_, nrow = transit_rows, ncol = 6L)
  t0 <- proc.time()[["elapsed"]]
  for (t in seq(tau)) {
    resultCx <- computeEnvironmentCx(
      tmpf, t - 1L,
      seq(0L, env$R),
      env$Q$vals,
      numThreads = N_THREADS)
    transit[(t - 1L) * env$nSdx * env$nAdx * env$nScen + seq(nrow(resultCx)), ] <- resultCx
  }
  cx_sec <- proc.time()[["elapsed"]] - t0

  ## -- DP (Bellman recursion) --------------------------------
  t0    <- proc.time()[["elapsed"]]
  dp    <- dynamic_programming(env, transit)
  dp_sec <- proc.time()[["elapsed"]] - t0

  mem_mb <- as.numeric(
    object.size(transit) + object.size(dp$V) + object.size(dp$Q)) / 1024^2

  cat(sprintf(fmt_row, label, nI*nJ, rate, tau,
    nSdx, nAdx, nScen, format(transit_rows, big.mark=","),
    cx_sec, dp_sec, mem_mb))

  unlink(tmpf)
  rm(transit, dp, env, model)
  gc(verbose = FALSE)
}

cat("\n", strrep("-", 105), "\n", sep="")
cat(sprintf("Threads: %d  |  Transit row limit: %s\n",
  N_THREADS, format(MAX_TRANSIT_ROWS, big.mark=",")))
