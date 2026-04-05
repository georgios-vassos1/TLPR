## ============================================================
## 04_simulation.R — Multi-policy simulation (tau=12 instance)
##
## Demonstrates simulate_system() by comparing four allocation
## policies over N Monte Carlo runs on the tau=12 instance.
## This script is standalone and does not depend on 01-03.
##
## Outputs:
##   result  — single-run simulation result (last run)
##   stats_  — (3 x npi) matrix: mean, lower CI, upper CI per policy
## ============================================================

suppressPackageStartupMessages(library(TLPR))

## -- Config --------------------------------------------------
TLPR_ROOT <- path.expand("~/drayage/TLPR")
JSON_PATH_12 <- file.path(TLPR_ROOT, "src/instances/instance1x1_12_001.json")

## -- Instance ------------------------------------------------
env <- new.env()
jsonlite::fromJSON(JSON_PATH_12) |> list2env(envir = env)
env$from_i <- apply(t(env$from_i), 2L, as.integer, simplify = FALSE)
env$to_j   <- apply(t(env$to_j),   2L, as.integer, simplify = FALSE)

cat(sprintf("[04] Instance: nI=%d  nJ=%d  tau=%d  R=%d  nCS=%d  nCO=%d\n",
            env$nI, env$nJ, env$tau, env$R, env$nCS, env$nCO))

## -- Single-period LP model ----------------------------------
model <- create_model(env)

## -- Allocation policies -------------------------------------
policy <- list(
  "myopic"             = optimal_assignment,
  "random"             = random_assignment,
  "capacitated_random" = capacitated_random_assignment,
  "heuristic"          = heuristic_assignment
)
npi <- length(policy)
cat(sprintf("[04] Policies: %s\n", paste(names(policy), collapse = ", ")))

## -- Simulation setup ----------------------------------------
## Sample a fixed spot-market cost matrix for all runs
set.seed(42L)
env$CTo <- matrix(
  TLPR:::sample_exogenous_state(
    sample,
    list(x = env$W$vals, size = env$tau * env$nCO,
         prob = env$W$prob, replace = TRUE)),
  ncol = env$nCO)

args <- list(
  model = model,
  obj_  = NULL,
  rhs_  = NULL,
  n     = NULL,
  k     = env$nvars,
  edx   = nrow(model$A)
)

exog <- list(
  Q = list(func   = sample,
           params = list(x = env$Q$vals, size = env$tau * env$nI,
                         prob = env$Q$prob, replace = TRUE)),
  D = list(func   = sample,
           params = list(x = env$D$vals, size = env$tau * env$nJ,
                         prob = env$D$prob, replace = TRUE))
)

## -- N Monte Carlo simulations --------------------------------
N     <- 50L
costs <- matrix(NA, nrow = env$tau, ncol = npi * N)
i     <- 1L
while (i <= N) {
  result <- simulate_system(env, random_volume, policy, args, exog = exog)
  if (isTRUE(result$status)) {
    costs[, (i - 1L) * npi + seq(npi)] <- result$cost
    i <- i + 1L
  }
}

## -- Compute summary statistics ------------------------------
stats_ <- do.call(cbind, lapply(seq(npi), function(p) {
  compute_stats(costs[, seq(p, npi * N, by = npi)], N = 1L)
}))

## -- Report --------------------------------------------------
cat("[04] Total cumulative cost (mean ± 95% CI):\n")
for (p in seq(npi)) {
  pcosts  <- costs[, seq(p, npi * N, by = npi)]
  totcost <- colSums(pcosts, na.rm = TRUE)
  cat(sprintf("  %-22s  mean=%7.2f  sd=%6.2f\n",
              names(policy)[p], mean(totcost), sd(totcost)))
}
