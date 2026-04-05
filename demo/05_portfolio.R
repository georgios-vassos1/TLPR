## ============================================================
## 05_portfolio.R — Portfolio contract evaluation (1x1, tau=4)
##
## Evaluates the per-period assignment cost over a grid of
## order quantities A_ using eval_portfolio() on the paper
## instance (1x1, tau=4).  No paper reference values are
## checked here; this script is a functional demonstration.
##
## Outputs:
##   cpx — list(states, cost, assignment) from eval_portfolio()
## ============================================================

suppressPackageStartupMessages(library(TLPR))

## -- Config (reuse paper instance if already loaded) ---------
if (!exists("JSON_PATH")) {
  TLPR_ROOT <- path.expand("~/drayage/TLPR")
  JSON_PATH <- file.path(TLPR_ROOT, "src/instances/instance1x1_4_001.json")
}

## -- Instance ------------------------------------------------
env <- new.env()
jsonlite::fromJSON(JSON_PATH) |> list2env(envir = env)
env$from_i <- apply(t(env$from_i), 2L, as.integer, simplify = FALSE)
env$to_j   <- apply(t(env$to_j),   2L, as.integer, simplify = FALSE)

cat(sprintf("[05] Instance: nI=%d  nJ=%d  tau=%d  R=%d\n",
            env$nI, env$nJ, env$tau, env$R))

## -- Single-period LP model ----------------------------------
model <- create_model(env)

## -- Portfolio action grid -----------------------------------
A_ <- seq(0L, env$R, by = 2L)
nA <- length(A_)
cat(sprintf("[05] Evaluating %d action levels: %s\n", nA, paste(A_, collapse = " ")))

## -- Exogenous scenario (fixed for reproducibility) ----------
set.seed(99L)
S0 <- matrix(c(0L, 8L), nrow = 1L)   # same initial state as DP experiments
Q  <- sample(env$Q$vals, env$tau * env$nI, replace = TRUE, prob = env$Q$prob)
D  <- sample(env$D$vals, env$tau * env$nJ, replace = TRUE, prob = env$D$prob)
env$CTo <- matrix(
  sample(env$W$vals, env$tau * env$nCO, replace = TRUE, prob = env$W$prob),
  ncol = env$nCO)

## -- Evaluate portfolio --------------------------------------
cpx <- eval_portfolio(env, model, A_, S0, Q, D)

## -- Report --------------------------------------------------
cat("[05] Mean per-period cost over action grid:\n")
for (t in seq(env$tau)) {
  avg_cost <- mean(cpx$cost[[t]], na.rm = TRUE)
  cat(sprintf("  t=%d  mean cost=%.4f\n", t, avg_cost))
}
cat(sprintf("[05] Total states explored: %d\n", nrow(cpx$states[[env$tau + 1L]])))
