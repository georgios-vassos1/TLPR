## ============================================================
## 03_capacity_opt.R — Capacity strategy optimisation
##
## Reproduces Sections 5.2.2 and 5.2.3 of the paper:
##   • Exact L-BFGS-B optimisation over the fixed Table 1 scenario
##   • 1,000,000 Monte Carlo capacity vectors (Table 3)
##   • EV-based approximate optimisation + in/out-of-sample regret
##
## All gurobi:: calls replaced with solve_lp() (HiGHS).
##
## Outputs:
##   v0       — cost under initial capacity x0  (≈ 557.2)
##   x        — optim result under fixed scenario
##   opcosts  — 1M CostPerTEU values (Table 3)
##   x.ev     — EV-approximate optimal capacity
##   regrets  — in-sample regrets
##   oob.regrets — out-of-sample regrets
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(parallel)
})

## -- Config --------------------------------------------------
if (!exists("JSON_PATH")) source(file.path(dirname(sys.frame(1)$ofile), "config/instance1x1_4.R"))

## -- Instance (capacity experiment uses CTb/5) ---------------
env <- new.env()
jsonlite::fromJSON(JSON_PATH) |> list2env(envir = env)

env$alpha  <- ALPHA_SCALE * env$alpha
env$from_i <- apply(t(env$from_i), 2L, as.integer, simplify = FALSE)
env$to_j   <- apply(t(env$to_j),   2L, as.integer, simplify = FALSE)
env$CTb    <- env$CTb / 5.0   # framework agreement rate (≈ $2.94/TEU)

## -- Reservation costs (Table 2 column) ----------------------
env$nSources <- length(env$L_) + env$nL * env$nCO
v      <- rep(0.0, env$nSources * env$tau)
vdx    <- c(outer(seq(env$nSources - 1L),
                  seq(0L, env$nSources * (env$tau - 1L), by = env$nSources), "+"))
v[vdx] <- RESERVATION_COSTS

## -- Scenario space ------------------------------------------
env$scndx <- do.call(CartesianProductX, c(
  replicate(env$nI,  seq(env$nQ), simplify = FALSE),
  replicate(env$nJ,  seq(env$nD), simplify = FALSE),
  replicate(env$nCO, seq(env$nW), simplify = FALSE)))
env$scnpb <- apply(env$scndx, 1L, function(x)
  prod(env$Q$prob[x[seq(env$nI)]],
       env$D$prob[x[env$nI + seq(env$nJ)]],
       env$W$prob[x[env$nI + env$nJ + 1L]]))

## -- Fixed scenario: Table 1 ---------------------------------
Q_fixed      <- env$Q$vals[env$scndx[FIXED_VARPHIDX, 1L]]
D_fixed      <- env$D$vals[env$scndx[FIXED_VARPHIDX, 2L]]
env$CTo[]    <- env$W$vals[env$scndx[FIXED_VARPHIDX, 3L]]

## -- Build multi-period LP -----------------------------------
ccx  <- carrier_capacity_padded(env)
tlx  <- transition_logic(env, q = Q_fixed[seq(env$nI)], d = D_fixed[seq(env$nJ)])
slx  <- storage_limits(env,  q = Q_fixed[seq(env$nI)])
obj_ <- c(env$alpha, env$CTb, env$CTo[1L, ], env$alpha)

A   <- rbind(ccx$A, tlx$A, slx$A)
rhs <- c(ccx$rhs, tlx$rhs, slx$rhs)
sns <- c(ccx$sense, tlx$sense, slx$sense)

mp_model             <- multiperiod_expansion(env, Q_fixed, D_fixed, A, obj_, rhs, sns)
mp_model$modelsense  <- "min"
mp_model$vtype       <- rep("C", ncol(mp_model$A))

## Capacity constraint RHS indices (1x1: nCS+nCO=2 carriers, tau=4 periods)
capdx <- c(outer(seq(env$nCS + env$nCO), seq(0L, 6L * (env$tau - 1L), 6L), "+"))

## -- Objective: transport cost + reservation cost ------------
f <- function(model, x, v, capdx) {
  model$rhs[capdx] <- x
  opt <- solve_lp(model)
  if (is.null(opt$objval)) return(Inf)
  opt$objval + sum(v * x)
}

CostPerTEU <- function(model, x, v, capdx) {
  model$rhs[capdx] <- x
  opt <- solve_lp(model)
  if (is.null(opt$objval) || sum(opt$x) == 0) return(NA_real_)
  (opt$objval + sum(v * x)) / sum(opt$x)
}

## -- Section 5.2.2: exact capacity optimisation --------------
x0 <- c(t(cbind(env$Cb, env$Co)))
v0 <- f(mp_model, x0, v, capdx)
cat(sprintf("[03] Initial cost V1*(x0) = %.2f  (paper: %.1f)\n", v0, PAPER_COST_INITIAL))

ctrl <- list(maxit = 1000L, factr = 1e3, pgtol = 1e-9, parscale = rep(1.0, length(x0)))
x    <- optim(x0, f, model = mp_model, v = v, capdx = capdx,
              method = "L-BFGS-B", lower = 0.0, upper = 10.0, control = ctrl)

v_opt <- f(mp_model, round(x$par), v, capdx)
cat(sprintf("[03] Optimal cost  V1*(x*) = %.2f  (paper: %.1f)\n", v_opt, PAPER_COST_OPTIMAL))
cat(sprintf("[03] Reduction     = %.1f%%  (paper: %.1f%%)\n",
            (1 - v_opt / v0) * 100, PAPER_REDUCTION_PCT))
cat(sprintf("[03] Optimal capacities (rounded): %s\n", paste(round(x$par), collapse = " ")))

## -- Section 5.2.2: CostPerTEU comparison -------------------
cpu_initial <- CostPerTEU(mp_model, x0,         v, capdx)
cpu_optimal <- CostPerTEU(mp_model, round(x$par), v, capdx)
cat(sprintf("[03] CostPerTEU: initial=%.4f  optimal=%.4f\n", cpu_initial, cpu_optimal))

## -- Table 3: 1,000,000 random capacity vectors --------------
cat("[03] Sampling 1,000,000 random capacity vectors...\n")
M          <- 1000000L
capacities <- matrix(
  sample(seq(0L, 10L), (env$nCS + env$nCO) * env$tau * M, replace = TRUE),
  nrow = M)

cl <- makeCluster(detectCores() - 2L)
clusterExport(cl, c("mp_model", "v", "capdx", "CostPerTEU", "solve_lp"))
clusterEvalQ(cl, { library(TLPR); library(highs); library(methods) })

timex   <- Sys.time()
opcosts <- parApply(cl, capacities, 1L, CostPerTEU, model = mp_model, v = v, capdx = capdx)
elapsed <- as.numeric(Sys.time() - timex, units = "secs")
stopCluster(cl)

cat(sprintf("[03] Monte Carlo done in %.1f s\n", elapsed))
cat("[03] CostPerTEU distribution (Table 3):\n")
print(summary(opcosts))

## -- Section 5.2.3: EV-based approximate optimisation --------
EV <- function(cl, env, model, x, v, scnmat, qdx, ddx, wdx, capdx) {
  prob <- apply(scnmat, 1L, function(idx) prod(env$scnpb[idx]))
  N    <- nrow(scnmat)

  job <- function(k) {
    vk       <- scnmat[k, ]
    Q_k      <- env$Q$vals[env$scndx[vk, 1L]]
    D_k      <- env$D$vals[env$scndx[vk, 2L]]
    w_k      <- env$W$vals[env$scndx[vk, 3L]]
    model$rhs[qdx]    <- rep(Q_k, each = 2L)
    model$rhs[ddx]    <- D_k
    model$obj[wdx]    <- w_k
    model$rhs[capdx]  <- x
    opt <- solve_lp(model)
    if (is.null(opt$objval)) NA_real_ else opt$objval
  }

  optx  <- unlist(parLapply(cl, seq(N), job))
  optcx <- sum(v * x) + optx
  sum(optcx * (prob / sum(prob)), na.rm = TRUE)
}

## Scenario-update indices in the multi-period model
qdx   <- c(outer(c(3L, 5L), seq(0L, 6L * (env$tau - 1L), 6L), "+"))
ddx   <- seq(4L, 6L * env$tau, by = 6L)
wdx   <- seq(5L, 6L * env$tau, by = 5L)
capdx2 <- c(outer(1L:2L, seq(0L, 6L * (env$tau - 1L), 6L), "+"))

set.seed(123L)
N      <- 1000L
scnmat <- unique(matrix(
  sample(nrow(env$scndx), env$tau * N, replace = TRUE, prob = env$scnpb), nrow = N))

cl <- makeCluster(detectCores() - 2L)
clusterExport(cl, c("env", "mp_model", "v", "scnmat", "qdx", "ddx", "wdx", "capdx2", "solve_lp"))
clusterEvalQ(cl, { library(TLPR); library(highs); library(methods) })

ctrl2 <- list(maxit = 1000L, factr = 1e3, pgtol = 1e-9)
x.ev  <- optim(x0, EV, cl = cl, env = env, model = mp_model,
               v = v, scnmat = scnmat, qdx = qdx, ddx = ddx,
               wdx = wdx, capdx = capdx2,
               method = "L-BFGS-B", lower = 0.0, upper = Inf, control = ctrl2)
stopCluster(cl)
cat(sprintf("[03] EV-optimal capacities: %s\n", paste(round(x.ev$par), collapse = " ")))

## -- Regret analysis -----------------------------------------
C.opt <- function(env, model, x, v, varphidx, qdx, ddx, wdx, capdx) {
  Q_k <- env$Q$vals[env$scndx[varphidx, 1L]]
  D_k <- env$D$vals[env$scndx[varphidx, 2L]]
  w_k <- env$W$vals[env$scndx[varphidx, 3L]]
  model$rhs[qdx]   <- rep(Q_k, each = 2L)
  model$rhs[ddx]   <- D_k
  model$obj[wdx]   <- w_k
  model$rhs[capdx] <- x
  opt <- solve_lp(model)
  if (is.null(opt$objval)) Inf else sum(v * x) + opt$objval
}

compute_regret <- function(env, model, x.ev, v, varphidx, qdx, ddx, wdx, capdx, ctrl) {
  x_opt <- optim(x.ev$par, C.opt, env = env, model = model, v = v,
                 varphidx = varphidx, qdx = qdx, ddx = ddx, wdx = wdx, capdx = capdx,
                 method = "L-BFGS-B", lower = 0.0, upper = 10.0, control = ctrl)
  C.opt(env, model, x.ev$par,          v, varphidx, qdx, ddx, wdx, capdx) -
  C.opt(env, model, round(x_opt$par),  v, varphidx, qdx, ddx, wdx, capdx)
}

cl <- makeCluster(detectCores() - 2L)
clusterExport(cl, c("C.opt", "compute_regret", "env", "mp_model",
                    "x.ev", "v", "scnmat", "qdx", "ddx", "wdx", "capdx2", "ctrl2", "solve_lp"))
clusterEvalQ(cl, { library(TLPR); library(highs); library(methods) })

regrets <- parApply(cl, scnmat, 1L, compute_regret,
                    env = env, model = mp_model, x.ev = x.ev, v = v,
                    qdx = qdx, ddx = ddx, wdx = wdx, capdx = capdx2, ctrl = ctrl2)
stopCluster(cl)
cat("[03] In-sample regret summary:\n"); print(summary(regrets))

## Out-of-sample
set.seed(456L)
oob    <- unique(matrix(
  sample(nrow(env$scndx), env$tau * N, replace = TRUE, prob = env$scnpb), nrow = N))

cl <- makeCluster(detectCores() - 2L)
clusterExport(cl, c("C.opt", "compute_regret", "env", "mp_model",
                    "x.ev", "v", "oob", "qdx", "ddx", "wdx", "capdx2", "ctrl2", "solve_lp"))
clusterEvalQ(cl, { library(TLPR); library(highs); library(methods) })

oob.regrets <- parApply(cl, oob, 1L, compute_regret,
                        env = env, model = mp_model, x.ev = x.ev, v = v,
                        qdx = qdx, ddx = ddx, wdx = wdx, capdx = capdx2, ctrl = ctrl2)
stopCluster(cl)
cat("[03] Out-of-sample regret summary:\n"); print(summary(oob.regrets))
