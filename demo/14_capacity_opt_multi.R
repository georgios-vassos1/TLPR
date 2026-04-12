## ============================================================
## 14_capacity_opt_multi.R — Capacity optimisation across topologies
##
## Runs dual-gradient L-BFGS-B capacity optimisation on four
## instances from the scalability table and reports cost
## reduction.  Target: new table in §5.2.2 of the paper.
##
## Configurations:
##   2×1 R=5   (|S| =    396)
##   1×2 R=5   (|S| =    726)
##   2×2 R=5   (|S| =  4,356)
##   2×2 R=10  (|S| = 53,361)  ← verify runtime first
##
## Usage:
##   Rscript demo/14_capacity_opt_multi.R
##   source("demo/14_capacity_opt_multi.R")   # from package root
##
## Outputs:
##   results  — list of per-config result lists
##   Prints a plain-text summary table and a LaTeX snippet.
## ============================================================

suppressPackageStartupMessages(library(TLPR))

## Handles both matrix-form (nCS small) and list-form (nCS large) from_i/to_j
as_int_list <- function(x)
  if (is.list(x)) lapply(x, as.integer) else apply(t(x), 2L, as.integer, simplify = FALSE)

ALPHA_SCALE <- 3.0
SEED        <- 42L

CONFIGS <- list(
  list(nI = 2L, nJ = 1L, rate = 0.5, label = "2x1, R=5",  nSdx_expected =    396L),
  list(nI = 1L, nJ = 2L, rate = 0.5, label = "1x2, R=5",  nSdx_expected =    726L),
  list(nI = 2L, nJ = 2L, rate = 0.5, label = "2x2, R=5",  nSdx_expected =  4356L),
  list(nI = 2L, nJ = 2L, rate = 1.0, label = "2x2, R=10", nSdx_expected = 53361L)
)

## ── Helper: build and solve capacity optimisation for one topology ──────────

run_cap_opt <- function(nI, nJ, rate, label, seed = SEED) {
  cat(sprintf("\n%s\n[%s]\n%s\n", strrep("=", 60), label, strrep("=", 60)))

  ## Generate instance (same parameters as 13_rtdp_benchmark.R)
  tmpf <- tempfile(fileext = ".json")
  generate_instance(nI = nI, nJ = nJ, tau = 4L,
                    nB = 5L, nCS = 10L, nCO = 1L,
                    rate = rate, regime = "balanced",
                    seed = seed, path = tmpf)

  env <- new.env()
  jsonlite::fromJSON(tmpf) |> list2env(envir = env)
  env$alpha  <- ALPHA_SCALE * env$alpha
  env$from_i <- as_int_list(env$from_i)
  env$to_j   <- as_int_list(env$to_j)
  env$CTb    <- env$CTb / 5.0          # framework agreement rate

  cat(sprintf("  nSdx=%d  nCS=%d  nCO=%d  tau=%d  R=%d\n",
              env$nSdx, env$nCS, env$nCO, env$tau, env$R))

  ## Scenario space
  env$scndx <- do.call(CartesianProductX, c(
    replicate(env$nI,  seq(env$nQ), simplify = FALSE),
    replicate(env$nJ,  seq(env$nD), simplify = FALSE),
    replicate(env$nCO, seq(env$nW), simplify = FALSE)))
  env$scnpb <- apply(env$scndx, 1L, function(idx)
    prod(env$Q$prob[idx[seq(env$nI)]],
         env$D$prob[idx[env$nI + seq(env$nJ)]],
         env$W$prob[idx[env$nI + env$nJ + 1L]]))

  ## Reference scenario: most probable outcome, same for all tau periods.
  ## Q_ref / D_ref must have length tau*nI and tau*nJ respectively so that
  ## multiperiod_expansion can index Q[(t-1)*nI + I_] and D[(t-1)*nJ + J_].
  mode_scn  <- which.max(env$scnpb)
  q_mode    <- env$Q$vals[env$scndx[mode_scn, seq(env$nI)]]
  d_mode    <- env$D$vals[env$scndx[mode_scn, env$nI + seq(env$nJ)]]
  Q_ref     <- rep(q_mode, env$tau)          # length tau*nI
  D_ref     <- rep(d_mode, env$tau)          # length tau*nJ
  env$CTo[] <- env$W$vals[env$scndx[mode_scn, env$nI + env$nJ + 1L]]

  ## Build single-period constraint matrix
  ccx  <- carrier_capacity_padded(env)
  tlx  <- transition_logic(env, q = Q_ref[seq(env$nI)], d = D_ref[seq(env$nJ)])
  slx  <- storage_limits(env,  q = Q_ref[seq(env$nI)])
  obj_ <- c(env$alpha, env$CTb, env$CTo[1L, ], env$alpha)

  A   <- rbind(ccx$A, tlx$A, slx$A)
  rhs <- c(ccx$rhs, tlx$rhs, slx$rhs)
  sns <- c(ccx$sense, tlx$sense, slx$sense)

  ## Expand to multi-period LP
  mp <- multiperiod_expansion(env, Q_ref, D_ref, A, obj_, rhs, sns)
  mp$modelsense <- "min"
  mp$vtype      <- rep("C", ncol(mp$A))

  ## Pin S0 = 0 (empty initial inventory for all nodes)
  offset <- env$nI + 2L * env$nJ
  n_col  <- ncol(mp$A)
  A_s0   <- cbind(Matrix::Diagonal(offset),
                  Matrix::Matrix(0L, nrow = offset, ncol = n_col - offset))
  mp$A     <- rbind(mp$A, A_s0)
  mp$rhs   <- c(mp$rhs,   numeric(offset))
  mp$sense <- c(mp$sense, rep("=", offset))

  ## Capacity constraint RHS indices (generalised: stride = nrow(A))
  n_per_period <- nrow(A)
  capdx <- c(outer(seq(env$nCS + env$nCO),
                   seq(0L, n_per_period * (env$tau - 1L), n_per_period), "+"))

  ## Initial capacity x0 from instance.
  ## Reservation costs: uniform per-unit charge equal to 2× the mean transport
  ## cost, creating a meaningful tradeoff between reservation spend and
  ## operational cost.  CTb may have more entries than (nCS+nCO) (one per lane),
  ## so we take its mean as the per-capacity-unit cost.
  x0 <- c(t(cbind(env$Cb, env$Co)))
  v  <- rep(mean(env$CTb) * 2.0, length(x0))

  ## Objective and dual-gradient functions
  dual_cache      <- new.env(parent = emptyenv())
  dual_cache$x    <- NULL
  dual_cache$dual <- NULL

  f <- function(model, x, v, capdx) {
    model$rhs[capdx] <- x
    opt <- solve_lp(model)
    if (is.null(opt$objval)) return(Inf)
    opt$objval + sum(v * x)
  }

  f_dual <- function(model, x, v, capdx) {
    model$rhs[capdx] <- x
    opt <- solve_lp(model)
    if (is.null(opt$objval)) return(Inf)
    dual_cache$x    <- x
    dual_cache$dual <- opt$dual
    opt$objval + sum(v * x)
  }

  gr_dual <- function(model, x, v, capdx) {
    if (!isTRUE(all.equal(dual_cache$x, x))) {
      model$rhs[capdx] <- x
      opt <- solve_lp(model)
      dual_cache$x    <- x
      dual_cache$dual <- opt$dual
    }
    dual_cache$dual[capdx] + v
  }

  ## Baseline cost at initial capacity
  t0 <- proc.time()["elapsed"]
  v0 <- f(mp, x0, v, capdx)
  cat(sprintf("  V(x0) = %.2f\n", v0))

  ## L-BFGS-B with analytical gradient
  ctrl <- list(maxit = 1000L, factr = 1e3, pgtol = 1e-9,
               parscale = rep(1.0, length(x0)))
  res  <- optim(x0, f_dual, gr = gr_dual,
                model = mp, v = v, capdx = capdx,
                method = "L-BFGS-B", lower = 0.0, upper = 10.0, control = ctrl)

  x_star <- round(res$par)
  v_opt  <- f(mp, x_star, v, capdx)
  redu   <- (1.0 - v_opt / v0) * 100.0
  elapsed <- proc.time()["elapsed"] - t0

  cat(sprintf("  V(x*) = %.2f  |  Reduction = %.1f%%  |  Time = %.1f s\n",
              v_opt, redu, elapsed))
  cat(sprintf("  Fn evals: %d  |  Gr evals: %d\n", res$counts[1L], res$counts[2L]))
  cat(sprintf("  x0    = %s\n", paste(x0,     collapse = " ")))
  cat(sprintf("  x*    = %s\n", paste(x_star, collapse = " ")))

  list(label    = label,
       nSdx     = env$nSdx,
       nCS      = env$nCS,
       nCO      = env$nCO,
       x0       = x0,
       x_star   = x_star,
       v0       = v0,
       v_opt    = v_opt,
       reduction = redu,
       elapsed  = elapsed,
       fn_evals = res$counts[1L])
}

## ── Run all configurations ─────────────────────────────────────────────────

results <- lapply(CONFIGS, function(cfg)
  run_cap_opt(cfg$nI, cfg$nJ, cfg$rate, cfg$label))

## ── Plain-text summary table ───────────────────────────────────────────────

cat("\n\n", strrep("=", 72), "\n", sep = "")
cat("CAPACITY OPTIMISATION ACROSS TOPOLOGIES\n")
cat(strrep("=", 72), "\n")
cat(sprintf("%-14s %8s %8s %8s %9s %7s\n",
            "Config", "|S|", "V(x0)", "V(x*)", "Reduc%", "Time(s)"))
cat(strrep("-", 72), "\n")
for (r in results)
  cat(sprintf("%-14s %8s %8.1f %8.1f %8.1f%% %7.1f\n",
              r$label,
              format(r$nSdx, big.mark = ","),
              r$v0, r$v_opt, r$reduction, r$elapsed))
cat(strrep("=", 72), "\n")

## ── LaTeX table snippet ────────────────────────────────────────────────────

cat("\n% ── LaTeX snippet for §5.2.2 ────────────────────────────\n")
cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\caption{Capacity optimisation results across topologies.\n")
cat("$V(x^0)$ and $V(x^*)$ are the multi-period LP costs at the\n")
cat("initial and optimised integer capacities respectively;\n")
cat("$S_0 = \\mathbf{0}$ throughout.}\n")
cat("\\label{tab:cap_multi}\n")
cat("\\begin{tabular}{lrrrr}\n")
cat("\\hline\n")
cat(paste0("Topology & $|\\mathcal{S}|$ & $V(x^0)$ & $V(x^*)$",
           " & Reduction \\\\\n"))
cat("\\hline\n")
for (r in results) {
  topo <- sub("x", "\\\\times ", sub(", R=(\\d+)", " ($R=\\1$)", r$label))
  cat(sprintf("%s & %s & %.1f & %.1f & %.1f\\%% \\\\\n",
              topo,
              format(r$nSdx, big.mark = ","),
              r$v0, r$v_opt, r$reduction))
}
cat("\\hline\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
