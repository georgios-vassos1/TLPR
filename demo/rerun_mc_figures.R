## ============================================================
## rerun_mc_figures.R — Recompute Table 3, SAA/regret, and figures
##
## Run from the package root:
##   Rscript demo/rerun_mc_figures.R
##
## Set FIGURES_ONLY=1 to skip computation and just regenerate
## figures from cached data:
##   FIGURES_ONLY=1 Rscript demo/rerun_mc_figures.R
##
## Outputs:
##   report_tlpr/figs/density_plots.png
##   report_tlpr/figs/regret_density.png
##   report_tlpr/figs/invsoutregret.png
##   /tmp/tlpr_mc_full.rds   (full vectors for figure reruns)
##   /tmp/tlpr_mc_results.rds  (numerical summaries)
## ============================================================

suppressPackageStartupMessages({
  library(TLPR)
  library(parallel)
  library(ggplot2)
})

FIGURES_ONLY <- identical(Sys.getenv("FIGURES_ONLY"), "1")

## ---- Config -------------------------------------------------
DEMO_DIR      <- file.path(path.expand("~/drayage/TLPR"), "demo")
FIGS_DIR      <- file.path(path.expand("~/drayage/report_tlpr"), "figs")
FULL_RDS      <- "/tmp/tlpr_mc_full.rds"
RESULT_RDS    <- "/tmp/tlpr_mc_results.rds"

## Common theme: 11pt tick labels, 12pt axis titles and legend text
theme10 <- theme_bw() + theme(
  axis.text   = element_text(size = 11),
  axis.title  = element_text(size = 12),
  legend.text = element_text(size = 12),
  plot.title  = element_text(size = 12)
)

if (FIGURES_ONLY && file.exists(FULL_RDS)) {
  cat("[Mode] Loading cached data — skipping computation.\n")
  cache       <- readRDS(FULL_RDS)
  opcosts     <- cache$opcosts
  totcosts    <- cache$totcosts
  regrets     <- cache$regrets
  oob.regrets <- cache$oob.regrets
} else {

  source(file.path(DEMO_DIR, "config/instance1x1_4.R"))

  ## ---- Instance -----------------------------------------------
  env <- new.env()
  jsonlite::fromJSON(JSON_PATH) |> list2env(envir = env)

  env$alpha  <- ALPHA_SCALE * env$alpha
  env$from_i <- apply(t(env$from_i), 2L, as.integer, simplify = FALSE)
  env$to_j   <- apply(t(env$to_j),   2L, as.integer, simplify = FALSE)
  env$CTb    <- env$CTb / 5.0

  ## ---- Reservation costs --------------------------------------
  env$nSources <- length(env$L_) + env$nL * env$nCO
  v      <- rep(0.0, env$nSources * env$tau)
  vdx    <- c(outer(seq(env$nSources - 1L),
                    seq(0L, env$nSources * (env$tau - 1L), by = env$nSources), "+"))
  v[vdx] <- RESERVATION_COSTS

  ## ---- Scenario space -----------------------------------------
  env$scndx <- do.call(CartesianProductX, c(
    replicate(env$nI,  seq(env$nQ), simplify = FALSE),
    replicate(env$nJ,  seq(env$nD), simplify = FALSE),
    replicate(env$nCO, seq(env$nW), simplify = FALSE)))
  env$scnpb <- apply(env$scndx, 1L, function(x)
    prod(env$Q$prob[x[seq(env$nI)]],
         env$D$prob[x[env$nI + seq(env$nJ)]],
         env$W$prob[x[env$nI + env$nJ + 1L]]))

  ## ---- Fixed scenario -----------------------------------------
  Q_fixed   <- env$Q$vals[env$scndx[FIXED_VARPHIDX, 1L]]
  D_fixed   <- env$D$vals[env$scndx[FIXED_VARPHIDX, 2L]]
  env$CTo[] <- env$W$vals[env$scndx[FIXED_VARPHIDX, 3L]]

  ## ---- Build multi-period LP ----------------------------------
  ccx  <- carrier_capacity_padded(env)
  tlx  <- transition_logic(env, q = Q_fixed[seq(env$nI)], d = D_fixed[seq(env$nJ)])
  slx  <- storage_limits(env,  q = Q_fixed[seq(env$nI)])
  obj_ <- c(env$alpha, env$CTb, env$CTo[1L, ], env$alpha)

  A   <- rbind(ccx$A, tlx$A, slx$A)
  rhs <- c(ccx$rhs, tlx$rhs, slx$rhs)
  sns <- c(ccx$sense, tlx$sense, slx$sense)

  mp_model            <- multiperiod_expansion(env, Q_fixed, D_fixed, A, obj_, rhs, sns)
  mp_model$modelsense <- "min"
  mp_model$vtype      <- rep("C", ncol(mp_model$A))

  ## ---- Fix initial state S0 = (0, 8) -------------------------
  offset   <- env$nI + 2L * env$nJ
  s0_entry <- S0[seq(env$nI)]
  s0_j_pos <- pmax(S0[env$nI + seq(env$nJ)],  0L)
  s0_j_neg <- pmax(-S0[env$nI + seq(env$nJ)], 0L)
  s0_rhs   <- c(s0_entry, as.vector(rbind(s0_j_pos, s0_j_neg)))
  n_col    <- ncol(mp_model$A)
  A_s0     <- cbind(Matrix::Diagonal(offset),
                    Matrix::Matrix(0L, nrow = offset, ncol = n_col - offset))
  mp_model$A     <- rbind(mp_model$A, A_s0)
  mp_model$rhs   <- c(mp_model$rhs,   s0_rhs)
  mp_model$sense <- c(mp_model$sense, rep("=", offset))

  ## ---- Helpers ------------------------------------------------
  capdx <- c(outer(seq(env$nCS + env$nCO), seq(0L, 6L * (env$tau - 1L), 6L), "+"))

  CostPerTEU <- function(model, x, v, capdx) {
    model$rhs[capdx] <- x
    opt <- solve_lp(model)
    if (is.null(opt$objval) || sum(opt$x) == 0) return(NA_real_)
    (opt$objval + sum(v * x)) / sum(opt$x)
  }

  TotalCost <- function(model, x, v, capdx) {
    model$rhs[capdx] <- x
    opt <- solve_lp(model)
    if (is.null(opt$objval)) return(NA_real_)
    opt$objval + sum(v * x)
  }

  ## ==========================================================
  ## Part 1: Monte Carlo — 1,000,000 random capacity vectors
  ## ==========================================================
  cat("[MC] Sampling 1,000,000 random capacity vectors...\n")
  set.seed(0L)
  M          <- 1000000L
  capacities <- matrix(
    sample(seq(0L, 10L), (env$nCS + env$nCO) * env$tau * M, replace = TRUE),
    nrow = M)

  n_cores <- max(1L, detectCores() - 2L)
  cat(sprintf("[MC] Using %d cores\n", n_cores))
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("mp_model", "v", "capdx", "CostPerTEU", "TotalCost", "solve_lp"))
  clusterEvalQ(cl, { library(TLPR); library(highs); library(methods) })

  t0       <- Sys.time()
  opcosts  <- parApply(cl, capacities, 1L, CostPerTEU, model = mp_model, v = v, capdx = capdx)
  totcosts <- parApply(cl, capacities, 1L, TotalCost,  model = mp_model, v = v, capdx = capdx)
  elapsed  <- as.numeric(Sys.time() - t0, units = "secs")
  stopCluster(cl); on.exit(NULL)
  cat(sprintf("[MC] Done in %.1f s\n", elapsed))

  mc_cpu_summary <- summary(opcosts)
  mc_tot_summary <- summary(totcosts)
  cat("[MC] Cost-per-TEU summary:\n"); print(mc_cpu_summary)
  cat("[MC] Total Cost summary:\n");   print(mc_tot_summary)

  ## ==========================================================
  ## Part 2: EV-based SAA optimisation
  ## ==========================================================
  qdx    <- c(outer(c(3L, 5L), seq(0L, 6L * (env$tau - 1L), 6L), "+"))
  ddx    <- seq(4L, 6L * env$tau, by = 6L)
  wdx    <- seq(5L, 6L * env$tau, by = 5L)
  capdx2 <- c(outer(1L:2L, seq(0L, 6L * (env$tau - 1L), 6L), "+"))

  EV <- function(cl, env, model, x, v, scnmat, qdx, ddx, wdx, capdx) {
    prob <- apply(scnmat, 1L, function(idx) prod(env$scnpb[idx]))
    N    <- nrow(scnmat)
    job  <- function(k) {
      vk          <- scnmat[k, ]
      Q_k         <- env$Q$vals[env$scndx[vk, 1L]]
      D_k         <- env$D$vals[env$scndx[vk, 2L]]
      w_k         <- env$W$vals[env$scndx[vk, 3L]]
      model$rhs[qdx]   <- rep(Q_k, each = 2L)
      model$rhs[ddx]   <- D_k
      model$obj[wdx]   <- w_k
      model$rhs[capdx] <- x
      opt <- solve_lp(model)
      if (is.null(opt$objval)) NA_real_ else opt$objval
    }
    optx  <- unlist(parLapply(cl, seq(N), job))
    optcx <- sum(v * x) + optx
    sum(optcx * (prob / sum(prob)), na.rm = TRUE)
  }

  x0    <- c(t(cbind(env$Cb, env$Co)))
  ctrl2 <- list(maxit = 1000L, factr = 1e3, pgtol = 1e-9)

  set.seed(123L)
  N      <- 1000L
  scnmat <- unique(matrix(
    sample(nrow(env$scndx), env$tau * N, replace = TRUE, prob = env$scnpb), nrow = N))

  cat(sprintf("[SAA] Optimising over %d scenarios...\n", nrow(scnmat)))
  cl   <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("env", "mp_model", "v", "scnmat", "qdx", "ddx", "wdx", "capdx2", "solve_lp"))
  clusterEvalQ(cl, { library(TLPR); library(highs); library(methods) })
  t0   <- Sys.time()
  x.ev <- optim(x0, EV, cl = cl, env = env, model = mp_model,
                v = v, scnmat = scnmat, qdx = qdx, ddx = ddx,
                wdx = wdx, capdx = capdx2,
                method = "L-BFGS-B", lower = 0.0, upper = Inf, control = ctrl2)
  stopCluster(cl); on.exit(NULL)
  cat(sprintf("[SAA] Done in %.1f s; EV-optimal: %s\n",
              as.numeric(Sys.time() - t0, units = "secs"),
              paste(round(x.ev$par), collapse = " ")))

  ## ==========================================================
  ## Part 3: Regret analysis
  ## ==========================================================
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
    C.opt(env, model, x.ev$par,         v, varphidx, qdx, ddx, wdx, capdx) -
    C.opt(env, model, round(x_opt$par), v, varphidx, qdx, ddx, wdx, capdx)
  }

  cat("[Regret] In-sample...\n")
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("C.opt", "compute_regret", "env", "mp_model",
                      "x.ev", "v", "scnmat", "qdx", "ddx", "wdx", "capdx2", "ctrl2", "solve_lp"))
  clusterEvalQ(cl, { library(TLPR); library(highs); library(methods) })
  t0      <- Sys.time()
  regrets <- parApply(cl, scnmat, 1L, compute_regret,
                      env = env, model = mp_model, x.ev = x.ev, v = v,
                      qdx = qdx, ddx = ddx, wdx = wdx, capdx = capdx2, ctrl = ctrl2)
  stopCluster(cl); on.exit(NULL)
  cat(sprintf("[Regret] In-sample done in %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))
  cat("[Regret] In-sample summary:\n"); print(summary(regrets))

  set.seed(456L)
  oob <- unique(matrix(
    sample(nrow(env$scndx), env$tau * N, replace = TRUE, prob = env$scnpb), nrow = N))

  cat("[Regret] Out-of-sample...\n")
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("C.opt", "compute_regret", "env", "mp_model",
                      "x.ev", "v", "oob", "qdx", "ddx", "wdx", "capdx2", "ctrl2", "solve_lp"))
  clusterEvalQ(cl, { library(TLPR); library(highs); library(methods) })
  t0          <- Sys.time()
  oob.regrets <- parApply(cl, oob, 1L, compute_regret,
                          env = env, model = mp_model, x.ev = x.ev, v = v,
                          qdx = qdx, ddx = ddx, wdx = wdx, capdx = capdx2, ctrl = ctrl2)
  stopCluster(cl); on.exit(NULL)
  cat(sprintf("[Regret] OOB done in %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))
  cat("[Regret] Out-of-sample summary:\n"); print(summary(oob.regrets))

  ## Cache full vectors so figure reruns skip computation
  saveRDS(list(opcosts = opcosts, totcosts = totcosts,
               regrets = regrets, oob.regrets = oob.regrets),
          FULL_RDS)
  cat(sprintf("[Cache] Full vectors saved to %s\n", FULL_RDS))

  ## Numerical summaries
  mc_cpu_summary <- summary(opcosts)
  mc_tot_summary <- summary(totcosts)
  saveRDS(list(mc_cpu_summary = mc_cpu_summary, mc_tot_summary = mc_tot_summary,
               regret_summary = summary(regrets), oob_summary = summary(oob.regrets)),
          RESULT_RDS)
  fmt <- function(s) {
    vals <- as.numeric(s[c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")])
    sprintf("%.2f & %.2f & %.2f & %.2f & %.2f & %.2f",
            vals[1], vals[2], vals[3], vals[4], vals[5], vals[6])
  }
  cat("\n=== PASTE INTO PAPER ===\n")
  cat(sprintf("Cost-per-TEU & %s \\\\\n", fmt(mc_cpu_summary)))
  cat(sprintf("Total Cost   & %s \\\\\n", fmt(mc_tot_summary)))
}

## ============================================================
## Figures (runs in both modes)
## ============================================================

## Figure 5: density plots (Total Cost and Cost-per-TEU)
df_mc <- data.frame(
  CostPerTEU = opcosts[!is.na(opcosts)],
  TotalCost  = totcosts[!is.na(totcosts)]
)

p_tc <- ggplot(df_mc, aes(x = TotalCost)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(x = "Total Cost ($)", y = "Density") +
  theme10

p_cpu <- ggplot(df_mc, aes(x = CostPerTEU)) +
  geom_density(fill = "coral", alpha = 0.5) +
  labs(x = "Cost per TEU ($/TEU)", y = "Density") +
  theme10

## Figure widths match their \includegraphics display width in the LaTeX document
## (elsarticle preprint, textwidth ≈ 6.5 in) so fonts appear at their true pt size.
## density_plots: width=\textwidth=6.5 in
fig_density <- file.path(FIGS_DIR, "density_plots.png")
tryCatch({
  if (requireNamespace("cowplot", quietly = TRUE)) {
    p_combined <- cowplot::plot_grid(p_tc, p_cpu, ncol = 2L)
    cowplot::save_plot(fig_density, p_combined, ncol = 2, base_width = 3.25, base_height = 3.0)
  } else if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    ggsave(fig_density, p_tc + p_cpu, width = 6.5, height = 3.0, dpi = 300)
  } else {
    ggsave(fig_density, p_tc, width = 3.25, height = 3.0, dpi = 300)
    message("Install cowplot or patchwork for the combined figure")
  }
  cat(sprintf("[Figures] density_plots.png saved\n"))
}, error = function(e) message("density_plots error: ", e$message))

## Figure 6: regret density (in-sample vs out-of-sample)
df_reg <- data.frame(
  type   = c(rep("In-sample", length(regrets)), rep("Out-of-sample", length(oob.regrets))),
  regret = c(regrets, oob.regrets)
)

## regret_density: width=0.75\textwidth=4.875 in
fig_regdensity <- file.path(FIGS_DIR, "regret_density.png")
p_rd <- ggplot(df_reg, aes(x = regret, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("In-sample" = "steelblue", "Out-of-sample" = "coral")) +
  labs(x = "Regret ($)", y = "Density", fill = NULL) +
  theme10 + theme(legend.position = "bottom")
ggsave(fig_regdensity, p_rd, width = 4.875, height = 3.0, dpi = 300)
cat(sprintf("[Figures] regret_density.png saved\n"))

## Figure 7: in-sample vs out-of-sample regret — sorted Q-Q comparison.
## Both vectors come from independent scenario draws, so pairing by index
## is arbitrary.  Sorting both produces a distributional comparison: if
## in-sample and out-of-sample regret distributions match (good generalisation),
## the points fall near the identity line.
n_min   <- min(length(regrets), length(oob.regrets))
df_reg2 <- data.frame(ins = sort(regrets)[seq(n_min)],
                       oob = sort(oob.regrets)[seq(n_min)])

## invsoutregret: width=0.5\textwidth=3.25 in; equal x/y axis range
ax_max <- ceiling(max(df_reg2$ins, df_reg2$oob) / 10) * 10
fig_invsout <- file.path(FIGS_DIR, "invsoutregret.png")
p_io <- ggplot(df_reg2, aes(x = ins, y = oob)) +
  geom_point(alpha = 0.2, size = 1.0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  coord_equal(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
  labs(x = "In-sample regret ($)", y = "Out-of-sample regret ($)") +
  theme10
ggsave(fig_invsout, p_io, width = 3.25, height = 3.25, dpi = 300)
cat(sprintf("[Figures] invsoutregret.png saved\n"))

cat("\n[Done]\n")
