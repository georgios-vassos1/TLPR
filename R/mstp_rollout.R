# ─────────────────────────────────────────────────────────────────────────────
# mstp_rollout.R — Evaluate an RTDP policy on MSTP instance noise paths
#
# Bridges TLPR (RTDP, VFA) with MSTP instance cost accounting, enabling
# apples-to-apples cost comparison against SDDP on shared noise trajectories.
#
# Requires the MSTP package for build_env() and the carrier-allocation LP
# structure.  The dependency is explicit: all MSTP calls use MSTP::.
# ─────────────────────────────────────────────────────────────────────────────

# Internal helper: solve one stage of the RTDP-guided policy.
#
# Identical constraint structure to MSTP:::.stage_lp().  The LP objective
# replaces the myopic holding-cost terms with VFA future-value gradients:
#
#   LP_coef[entry_out[k]] = VFA_slope[k] + entry_store_coef[k]
#
# When theta is all-zero this collapses exactly to MSTP:::.stage_lp().
#
# Action selection uses the VFA objective; the returned 'cost' is the ACTUAL
# stage cost (transport + real holding/shortage coefficients) so accumulated
# costs are comparable to MSTP::compute_vss() and MSTP::mstp_simulate()$obj.
#
# env      — MSTP env from MSTP::build_env(sim)
# theta    — VFA parameter from vfa_fit(V_cache[, t+1], tlpr_env)
# tlpr_env — TLPR dp_config env (needs tlpr_env$R)
# sim      — MSTP instance list
.stage_lp_rtdp <- function(t, entry_in, exitp_in, exitm_in, Q_t, D_t,
                             env, sim, theta, tlpr_env) {
  nI      <- env$nI
  nJ      <- env$nJ
  nCS     <- env$nCS
  nCO     <- env$nCO
  nL      <- env$nL
  nL_     <- env$nL_
  nMoves  <- nL_ + nCO * nL
  nVars   <- nMoves + nI + 2L * nJ
  tau_sim <- env$tau + 1L

  col_eo <- nMoves + seq_len(nI)
  col_xp <- nMoves + nI + seq_len(nJ)
  col_xm <- nMoves + nI + nJ + seq_len(nJ)

  p_vec  <- seq_len(nCO * nL)
  c_vec  <- ((p_vec - 1L) %/% nL) + 1L
  l_vec  <- ((p_vec - 1L) %%  nL) + 1L
  spot_t <- sim$spot_coef[(c_vec - 1L) * nL * tau_sim + (l_vec - 1L) * tau_sim + t]

  # ── VFA gradient coefficients ───────────────────────────────────────────────
  # The TLPR VFA includes a terminal salvage term (-alpha * s_final), making
  # raw slopes negative.  We add the myopic holding cost as a baseline so that
  # a zero VFA reproduces the myopic (.stage_lp) objective exactly, and any
  # positive look-ahead signal from RTDP is added on top.
  #
  # NOTE: a proper apples-to-apples comparison requires training RTDP under the
  # MSTP-compatible cost model (no salvage — task D1).  Until then, the returned
  # costs reflect the TLPR-trained policy evaluated under MSTP cost accounting.
  R    <- tlpr_env$R
  nR_o <- R + 1L
  nR_d <- 2L * R + 1L

  vfa_entry <- vapply(seq_len(nI), function(k) {
    r     <- min(as.integer(round(entry_in[k])) + 1L, nR_o)
    slope <- if (r < nR_o) theta$orig[r + 1L, k] - theta$orig[r, k]
             else          theta$orig[nR_o, k]    - theta$orig[nR_o - 1L, k]
    slope + sim$entry_store_coef[k]
  }, numeric(1L))

  vfa_exitp <- numeric(nJ)
  vfa_exitm <- numeric(nJ)
  for (j in seq_len(nJ)) {
    net   <- as.integer(round(exitp_in[j])) - as.integer(round(exitm_in[j]))
    net   <- max(-R, min(R, net))
    r     <- net + R + 1L
    slope <- if (r < nR_d) theta$dest[r + 1L, j] - theta$dest[r, j]
             else          theta$dest[nR_d, j]    - theta$dest[nR_d - 1L, j]
    vfa_exitp[j] <- slope + sim$exit_store_coef[j]
    vfa_exitm[j] <- -slope + sim$exit_short_coef[j]
  }

  obj_vfa <- c(sim$transport_coef, spot_t, vfa_entry, vfa_exitp, vfa_exitm)

  # Constraint matrix (identical to MSTP:::.stage_lp)
  nLc   <- env$nLc
  nCons <- nCS + nCO + nI + nJ + nJ
  is_ <- integer(0L); js_ <- integer(0L); xs_ <- numeric(0L)
  add <- function(i, j, x = 1.0) { is_ <<- c(is_, i); js_ <<- c(js_, j); xs_ <<- c(xs_, x) }

  for (k in seq_len(nCS))
    for (col in (nLc[k] + 1L):nLc[k + 1L]) add(k, col)
  for (cc in seq_len(nCO))
    for (col in nL_ + (cc - 1L) * nL + seq_len(nL)) add(nCS + cc, col)

  r0 <- nCS + nCO
  for (i in seq_len(nI)) {
    add(r0 + i, col_eo[i])
    for (col in env$from_i[[i]]) add(r0 + i, col)
  }
  r0 <- nCS + nCO + nI
  for (j in seq_len(nJ))
    for (col in env$to_j[[j]]) add(r0 + j, col)
  r0 <- nCS + nCO + nI + nJ
  for (j in seq_len(nJ)) {
    add(r0 + j, col_xp[j],  1.0)
    add(r0 + j, col_xm[j], -1.0)
    for (col in env$to_j[[j]]) add(r0 + j, col, -1.0)
  }

  A <- methods::as(
    Matrix::sparseMatrix(i = is_, j = js_, x = xs_, dims = c(nCons, nVars)),
    "CsparseMatrix"
  )
  sense <- c(rep("<=", nCS + nCO), rep("=", nI), rep("<=", nJ), rep("=", nJ))
  rhs   <- c(env$Cb[t, seq_len(nCS)],
             env$Co[t, seq_len(nCO)],
             entry_in + Q_t,
             sim$exit_capacity - exitp_in,
             exitp_in - exitm_in - D_t)

  res <- solve_lp(list(obj = obj_vfa, A = A, rhs = rhs,
                        sense = sense, modelsense = "min", vtype = NULL))
  x <- pmax(0, res$x)

  state_cost   <- sum(sim$entry_store_coef * entry_in) +
                  sum(sim$exit_store_coef  * exitp_in) +
                  sum(sim$exit_short_coef  * exitm_in)
  move_cost    <- sum(c(sim$transport_coef, spot_t) * x[seq_len(nMoves)])
  next_holding <- sum(sim$entry_store_coef * x[col_eo]) +
                  sum(sim$exit_store_coef  * x[col_xp]) +
                  sum(sim$exit_short_coef  * x[col_xm])

  list(cost      = state_cost + move_cost + next_holding,
       entry_out = x[col_eo],
       exitp_out = x[col_xp],
       exitm_out = x[col_xm])
}

#' Evaluate an RTDP policy on MSTP noise trajectories
#'
#' Simulates the RTDP greedy policy on the same noise paths produced by
#' \code{MSTP::mstp_simulate()}, enabling direct cost comparison against SDDP
#' on identical out-of-sample scenarios.
#'
#' At each stage the function solves the carrier-allocation LP with VFA
#' future-value gradient coefficients substituted for the myopic holding-cost
#' terms.  Because the VFA is separable and piecewise-linear, the gradient at
#' an integer grid point is a constant — the LP stays linear.  The reported
#' cost uses actual holding/shortage coefficients so the returned vector is on
#' the same scale as \code{MSTP::mstp_simulate()$obj}.
#'
#' When \code{theta} is all-zero the function reproduces the myopic
#' (\code{MSTP:::compute_vss}) rollout exactly.
#'
#' @param rtdp_res  Result from \code{rolling_dp_rtdp_p2()} or
#'   \code{rolling_dp_rtdp_p2_heur()}; must contain \code{V_cache}
#'   (\eqn{n_{Sdx} \times (\tau+1)}).
#' @param tlpr_env  TLPR env from \code{dp_config()} (needs \code{R}).
#' @param sim       MSTP instance list (from \code{MSTP::load_instances()} or
#'   \code{tlpr_env_to_mstp_inst()}).
#' @param noise     Noise matrix from \code{MSTP::mstp_simulate()$noise}:
#'   shape \eqn{(n\_trials \times \tau) \times (n_I + n_J)}.
#' @param n_trials  Number of trials (default: inferred from \code{nrow(noise)}).
#' @param n_cores   Parallel workers (default: all but 2).
#' @return Numeric vector of length \code{n_trials} — per-trial cumulative
#'   realised costs comparable to \code{MSTP::mstp_simulate()$obj}.
#' @export
rtdp_rollout_on_noise <- function(rtdp_res, tlpr_env, sim, noise,
                                   n_trials = NULL, n_cores = NULL) {
  env  <- MSTP::build_env(sim)
  tau  <- env$tau
  taux <- tau + 1L

  if (is.null(n_trials)) n_trials <- nrow(noise) %/% taux
  if (is.null(n_cores))  n_cores  <- max(1L, parallel::detectCores() - 2L)

  # Per-period VFA theta: theta_list[[t]] from V_cache[, t+1] =
  # value-to-go entering period t+1, used as future value at period t.
  theta_list <- vector("list", tau)
  for (t in seq_len(tau))
    theta_list[[t]] <- vfa_fit(rtdp_res$V_cache[, t + 1L], tlpr_env)

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  .stage_lp_rtdp_local <- .stage_lp_rtdp
  solve_lp_local        <- solve_lp

  parallel::clusterExport(cl,
    c("sim", "env", "tlpr_env", "theta_list", "noise", "tau", "taux"),
    envir = environment())
  parallel::clusterExport(cl,
    c(".stage_lp_rtdp_local", "solve_lp_local"),
    envir = environment())
  parallel::clusterEvalQ(cl, {
    library(Matrix)
    .stage_lp_rtdp <- .stage_lp_rtdp_local
    solve_lp       <- solve_lp_local
  })

  job <- function(trial) {
    entry_in <- as.numeric(sim$entry_stock_0)
    exitp_in <- as.numeric(sim$exit_stock_0)
    exitm_in <- as.numeric(sim$exit_short_0)
    total    <- 0.0
    for (t in seq_len(tau)) {
      row <- (trial - 1L) * taux + t
      Q_t <- noise[row, env$I_]
      D_t <- noise[row, env$nI + env$J_]
      s   <- .stage_lp_rtdp(t, entry_in, exitp_in, exitm_in, Q_t, D_t,
                              env, sim, theta_list[[t]], tlpr_env)
      total    <- total + s$cost
      entry_in <- s$entry_out
      exitp_in <- s$exitp_out
      exitm_in <- s$exitm_out
    }
    total + sum(sim$entry_store_coef * entry_in) +
            sum(sim$exit_store_coef  * exitp_in) +
            sum(sim$exit_short_coef  * exitm_in)
  }

  parallel::parLapply(cl, seq_len(n_trials), job) |> as.numeric()
}
