# ─────────────────────────────────────────────────────────────────────────────
# rtdp.R — Iterative Cache-Based ADP (RTDP Phases 1, 2 & 3)
#
# Phase 1 — rolling_dp_rtdp:
#   Uniform random state sampling.  K states per period per iteration.
#   Validates that the cache/blend mechanism works; no grid removal.
#
# Phase 2 — rolling_dp_rtdp_p2:
#   LP-accurate forward simulation via simulateStepPtr().  n_traj trajectories
#   per iteration visit O(tau) states each; only those states are Bellman-
#   updated.  Removes the explicit grid dependency: visited-state count is
#   O(n_traj × tau), independent of |S|.
#
# Phase 3 (convergence + exploration):
#   Both functions support early stopping using a patience-windowed relative
#   change in V_approx[1,] (the value function at t=1).  The loop stops when
#   the maximum relative v_delta over the last `patience` iterations drops
#   below `tol`.  rolling_dp_rtdp_p2 additionally accepts epsilon for ε-greedy
#   trajectory exploration: with probability epsilon a random starting state is
#   drawn even when s0 is fixed, decaying linearly to 0 over n_iter.
#
# Convergence metric rationale:
#   theta_delta (absolute VFA parameter change) is scale-dependent and dominated
#   by sampling noise — it oscillates rather than monotonically decreasing.
#   v_delta = mean |V_approx[1,t] - V_approx[1,t-1]| / mean |V_approx[1,t-1]|
#   is scale-invariant and directly measures whether the value-function estimate
#   has stabilised.  The patience window requires all of the last `patience`
#   iterations to be below tol, preventing premature stopping on a lucky
#   low-noise iteration.
# ─────────────────────────────────────────────────────────────────────────────

# ── Internal helper: relative V_approx[1,] change ────────────────────────────
.v_delta <- function(v_cur, v_prev) {
  ok <- is.finite(v_cur) & is.finite(v_prev)
  if (!any(ok)) return(Inf)
  mean(abs(v_cur[ok] - v_prev[ok])) / (mean(abs(v_prev[ok])) + 1e-8)
}

#' Iterative Cache-Based ADP (RTDP Phase 1)
#'
#' Runs up to \code{n_iter} backward-induction passes with early stopping when
#' the value-function estimate has stabilised.  Each pass solves \code{K}
#' randomly selected states per period, caches the exact Bellman values, and
#' refits the separable VFA from all cached values.  The blended
#' \eqn{\hat{V}_\text{next}} (cache + VFA fallback) improves with each pass as
#' the cache grows and the VFA benefits from more exact targets.
#'
#' @param env Environment produced by \code{dp_config}.
#' @param jsonFile Path to the instance JSON file.
#' @param n_iter Maximum number of backward-sweep iterations.  Default 10.
#' @param K States to solve per period per iteration.  \code{NULL} = full grid.
#' @param tol Convergence tolerance.  The loop stops early once the maximum
#'   relative change in \code{V_approx[1,]} over the last \code{patience}
#'   iterations falls below \code{tol}.  Default \code{1e-2} (1\%).
#' @param min_iter Minimum iterations before early stopping is checked.
#'   Default \code{3L}.
#' @param patience Number of consecutive stable iterations required before
#'   stopping.  Default \code{3L}.
#' @param numThreads OMP threads passed to \code{bellmanUpdatePtr}.  Default 8.
#' @param seed Integer seed for reproducibility of state sampling.
#' @return A named list:
#'   \describe{
#'     \item{V_cache}{\code{nSdx × (tau+1)} matrix of cached exact Bellman
#'       values.  \code{NA} entries are states never visited.}
#'     \item{V_approx}{\code{(tau+1) × nSdx} VFA approximation after the final
#'       iteration.}
#'     \item{pi_star}{Policy matrix from the final iteration.}
#'     \item{cache_coverage}{Integer vector of length \code{n_done}.}
#'     \item{v_delta}{Numeric vector of length \code{n_done}: mean relative
#'       change in \code{V_approx[1,]} at each iteration.  The convergence
#'       criterion is \code{max(tail(v_delta, patience)) < tol}.}
#'     \item{n_done}{Actual iterations completed.}
#'     \item{V_history_t1}{Matrix \code{n_done × nSdx}: V_approx[1,] per iter.}
#'     \item{pi_history_t1}{Matrix \code{n_done × (nSdx*nAdx)}: pi_star[1,] per iter.}
#'   }
#' @export
rolling_dp_rtdp <- function(env, jsonFile,
                            n_iter     = 10L,
                            K          = NULL,
                            tol        = 1e-2,
                            min_iter   = 3L,
                            patience   = 3L,
                            numThreads = 8L,
                            seed       = NULL) {

  if (!is.null(seed)) set.seed(seed)

  prob     <- loadProblemDataCx(jsonFile)
  use_full <- is.null(K) || K >= env$nSdx
  if (!use_full && K < 1L) stop("K must be >= 1")

  # ── Sparse cache ─────────────────────────────────────────────────────────────
  V_cache <- matrix(NA_real_, nrow = env$nSdx, ncol = env$tau + 1L)

  # ── Terminal values ──────────────────────────────────────────────────────────
  V_term <- -c(
    cbind(
      apply(env$Sdx[, env$I_,          drop = FALSE], 2L,
            function(sdx) env$stateSupport[sdx]),
      apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
            function(sdx) pmax(env$extendedStateSupport[sdx], 0L)),
     -apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
            function(sdx) pmin(env$extendedStateSupport[sdx], 0L))
    ) %*% c(env$alpha))

  V_cache[, env$tau + 1L] <- V_term

  # ── Initial VFA ──────────────────────────────────────────────────────────────
  theta <- vfa_fit(V_term, env)

  # ── Output accumulators ───────────────────────────────────────────────────────
  V_approx <- matrix(NA_real_, nrow = env$tau + 1L, ncol = env$nSdx)
  pi_star  <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  V_approx[env$tau + 1L, ] <- V_term

  cache_coverage <- integer(0L)
  v_delta        <- numeric(0L)
  V_history_t1   <- list()
  pi_history_t1  <- list()

  iter      <- 0L
  V_prev_t1 <- NULL

  while (iter < n_iter) {
    iter <- iter + 1L

    state_idx <- if (use_full) integer(0L) else sort(sample(env$nSdx, K))

    for (t in seq(env$tau, 1L)) {

      V_next      <- vfa_eval(theta, env)
      cached_next <- V_cache[, t + 1L]
      known       <- !is.na(cached_next)
      if (any(known)) V_next[known] <- cached_next[known]

      step <- bellmanUpdatePtr(
        problem_ptr    = prob,
        t              = t - 1L,
        stateSupport   = as.double(env$stateSupport),
        flowSupport    = as.double(env$Q$vals),
        scnpb          = env$scnpb,
        alpha          = env$alpha,
        V_next         = V_next,
        numThreads     = numThreads,
        stateSubset    = state_idx
      )

      v_t <- step$V_t
      if (use_full) {
        finite_ok <- is.finite(v_t)
        if (any(finite_ok)) V_cache[finite_ok, t] <- v_t[finite_ok]
      } else {
        finite_ok <- state_idx[is.finite(v_t[state_idx])]
        if (length(finite_ok) > 0L) V_cache[finite_ok, t] <- v_t[finite_ok]
      }

      theta <- vfa_fit(V_cache[, t], env)

      V_approx[t, ] <- vfa_eval(theta, env)
      pi_star[t, ]  <- step$pi_star_t
    }

    # ── Per-iteration diagnostics ─────────────────────────────────────────────
    cache_coverage <- c(cache_coverage, sum(!is.na(V_cache[, 1L])))

    vd <- if (is.null(V_prev_t1)) Inf else .v_delta(V_approx[1L, ], V_prev_t1)
    v_delta   <- c(v_delta, vd)
    V_prev_t1 <- V_approx[1L, ]

    V_history_t1 [[iter]] <- V_approx[1L, ]
    pi_history_t1[[iter]] <- pi_star [1L, ]

    # ── Convergence check ─────────────────────────────────────────────────────
    if (iter >= min_iter && length(v_delta) >= patience) {
      if (max(tail(v_delta, patience)) < tol) break
    }
  }

  list(
    V_cache        = V_cache,
    V_approx       = V_approx,
    pi_star        = pi_star,
    cache_coverage = cache_coverage,
    v_delta        = v_delta,
    n_done         = iter,
    V_history_t1   = do.call(rbind, V_history_t1),
    pi_history_t1  = do.call(rbind, pi_history_t1)
  )
}

#' LP-Accurate RTDP with Trajectory-Guided State Sampling (Phase 2 / Phase 3)
#'
#' Replaces the uniform random state selection of Phase 1 with LP-optimal
#' forward simulation.  Each iteration runs \code{n_traj} forward trajectories
#' via \code{simulateStepPtr}: at each step the LP is solved for all actions
#' and the greedy-optimal action under the current blended \eqn{\hat{V}} is
#' chosen.  Only the states visited on these trajectories receive a Bellman
#' update in the subsequent backward sweep — eliminating the explicit grid
#' dependency.
#'
#' Phase 3 extensions:
#' \itemize{
#'   \item \strong{Convergence stopping}: the loop terminates early when the
#'     maximum relative change in \code{V_approx[1,]} over the last
#'     \code{patience} iterations falls below \code{tol}.
#'   \item \strong{ε-greedy exploration}: with probability \code{epsilon *
#'     (1 - (iter-1)/n_iter)} a fresh random starting state is drawn even when
#'     \code{s0} is fixed, decaying linearly to 0.  When \code{s0 = NULL}
#'     exploration is always active.
#' }
#'
#' @param env Environment produced by \code{dp_config}.
#' @param jsonFile Path to the instance JSON file.
#' @param n_iter Maximum number of iterations.  Default 10.
#' @param n_traj Forward trajectories per iteration.  Default 50.
#' @param s0 Initial state index (1-based).  \code{NULL} = random start.
#' @param epsilon ε-greedy exploration rate (Phase 3).  Default \code{0.2}.
#' @param tol Convergence tolerance on relative \code{V_approx[1,]} change.
#'   Default \code{1e-2} (1\%).
#' @param min_iter Minimum iterations before early stopping.  Default \code{3L}.
#' @param patience Consecutive stable iterations required before stopping.
#'   Default \code{3L}.
#' @param numThreads OMP threads for \code{bellmanUpdatePtr}.  Default 8.
#' @param seed Integer seed.  Default \code{NULL}.
#' @return Named list with \code{V_cache}, \code{V_approx}, \code{pi_star},
#'   \code{cache_coverage}, \code{v_delta}, \code{n_done},
#'   \code{V_history_t1}, \code{pi_history_t1}, \code{n_visited}.
#' @export
rolling_dp_rtdp_p2 <- function(env, jsonFile,
                                n_iter     = 10L,
                                n_traj     = 50L,
                                s0         = NULL,
                                epsilon    = 0.2,
                                tol        = 1e-2,
                                min_iter   = 3L,
                                patience   = 3L,
                                numThreads = 8L,
                                seed       = NULL) {

  if (!is.null(seed)) set.seed(seed)

  prob <- loadProblemDataCx(jsonFile)

  # ── Terminal values ──────────────────────────────────────────────────────────
  V_term <- -c(
    cbind(
      apply(env$Sdx[, env$I_,          drop = FALSE], 2L,
            function(sdx) env$stateSupport[sdx]),
      apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
            function(sdx) pmax(env$extendedStateSupport[sdx], 0L)),
     -apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
            function(sdx) pmin(env$extendedStateSupport[sdx], 0L))
    ) %*% c(env$alpha))

  # ── Sparse cache ─────────────────────────────────────────────────────────────
  V_cache <- matrix(NA_real_, nrow = env$nSdx, ncol = env$tau + 1L)
  V_cache[, env$tau + 1L] <- V_term

  # ── Initial VFA ──────────────────────────────────────────────────────────────
  theta <- vfa_fit(V_term, env)

  # ── Output accumulators ───────────────────────────────────────────────────────
  V_approx <- matrix(NA_real_, nrow = env$tau + 1L, ncol = env$nSdx)
  pi_star  <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  V_approx[env$tau + 1L, ] <- V_term

  cache_coverage <- integer(0L)
  v_delta        <- numeric(0L)
  V_history_t1   <- list()
  pi_history_t1  <- list()
  n_visited_list <- list()

  iter      <- 0L
  V_prev_t1 <- NULL

  # Pre-compute blended V_next helper
  blend_v_next <- function(theta_cur, col) {
    V_next <- vfa_eval(theta_cur, env)
    cached <- V_cache[, col]
    known  <- !is.na(cached)
    if (any(known)) V_next[known] <- cached[known]
    V_next
  }

  while (iter < n_iter) {
    iter <- iter + 1L

    # ε-greedy decay: linearly anneals from epsilon to 0 (no-op when s0=NULL)
    epsilon_t <- if (is.null(s0)) 1.0 else epsilon * (1.0 - (iter - 1L) / n_iter)

    # ── Forward simulation ─────────────────────────────────────────────────────
    visited_t <- vector("list", env$tau)

    for (traj in seq(n_traj)) {
      s <- if (is.null(s0) || runif(1L) < epsilon_t) {
        sample(env$nSdx, 1L) - 1L
      } else {
        s0 - 1L
      }

      for (t in seq(env$tau)) {
        V_next_tp1 <- blend_v_next(theta, t + 1L)
        kdx_r      <- sample(env$nScen, 1L, prob = env$scnpb)

        step <- simulateStepPtr(
          problem_ptr  = prob,
          t            = t - 1L,
          state_idx    = s + 1L,
          scenario_kdx = kdx_r,
          stateSupport = as.double(env$stateSupport),
          flowSupport  = as.double(env$Q$vals),
          alpha        = env$alpha,
          V_next       = V_next_tp1
        )

        visited_t[[t]] <- c(visited_t[[t]], s)
        s <- step$next_state_idx
      }
    }

    state_subsets <- lapply(visited_t, function(vs) sort(unique(vs)) + 1L)

    # ── Backward sweep ─────────────────────────────────────────────────────────
    nv_iter <- integer(env$tau)
    for (t in seq(env$tau, 1L)) {
      V_next    <- blend_v_next(theta, t + 1L)
      state_idx <- state_subsets[[t]]
      nv_iter[t] <- length(state_idx)

      step <- bellmanUpdatePtr(
        problem_ptr  = prob,
        t            = t - 1L,
        stateSupport = as.double(env$stateSupport),
        flowSupport  = as.double(env$Q$vals),
        scnpb        = env$scnpb,
        alpha        = env$alpha,
        V_next       = V_next,
        numThreads   = numThreads,
        stateSubset  = state_idx
      )

      v_t       <- step$V_t
      ok_states <- state_idx[is.finite(v_t[state_idx])]
      if (length(ok_states) > 0L) V_cache[ok_states, t] <- v_t[ok_states]

      theta         <- vfa_fit(V_cache[, t], env)
      V_approx[t, ] <- vfa_eval(theta, env)
      pi_star[t, ]  <- step$pi_star_t
    }

    # ── Per-iteration diagnostics ─────────────────────────────────────────────
    cache_coverage     <- c(cache_coverage, sum(!is.na(V_cache[, 1L])))
    n_visited_list[[iter]] <- nv_iter

    vd <- if (is.null(V_prev_t1)) Inf else .v_delta(V_approx[1L, ], V_prev_t1)
    v_delta   <- c(v_delta, vd)
    V_prev_t1 <- V_approx[1L, ]

    V_history_t1 [[iter]] <- V_approx[1L, ]
    pi_history_t1[[iter]] <- pi_star [1L, ]

    # ── Convergence check ─────────────────────────────────────────────────────
    if (iter >= min_iter && length(v_delta) >= patience) {
      if (max(tail(v_delta, patience)) < tol) break
    }
  }

  list(
    V_cache        = V_cache,
    V_approx       = V_approx,
    pi_star        = pi_star,
    cache_coverage = cache_coverage,
    v_delta        = v_delta,
    n_done         = iter,
    V_history_t1   = do.call(rbind, V_history_t1),
    pi_history_t1  = do.call(rbind, pi_history_t1),
    n_visited      = do.call(rbind, n_visited_list)
  )
}
