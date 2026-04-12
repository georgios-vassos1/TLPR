# ─────────────────────────────────────────────────────────────────────────────
# rtdp.R — Iterative Cache-Based ADP (RTDP Phases 1 & 2)
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
# Common structure (both phases):
#   Persistent V_cache[nSdx × (tau+1)] — NA = unvisited.
#   Blended V_next = cached exact values (priority) + VFA fallback.
#   Backward sweep: bellmanUpdatePtr(stateSubset) → update cache → refit VFA.
# ─────────────────────────────────────────────────────────────────────────────

#' Iterative Cache-Based ADP (RTDP Phase 1)
#'
#' Runs \code{n_iter} backward-induction passes.  Each pass solves \code{K}
#' randomly selected states per period, caches the exact Bellman values, and
#' refits the separable VFA from all cached values.  The blended
#' \eqn{\hat{V}_\text{next}} (cache + VFA fallback) improves with each pass as
#' the cache grows and the VFA benefits from more exact targets.
#'
#' @param env Environment produced by \code{dp_config}.
#' @param jsonFile Path to the instance JSON file.
#' @param n_iter Number of backward-sweep iterations.  Default 10.
#' @param K States to solve per period per iteration.  \code{NULL} = full grid
#'   (equivalent to repeated exact DP, useful as a convergence upper bound).
#' @param numThreads OMP threads passed to \code{bellmanUpdatePtr}.  Default 8.
#' @param seed Integer seed for reproducibility of state sampling.
#' @return A named list:
#'   \describe{
#'     \item{V_cache}{\code{nSdx × (tau+1)} matrix of cached exact Bellman
#'       values.  \code{NA} entries are states never visited.}
#'     \item{V_approx}{\code{(tau+1) × nSdx} VFA approximation after the final
#'       iteration (same layout as \code{rolling_dp_vfa}).}
#'     \item{pi_star}{Policy matrix from the final iteration (same layout as
#'       \code{rolling_dp_ptr}).}
#'     \item{cache_coverage}{Integer vector of length \code{n_iter}: number of
#'       states cached at \code{t=1} after each iteration.}
#'     \item{theta_delta}{Numeric vector of length \code{n_iter}: max absolute
#'       VFA parameter change at \code{t=1} relative to the previous iteration.
#'       Measures convergence of the approximation.}
#'   }
#' @export
rolling_dp_rtdp <- function(env, jsonFile,
                            n_iter     = 10L,
                            K          = NULL,
                            numThreads = 8L,
                            seed       = NULL) {

  if (!is.null(seed)) set.seed(seed)

  prob     <- loadProblemDataCx(jsonFile)
  use_full <- is.null(K) || K >= env$nSdx
  if (!use_full && K < 1L) stop("K must be >= 1")

  # ── Sparse cache ─────────────────────────────────────────────────────────────
  # V_cache[s, t]: exact Bellman value for state s at period t; NA = unvisited.
  # Column tau+1 = terminal values (fully populated from the start).
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

  # ── Output accumulators ──────────────────────────────────────────────────────
  V_approx <- matrix(NA_real_, nrow = env$tau + 1L, ncol = env$nSdx)
  pi_star  <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  V_approx[env$tau + 1L, ] <- V_term

  cache_coverage  <- integer(n_iter)
  theta_delta     <- numeric(n_iter)
  # Per-iteration snapshot of V_approx at t=1 — enables convergence tracking.
  V_history_t1    <- matrix(NA_real_, nrow = n_iter, ncol = env$nSdx)
  pi_history_t1   <- matrix(NA_real_, nrow = n_iter, ncol = env$nSdx * env$nAdx)

  # ── Main iteration loop ──────────────────────────────────────────────────────
  for (iter in seq(n_iter)) {

    theta_prev <- theta

    # States to solve this iteration (same sample used for all periods)
    state_idx <- if (use_full) integer(0L) else sort(sample(env$nSdx, K))

    for (t in seq(env$tau, 1L)) {

      # V_next: blend cached exact values with VFA (cache takes priority)
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

      # Write exact values into cache (finite only — infeasible states left NA)
      v_t <- step$V_t
      if (use_full) {
        finite_ok <- is.finite(v_t)
        if (any(finite_ok)) V_cache[finite_ok, t] <- v_t[finite_ok]
      } else {
        finite_ok <- state_idx[is.finite(v_t[state_idx])]
        if (length(finite_ok) > 0L) V_cache[finite_ok, t] <- v_t[finite_ok]
      }

      # Refit VFA from ALL cached values at this period (growing set each iter)
      theta <- vfa_fit(V_cache[, t], env)

      V_approx[t, ] <- vfa_eval(theta, env)
      pi_star[t, ]  <- step$pi_star_t
    }

    # ── Per-iteration diagnostics ─────────────────────────────────────────────
    cache_coverage[iter] <- sum(!is.na(V_cache[, 1L]))

    theta_delta[iter] <- max(
      max(abs(theta$orig - theta_prev$orig), na.rm = TRUE),
      max(abs(theta$dest - theta_prev$dest), na.rm = TRUE)
    )

    # Snapshot t=1 approximation and policy after this iteration
    V_history_t1 [iter, ] <- V_approx[1L, ]
    pi_history_t1[iter, ] <- pi_star [1L, ]
  }

  list(
    V_cache        = V_cache,
    V_approx       = V_approx,
    pi_star        = pi_star,
    cache_coverage = cache_coverage,
    theta_delta    = theta_delta,
    V_history_t1   = V_history_t1,
    pi_history_t1  = pi_history_t1
  )
}

#' LP-Accurate RTDP with Trajectory-Guided State Sampling (Phase 2)
#'
#' Replaces the uniform random state selection of Phase 1 with LP-optimal
#' forward simulation.  Each iteration runs \code{n_traj} forward trajectories
#' via \code{simulateStepPtr}: at each step the LP is solved for all actions
#' and the greedy-optimal action under the current blended \eqn{\hat{V}} is
#' chosen.  Only the states visited on these trajectories receive a Bellman
#' update in the subsequent backward sweep — eliminating the explicit grid
#' dependency.
#'
#' Key scaling property: the number of LP solves per iteration is
#' \eqn{O(n_\text{traj} \times \tau \times n_A)} for forward simulation plus
#' \eqn{O(n_\text{traj} \times \tau \times n_A \times n_\text{scen})} for the
#' backward update — both independent of \eqn{|S|}.
#'
#' @param env Environment produced by \code{dp_config}.  Must have
#'   \code{nQdx}, \code{nDdx}, \code{nWdx}, \code{nScen} (set by
#'   \code{dp_config}).
#' @param jsonFile Path to the instance JSON file.
#' @param n_iter Number of iterations.  Default 10.
#' @param n_traj Forward trajectories per iteration.  Default 50.
#' @param s0 Initial state index (1-based) for all trajectories.  \code{NULL}
#'   (default) samples a fresh random starting state each trajectory, providing
#'   broad state-space coverage.  Set to an integer to fix a single start.
#' @param numThreads OMP threads for \code{bellmanUpdatePtr}.  Default 8.
#' @param seed Integer seed.  Default \code{NULL}.
#' @return Same structure as \code{rolling_dp_rtdp}: \code{V_cache},
#'   \code{V_approx}, \code{pi_star}, \code{cache_coverage},
#'   \code{theta_delta}, \code{V_history_t1}, \code{pi_history_t1}.
#'   Additionally: \code{n_visited} — integer matrix \code{n_iter × tau}
#'   recording unique states visited per period per iteration.
#' @export
rolling_dp_rtdp_p2 <- function(env, jsonFile,
                                n_iter     = 10L,
                                n_traj     = 50L,
                                s0         = NULL,
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

  # ── Output accumulators ──────────────────────────────────────────────────────
  V_approx <- matrix(NA_real_, nrow = env$tau + 1L, ncol = env$nSdx)
  pi_star  <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  V_approx[env$tau + 1L, ] <- V_term

  cache_coverage <- integer(n_iter)
  theta_delta    <- numeric(n_iter)
  V_history_t1   <- matrix(NA_real_, nrow = n_iter, ncol = env$nSdx)
  pi_history_t1  <- matrix(NA_real_, nrow = n_iter, ncol = env$nSdx * env$nAdx)
  n_visited      <- matrix(0L,        nrow = n_iter, ncol = env$tau)

  # Pre-compute blended V_next helper
  blend_v_next <- function(theta_cur, col) {
    V_next <- vfa_eval(theta_cur, env)
    cached <- V_cache[, col]
    known  <- !is.na(cached)
    if (any(known)) V_next[known] <- cached[known]
    V_next
  }

  # ── Main iteration loop ──────────────────────────────────────────────────────
  for (iter in seq(n_iter)) {

    theta_prev <- theta

    # ── Forward simulation: n_traj trajectories ────────────────────────────────
    # visited_t[[t]]: 0-based state indices visited at period t across all trajs
    visited_t <- vector("list", env$tau)

    for (traj in seq(n_traj)) {
      # Starting state (0-based)
      s <- if (is.null(s0)) sample(env$nSdx, 1L) - 1L else s0 - 1L

      for (t in seq(env$tau)) {
        # Blended V_next for period t+1
        V_next_tp1 <- blend_v_next(theta, t + 1L)

        # Sample scenario (1-based kdx for simulateStepPtr)
        kdx_r <- sample(env$nScen, 1L, prob = env$scnpb)

        # LP-optimal forward step: returns 0-based next_state_idx, action_idx
        step <- simulateStepPtr(
          problem_ptr  = prob,
          t            = t - 1L,          # 0-based period for C++
          state_idx    = s + 1L,          # 1-based for C++ wrapper
          scenario_kdx = kdx_r,           # 1-based for C++ wrapper
          stateSupport = as.double(env$stateSupport),
          flowSupport  = as.double(env$Q$vals),
          alpha        = env$alpha,
          V_next       = V_next_tp1
        )

        visited_t[[t]] <- c(visited_t[[t]], s)  # record 0-based state at period t
        s <- step$next_state_idx                 # 0-based next state
      }
    }

    # Deduplicate and convert to 1-based for bellmanUpdatePtr
    state_subsets <- lapply(visited_t, function(vs) sort(unique(vs)) + 1L)

    # ── Backward sweep over visited states ─────────────────────────────────────
    for (t in seq(env$tau, 1L)) {
      V_next <- blend_v_next(theta, t + 1L)

      state_idx <- state_subsets[[t]]
      n_visited[iter, t] <- length(state_idx)

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

      # Write finite values into cache
      v_t      <- step$V_t
      ok_states <- state_idx[is.finite(v_t[state_idx])]
      if (length(ok_states) > 0L) V_cache[ok_states, t] <- v_t[ok_states]

      theta         <- vfa_fit(V_cache[, t], env)
      V_approx[t, ] <- vfa_eval(theta, env)
      pi_star[t, ]  <- step$pi_star_t
    }

    # ── Per-iteration diagnostics ─────────────────────────────────────────────
    cache_coverage[iter] <- sum(!is.na(V_cache[, 1L]))

    theta_delta[iter] <- max(
      max(abs(theta$orig - theta_prev$orig), na.rm = TRUE),
      max(abs(theta$dest - theta_prev$dest), na.rm = TRUE)
    )

    V_history_t1 [iter, ] <- V_approx[1L, ]
    pi_history_t1[iter, ] <- pi_star [1L, ]
  }

  list(
    V_cache        = V_cache,
    V_approx       = V_approx,
    pi_star        = pi_star,
    cache_coverage = cache_coverage,
    theta_delta    = theta_delta,
    V_history_t1   = V_history_t1,
    pi_history_t1  = pi_history_t1,
    n_visited      = n_visited
  )
}
