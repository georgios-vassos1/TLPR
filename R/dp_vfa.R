# ─────────────────────────────────────────────────────────────────────────────
# dp_vfa.R — Approximate Value Iteration with separable VFA
#
# Drop-in companion to rolling_dp_ptr / rolling_dp_cx.  Uses the same
# bellmanUpdatePtr call but replaces the exact continuation vector V[t+1,]
# with the separable approximation V̂(s) = Φθ evaluated on the full grid.
#
# Backward sweep (t = tau … 1):
#   1. Evaluate V̂ on full grid from current theta         O(nSdx · (nI+nJ))
#   2. Call bellmanUpdatePtr with V_next = V̂             (same as exact DP)
#   3. Fit new theta from the returned V_t targets         O(nSdx · p)
#   4. Optionally project each vₖ to be concave
#
# The only difference from exact DP is step 1: the continuation value is the
# VFA approximation rather than the previously computed exact period vector.
# No C++ changes are required for this prototype.
# ─────────────────────────────────────────────────────────────────────────────

#' Approximate Value Iteration with Separable VFA
#'
#' Runs backward induction like \code{rolling_dp_ptr} but approximates the
#' value function at each period with a separable piecewise-constant function
#' \eqn{\hat{V}(s) = \sum_k v_k(s_k)} fit by OLS.  The full state grid is
#' solved each period (no state sampling); the approximation lies in the
#' \emph{continuation} value passed to the Bellman operator.
#'
#' @param env Environment produced by \code{dp_config}.
#' @param jsonFile Path to the instance JSON file.
#' @param numThreads Number of OMP threads passed to \code{bellmanUpdatePtr}.
#'   Default 8.
#' @param project_concavity Logical; if \code{TRUE} (default) each 1-D value
#'   function is projected onto concave piecewise-linear functions after
#'   fitting.  Disable to isolate the effect of the concavity constraint.
#' @param traversalOrder Passed through to \code{bellmanUpdatePtr}.
#' @param chunkSize Passed through to \code{bellmanUpdatePtr}.
#' @return A named list:
#'   \describe{
#'     \item{V_approx}{Matrix (\code{tau+1}) × \code{nSdx} of approximated
#'       value functions.  Row 1 = period 1 approximation; row \code{tau+1} =
#'       terminal values (exact).}
#'     \item{theta_list}{List of length \code{tau+1}; element \code{t} holds
#'       the \code{theta} fitted at period \code{t} (terminal = element
#'       \code{tau+1}).}
#'     \item{Q}{Matrix (\code{tau}) × (\code{nSdx * nAdx}) of Q-values
#'       computed under the VFA continuation — same structure as
#'       \code{rolling_dp_ptr} output.}
#'     \item{pi_star}{Greedy policy matrix — same structure as
#'       \code{rolling_dp_ptr}.}
#'     \item{pi_rand}{Softmax policy matrix.}
#'   }
#' @export
rolling_dp_vfa <- function(env, jsonFile,
                           K                 = NULL,
                           numThreads        = 8L,
                           project_concavity = FALSE,
                           traversalOrder    = "lexicographic",
                           chunkSize         = 32L,
                           seed              = NULL) {

  if (!is.null(seed)) set.seed(seed)

  prob <- loadProblemDataCx(jsonFile)

  # Full grid if K is NULL or K >= nSdx
  use_full <- is.null(K) || K >= env$nSdx
  if (!use_full && K < 1L)
    stop("K must be >= 1")

  V_approx   <- matrix(NA_real_, nrow = env$tau + 1L, ncol = env$nSdx)
  Q          <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  pi_star    <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  pi_rand    <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  theta_list <- vector("list", env$tau + 1L)

  # ── Terminal values (identical to rolling_dp_ptr) ──────────────────────────
  V_term <- -c(
    cbind(
      apply(env$Sdx[, env$I_,          drop = FALSE], 2L,
            function(sdx) env$stateSupport[sdx]),
      apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
            function(sdx) pmax(env$extendedStateSupport[sdx], 0L)),
     -apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
            function(sdx) pmin(env$extendedStateSupport[sdx], 0L))
    ) %*% c(env$alpha))

  V_approx[env$tau + 1L, ] <- V_term

  # Fit initial theta from terminal values
  theta <- vfa_fit(V_term, env)
  if (project_concavity) theta <- vfa_project_concavity(theta)
  theta_list[[env$tau + 1L]] <- theta

  # ── Backward sweep ─────────────────────────────────────────────────────────
  for (t in seq(env$tau, 1L)) {

    # Step 1: evaluate separable approximation on full grid (cheap: O(nSdx*(nI+nJ)))
    V_next <- vfa_eval(theta, env)

    # Step 2: select states to solve LPs for
    state_idx <- if (use_full) integer(0L) else sort(sample(env$nSdx, K))

    # Step 3: Bellman update — LP solves only for state_idx (or all if empty)
    step <- bellmanUpdatePtr(
      problem_ptr    = prob,
      t              = t - 1L,
      stateSupport   = as.double(env$stateSupport),
      flowSupport    = as.double(env$Q$vals),
      scnpb          = env$scnpb,
      alpha          = env$alpha,
      V_next         = V_next,
      numThreads     = numThreads,
      traversalOrder = traversalOrder,
      chunkSize      = chunkSize,
      stateSubset    = state_idx   # empty = all states; 1-based
    )

    Q      [t, ] <- step$Q_t
    pi_star[t, ] <- step$pi_star_t
    pi_rand[t, ] <- step$pi_rand_t

    # Step 4: fit new theta from Bellman targets (NA for unsampled states ignored)
    theta <- vfa_fit(step$V_t, env)

    # Step 5: optional concavity projection
    if (project_concavity) theta <- vfa_project_concavity(theta)

    theta_list[[t]]  <- theta
    V_approx[t, ]    <- vfa_eval(theta, env)
  }

  list(
    V_approx   = V_approx,
    theta_list = theta_list,
    Q          = Q,
    pi_star    = pi_star,
    pi_rand    = pi_rand
  )
}
