# ─────────────────────────────────────────────────────────────────────────────
# vfa.R — Separable Value Function Approximation (VFA) for TLPR
#
# Represents V(s) ≈ Σᵢ vᵢ(sᵢ) + Σⱼ vⱼ(sⱼ) where each vₖ is a 1-D
# piecewise-constant function on the discrete inventory support.
#
# theta layout:
#   theta$orig  — (R+1)   × nI  matrix; row r = vₖ(stateSupport[r])
#   theta$dest  — (2R+1)  × nJ  matrix; row r = vₖ(extendedStateSupport[r])
#
# Both matrices are indexed in the same 1-based convention as env$Sdx, so
# theta$orig[env$Sdx[s, k], k] is the value of the k-th origin function at
# the inventory level of state s, and analogously for destinations.
# ─────────────────────────────────────────────────────────────────────────────

#' Initialise a zero-valued separable VFA parameter object
#'
#' @param env Environment produced by \code{dp_config}.
#' @return A named list with elements \code{orig} and \code{dest}.
#' @export
vfa_init <- function(env) {
  list(
    orig = matrix(0.0, nrow = env$R + 1L,      ncol = env$nI),
    dest = matrix(0.0, nrow = 2L * env$R + 1L, ncol = env$nJ)
  )
}

#' Evaluate the separable VFA on the full state grid
#'
#' Computes \eqn{\hat{V}(s) = \sum_k v_k(s_k)} for every state
#' \eqn{s \in \{1,\ldots,n_{Sdx}\}}.
#'
#' @param theta Parameter object returned by \code{vfa_init} or \code{vfa_fit}.
#' @param env Environment produced by \code{dp_config}.
#' @return Numeric vector of length \code{env$nSdx}.
#' @export
vfa_eval <- function(theta, env) {
  V <- numeric(env$nSdx)
  for (k in seq(env$nI))
    V <- V + theta$orig[env$Sdx[, k], k]
  for (k in seq(env$nJ))
    V <- V + theta$dest[env$Sdx[, env$nI + k], k]
  V
}

#' Fit a separable VFA from (state, target) pairs via OLS
#'
#' Builds a one-hot indicator feature matrix \eqn{\Phi} (one indicator column
#' per location-level combination) and solves the least-squares problem
#' \eqn{\min_\theta \|\Phi\theta - y\|^2}.  The design matrix is rank-deficient
#' by one constant (all rows sum to \eqn{n_I + n_J}); \code{lm.fit} handles
#' this gracefully via QR pivoting — any NA coefficients are set to zero.
#'
#' @param targets Numeric vector of length \code{env$nSdx}; NA entries are
#'   ignored (e.g. infeasible states).
#' @param env Environment produced by \code{dp_config}.
#' @param state_idx Integer vector of state indices that \code{targets}
#'   corresponds to.  Default: all states.
#' @return A \code{theta} parameter object.
#' @export
vfa_fit <- function(targets, env, state_idx = seq(env$nSdx)) {
  valid <- which(is.finite(targets))
  if (length(valid) == 0L) return(vfa_init(env))

  s <- state_idx[valid]
  y <- targets[valid]

  nR_orig <- env$R + 1L
  nR_dest <- 2L * env$R + 1L
  p       <- env$nI * nR_orig + env$nJ * nR_dest
  n       <- length(s)

  # Build one-hot indicator feature matrix (vectorised column assignment)
  Phi <- matrix(0.0, nrow = n, ncol = p)
  for (k in seq(env$nI)) {
    off  <- (k - 1L) * nR_orig
    Phi[cbind(seq(n), off + env$Sdx[s, k])] <- 1.0
  }
  for (k in seq(env$nJ)) {
    off  <- env$nI * nR_orig + (k - 1L) * nR_dest
    Phi[cbind(seq(n), off + env$Sdx[s, env$nI + k])] <- 1.0
  }

  coefs        <- lm.fit(Phi, y)$coefficients
  coefs[is.na(coefs)] <- 0.0

  theta <- vfa_init(env)
  for (k in seq(env$nI)) {
    off <- (k - 1L) * nR_orig
    theta$orig[, k] <- coefs[off + seq(nR_orig)]
  }
  for (k in seq(env$nJ)) {
    off <- env$nI * nR_orig + (k - 1L) * nR_dest
    theta$dest[, k] <- coefs[off + seq(nR_dest)]
  }
  theta
}

#' Pool-adjacent-violators for a non-increasing sequence
#'
#' Thin wrapper around base-R \code{isoreg}: negates input and output to
#' convert the non-decreasing isotonic regression into a non-increasing one.
#'
#' @param x Numeric vector.
#' @return Numeric vector of the same length — the closest non-increasing
#'   sequence to \code{x} in squared-error sense.
pava_nonincreasing <- function(x) -isoreg(-x)$yf

#' Project a 1-D value function onto concave piecewise-linear functions
#'
#' Enforces non-increasing finite differences
#' \eqn{v[r+1]-v[r] \leq v[r]-v[r-1]} via pool-adjacent-violators on the
#' slope vector, then reconstructs \eqn{v} preserving its mean.
#'
#' @param v Numeric vector representing \eqn{v(0), v(1), \ldots, v(n-1)}.
#' @return Projected vector of the same length.
concavity_project <- function(v) {
  if (length(v) < 3L) return(v)
  dv      <- diff(v)
  dv_proj <- pava_nonincreasing(dv)
  v_proj  <- cumsum(c(v[1L], dv_proj))
  v_proj + (mean(v) - mean(v_proj))   # preserve mean to avoid scale drift
}

#' Apply concavity projection to all 1-D functions in a theta object
#'
#' @param theta Parameter object from \code{vfa_init} or \code{vfa_fit}.
#' @return Updated \code{theta} with concave \code{orig} and \code{dest}
#'   columns.
#' @export
vfa_project_concavity <- function(theta) {
  theta$orig <- apply(theta$orig, 2L, concavity_project)
  theta$dest <- apply(theta$dest, 2L, concavity_project)
  theta
}
