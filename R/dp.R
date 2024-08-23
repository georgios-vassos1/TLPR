#' Generate Adjustment Weights
#'
#' This function generates a vector of adjustment weights based on the input environment `env`. 
#' The weights are computed using powers of `nSI` and `nSJ` (if `nJ` is greater than 1).
#'
#' @param env An environment containing the variables `nSI`, `nI`, `nJ`, and `nSJ`. 
#'   - `nSI` is an integer representing a scaling factor for the first group.
#'   - `nI` is an integer indicating the number of elements in the first group.
#'   - `nJ` is an integer representing the number of elements in the second group.
#'   - `nSJ` is an integer representing a scaling factor for the second group.
#'
#' @return A numeric vector containing the adjustment weights.
#'
#' @export
get_adjustment_weights <- function(env) {
  # Generate powers of nSI from 1 to nI
  si_powers <- with(env, nSI ^ (1L:nI))
  
  if (env$nJ == 1L) return(c(1L, si_powers))
  # Calculate the base value (nSI ^ nI)
  base_value <- with(env, nSI ^ nI)
  # Generate products of base_value and powers of nSJ from 0 to (nJ - 1)
  sj_powers <- base_value * with(env, nSJ ^ (1L:(nJ - 1L)))
  # Combine the results
  c(1L, si_powers, sj_powers)
}

#' Compute Environment Transitions and Rewards
#'
#' This function computes the transition matrix and rewards for a given environment over a specified time horizon.
#' The function iterates over possible states, actions, and scenarios, computing the optimal assignments and 
#' determining the resulting state transitions and associated rewards.
#'
#' @param env A list representing the environment, which includes all necessary data structures such as state, action, and scenario spaces.
#' @param model An optimization model object used for solving the assignment problem.
#' @param t An integer representing the current time step.
#' @param start An integer representing the starting index of the state space to be considered.
#' @param end An integer representing the ending index of the state space to be considered.
#'
#' @return A matrix where each row represents a transition, containing:
#' \describe{
#'   \item{next_i}{Index of the next state of the system.}
#'   \item{i}{The current state index.}
#'   \item{j}{The action index.}
#'   \item{kdx}{The scenario index, a combination of Q, D, and W indices.}
#' }
#'
#' @details
#' The function works by iterating over the state indices `i`, action indices `j`, and scenario indices (combinations of `Q`, `D`, and `W`).
#' For each combination, it computes the optimal assignment using the provided optimization model, and based on the results, 
#' it computes the index of the next state and the reward. If the optimization problem is infeasible, the corresponding transition is skipped.
#'
#' @examples
#' \dontrun{
#' result <- computeEnvironmentRx(env, model, t = 1, start = 1, end = 100)
#' }
#'
#' @export
computeEnvironmentRx <- function(env, model, t, start, end) {

  transit <- matrix(NA, nrow = (end - start + 1L) * env$nAdx * env$nScen, ncol = 6L)
  for (i in seq(start, end)) {
    for (j in seq(env$nAdx)) {
      for (k1 in seq(env$nQdx)) {
        q <- env$Q$vals[env$Qdx[k1,]]

        rhs <- c(
          # Carrier capacity
          env$Cb[t,], env$Co[t,], 
          # Storage limits
          env$R - env$extendedStateSupport[env$Sdx[i,env$nI+env$J_]], env$stateSupport[env$Sdx[i,env$I_]] + q, 
          # Transport volume
          env$actionSupport[j])

        for (k3 in seq(env$nWdx)) {
          w <- rep(env$W$vals[env$Wdx[k3,]], env$nL)

          optx <- optimal_assignment(
            model, obj_ = c(env$CTb, w), 
            rhs_ = rhs)

          if (optx$status == "INFEASIBLE") next

          # Obtain origin-destination assignment volumes
          xI <- unlist(lapply(env$from_i, function(l) sum(optx$x[l])))
          xJ <- unlist(lapply(env$to_j,   function(l) sum(optx$x[l])))

          for (k2 in seq(env$nDdx)) {
            d <- env$D$vals[env$Ddx[k2,]]

            # Compute the index of the next state of the system
            next_i <- sum(
              c(
                pmax(pmin(env$stateSupport[env$Sdx[i,env$I_]] + q - xI, env$R), 0L), 
                pmin(pmax(env$extendedStateSupport[env$Sdx[i,env$nI+env$J_]] - d + xJ, -env$R), env$R) + env$R
              ) * env$stateKeys) + 1L

            kdx <- ((k3 - 1L) * env$nQdx + (k2 - 1L)) * env$nDdx + k1
            # kdx <- sum(c(Qdx[k1,] - 1L, Ddx[k2,] - 1L, Wdx[k3,] - 1L) * flowKeys) + 1L

            # Store into the transition matrix
            transit[((i - start) * env$nAdx + (j - 1L)) * env$nScen + kdx,] <- c(
              next_i,
              optx$objval, 
              i, j, kdx, t)
          }
        }
      }
    }
  }

  transit
}

#' Endogenous state vector
#'
#' @param env Environment to host dynamic programming data
#' @param nI Number of origin hubs
#' @param nJ Number of destination hubs
#' @param max_ Maximum value
#' @param incr_ Increment value
#' @export
get_state_indices <- function(env, nI, nJ, max_, incr_ = 1L) {
  env$SI_ <- seq(0.0, max_, by = incr_)
  env$SJ_ <- seq(-max_, max_, by = incr_)
  env$nSI <- length(env$SI_)
  env$nSJ <- length(env$SJ_)

  env$nSdx <- (env$nSI^nI)*(env$nSJ^nJ)
  env$Sdx  <- consolidate_idx(c(replicate(nI, seq(env$nSI), simplify = FALSE), replicate(nJ, seq(env$nSJ), simplify = FALSE)))
}

#' Get scenario space
#'
#' @param Q Q parameter
#' @param D D parameter
#' @param W W parameter
#' @param hmg Boolean indicating whether to use hmg or not (default: TRUE)
#' @export
get_scenario_space <- function(env, nI, nJ, Q, D, W, hmg=TRUE) {
  nQ <- length(Q$vals)
  nD <- length(D$vals)
  nW <- length(W$vals)

  len_spot <- 1L
  if (!hmg) {
    len_spot <- env$nL
  }

  consolidate_idx(
    c(replicate(nI, seq(nQ), simplify = FALSE), 
      replicate(nJ, seq(nD), simplify = FALSE), 
      replicate(len_spot, seq(nW), simplify = FALSE))
  ) -> idx

  env$nOmega <- nrow(idx)

  I_ <- seq(nI)
  J_ <- seq(nJ)

  env$Phi.t <- cbind(
    matrix(Q$vals[idx[,I_]], ncol = nI), 
    matrix(D$vals[idx[,nI+J_]], ncol = nJ), 
    matrix(W$vals[idx[,nI+nJ+seq(len_spot)]], ncol = len_spot))

  env$PPhi.t <- apply(cbind(
    matrix(Q$prob[idx[,I_]], ncol = nI), 
    matrix(D$prob[idx[,nI+J_]], ncol = nJ), 
    matrix(W$prob[idx[,nI+nJ+seq(len_spot)]], ncol = len_spot)), 1L, prod)
}

#' Convert state variable to state index
#' 
#' @param env Environment object containing relevant data
#' @param Si Value representing some state information
#' @param Sj Value representing some state information
#' @param max.S Maximum value for a state
#' @param weights Numeric vector of weights
#' @return A numeric vector representing state index
stateIdx <- function(env, Si, Sj, max.S, weights, ...) {
  # Example weights: c(1L, nSI, nSI^2L, (nSI^2L)*nSJ)
  c(c(Si, Sj+max.S) %*% weights + 1L)
}

#' Run scenarios
#' 
#' @param env Environment object containing relevant data
#' @param env.dp Environment object containing dynamic programming data
#' @param t Value representing time
#' @param start Start index
#' @param end End index
#' @param scenaria A vector of scenario indexes
#' @param weights Numeric vector of weights
#' @param FUN The assignment algorithm
#' @return A matrix representing transit
#' @export
compute_environment <- function(env, env.dp, t, start, end, scenaria, weights, FUN, ...) {
  args  <- list(...)
  alpha <- args[["alpha"]]
  edx   <- nrow(env.dp$model$A)

  nScen   <- length(scenaria)
  transit <- matrix(NA, nrow = (end-start+1L)*env.dp$nA*nScen, ncol = 5L)
  for (i in seq(start, end)) {
    if(!(i%%100L)) print(i)
    for (j in seq(env.dp$nA)) {
      for (k in seq(nScen)) {
        # Exogenous variables
        scndx <- scenaria[k]
        q <- env.dp$Phi.t[scndx,env$I_]
        d <- env.dp$Phi.t[scndx,env$nI+env$J_]
        w <- rep(env.dp$Phi.t[scndx,env$nI+env$nJ+1L], env$nL)

        rhs <- c(
          # Carrier capacity
          env$Cb[t,], env$Co[t,], 
          # Storage limits
          env$R - env.dp$SJ_[env.dp$Sdx[i,env$nI+env$J_]], env.dp$SI_[env.dp$Sdx[i,env$I_]] + q, 
          # Transport volume
          env.dp$A_[j])

        optx <- FUN(
          env.dp$model, obj_ = c(env$CTb, w), 
          rhs_ = rhs, k = env$nvars, edx = edx)

        if (optx$status == "INFEASIBLE") next

        # Obtain origin-destination assignment volumes
        xI <- unlist(lapply(env$from_i, function(l) sum(optx$x[l])))
        xJ <- unlist(lapply(env$to_j,   function(l) sum(optx$x[l])))

        # Compute the index of the next state of the system
        next_i <- stateIdx(
          env, 
          pmax(pmin(env.dp$SI_[env.dp$Sdx[i,env$I_]] + q, env$R) - xI, 0L), 
          pmin(pmax(env.dp$SJ_[env.dp$Sdx[i,env$nI+env$J_]] - d, -env$R) + xJ, env$R), 
          env.dp$max.S, weights)

        # Store into the transition matrix
        transit[((i-start)*env.dp$nA+(j-1L))*nScen+k,] <- c(
          next_i,
          h.t(env, env.dp$SI_[env.dp$Sdx[next_i,env$I_]], env.dp$SJ_[env.dp$Sdx[next_i,env$nI+env$J_]], alpha = alpha) + optx$objval, 
          i, j, scndx)
      }
    }
  }
  transit
}
