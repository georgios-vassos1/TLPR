#' Sample Exogenous State
#' 
#' Samples the exogenous state based on a given function and parameters.
#' 
#' @param func The function for sampling.
#' @param params The parameters for the sampling function.
sample_exogenous_state <- function(func, params, ...) {
  do.call(func, params)
}

#' Random Volume Selection
#' 
#' Selects a random volume based on the given capacity and limits.
#' 
#' @param capacity A numeric vector representing the available capacities.
#' @param limit A numeric value representing the upper limit for the selected volume.
#' @param ... Additional arguments (currently not used).
#' 
#' @return A numeric value representing the selected volume, which is the minimum of a random sample from the capacity and the provided limit.
#' @export
random_volume <- function(capacity, limit, ...) {
  min(sample(capacity, 1L), limit)
}

#' Simulate System
#' 
#' Simulates the system based on policy functions for transition and allocation, exogenous states, and environment parameters.
#' 
#' @param env A list containing the environment variables including time horizon, resource limits, and cost parameters.
#' @param pi_trans A function representing the policy for the transition.
#' @param pi_alloc A list of functions representing the allocation policies to be used during the simulation.
#' @param args A list of additional arguments to be passed to the policy functions.
#' @param Q A matrix representing the exogenous state related to supply (optional). If NULL, it will be sampled based on the exog parameter.
#' @param D A matrix representing the exogenous state related to demand (optional). If NULL, it will be sampled based on the exog parameter.
#' @param ... Additional arguments (not used in the function directly).
#' @param exog A list of functions and parameters for sampling exogenous states if Q or D are not provided.
#' 
#' @return A list containing the simulation results, including costs, allocations, state matrices, and the status of the simulation.
#' @export
simulate_system <- function(env, pi_trans, pi_alloc, args, Q=NULL, D=NULL, ...) {
  npi <- length(pi_alloc)
  # Simulation metrics
  S.I  <- matrix(0.0, nrow = (env$tau + 1L) * env$nI, ncol = npi)
  S.J  <- matrix(0.0, nrow = (env$tau + 1L) * env$nJ, ncol = npi)
  X.I  <- matrix(0.0, nrow = env$tau * env$nI, ncol = npi)
  X.J  <- matrix(0.0, nrow = env$tau * env$nJ, ncol = npi)
  cost <- matrix(NA, nrow = env$tau, ncol = npi)
  allocation <- matrix(NA, env$tau * env$nvars, ncol = npi)
  q    <- numeric(env$tau)

  # Exogenous state variables
  if (is.null(Q)) {
    Q <- sample_exogenous_state(exog$Q$func, exog$Q$params)
  }
  if (is.null(D)) {
    D <- sample_exogenous_state(exog$D$func, exog$D$params)
  }

  # Simulation
  for (t in seq(env$tau)) {
    idx <- (t-1L)*env$nI+env$I_
    jdx <- (t-1L)*env$nJ+env$J_

    # Transportation policy
    args[["capacity"]] <- sum(env$Cb[t,])+sum(env$Co[t,])
    args[["limit"]]    <- colSums(S.I[idx,,drop=F]+Q[idx])

    q[t] <- do.call(pi_trans, args)

    # Transition mechanism
    args[["obj_"]] <- c(env$CTb, env$CTo[t,])
    args[["n"]]    <- q[t]
    args[["x"]]    <- c(env$Cb[t,], env$Co[t,])

    for (pdx in seq_along(pi_alloc)) {
      args[["rhs_"]] <- c(env$Cb[t,], env$Co[t,], env$R-S.J[jdx,pdx], S.I[idx,pdx]+Q[idx], q[t])

      # Optimize for policy index pdx
      optx <- do.call(pi_alloc[[pdx]], args)
      if (!is.null(optx[["status"]])) {
        if (optx[["status"]] == "INFEASIBLE") {
          return(list("status" = FALSE))
        }
      }
      # Store optimal cost
      cost[t,pdx] <- optx$objval + h.t(env, S.I[idx,pdx], S.J[jdx,pdx], env$alpha)
      # Store optimal allocation
      a.t <- optx$x
      allocation[(t-1L)*env$nvars+seq(env$nvars),pdx] <- optx$x
      # Store optimal decisions
      X.I[idx,pdx] <- unlist(lapply(env$from_i, function(k) sum(a.t[k])))
      X.J[jdx,pdx] <- unlist(lapply(env$to_j,   function(k) sum(a.t[k])))
      # Update state variables
      S.I[idx+env$nI,pdx] <- pmin(pmax(S.I[idx,pdx] + Q[idx] - X.I[idx,pdx], 0L), env$R)
      S.J[jdx+env$nJ,pdx] <- pmax(pmin(S.J[jdx,pdx] - D[jdx] + X.J[jdx,pdx], env$R), -env$R)
    }
  }
  # Return list
  list(
    "cost" = cost, "allocation" = allocation, 
    "S.I" = S.I, "S.J" = S.J, "X.I" = X.I, "X.J" = X.J, "Q" = Q, "D" = D, 
    "status" = TRUE
  )
}

#' Simulate the transition of states and policies in a dynamic system
#'
#' @description
#' The `system_transition` function simulates the transitions of states in a dynamic system based on given policies, transition probabilities, and initial states. The function calculates the system's state transitions, the resulting actions, and the associated costs over a specified time horizon.
#'
#' @param env A list containing environment parameters such as time horizon `tau`, number of states `nI`, `nJ`, state keys, action support, state support, and scenario indices.
#' @param transit A matrix representing the state transition probabilities and associated costs for each action.
#' @param pi A matrix representing the policy probabilities for each state-action pair over time.
#' @param varphidx A vector of scenario indices representing different possible outcomes at each time step.
#' @param init_s An optional initial state vector. If `NULL`, the initial state is assumed to be zero.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{S.I}{A matrix representing the states of type I over time.}
#'   \item{S.J}{A matrix representing the states of type J over time.}
#'   \item{X.I}{A matrix representing the differences in states of type I between time steps.}
#'   \item{X.J}{A matrix representing the differences in states of type J between time steps.}
#'   \item{cost}{The total cost incurred over the time horizon.}
#' }
#' 
#' @details
#' The function iteratively updates the state matrices `S.I` and `S.J` based on the given transition probabilities and policies. It computes the costs at each time step and accumulates them to return the total cost.
#'
#' @export
system_transition <- function(env, transit, pi, varphidx, init_s = NULL) {
  S.I  <- matrix(0.0, nrow = (env$tau + 1L) * env$nI, ncol = 1L)
  S.J  <- matrix(0.0, nrow = (env$tau + 1L) * env$nJ, ncol = 1L)
  X.I  <- matrix(0.0, nrow = env$tau * env$nI, ncol = 1L)
  X.J  <- matrix(0.0, nrow = env$tau * env$nJ, ncol = 1L)
  cost <- matrix(NA,  nrow = env$tau + 1L, ncol = 1L)
  q    <- numeric(env$tau)

  if (!is.null(init_s)) {
    S.I[1L,] <- init_s[env$I_]
    S.J[1L,] <- init_s[env$nI + env$J_]
  }

  # Simulation loop
  for (t in seq(env$tau)) {
    idx <- (t - 1L) * env$nI + env$I_
    jdx <- (t - 1L) * env$nJ + env$J_

    i <- sum(c(S.I[idx, 1L], env$R + S.J[jdx, 1L]) * env$stateKeys) + 1L
    j <- sample(seq(env$nAdx), 1L, prob = pi[t, (i - 1L) * env$nAdx + seq(env$nAdx)])
    k <- varphidx[t]

    p <- (((t - 1L) * env$nSdx + (i - 1L)) * env$nAdx + (j - 1L)) * env$nScen + k

    q[t] <- env$actionSupport[j]
    cost[t] <- transit[p, 2L] + h.t(env, S.I[idx, 1L], S.J[jdx, 1L], env$alpha)
    next_i  <- transit[p, 1L]

    S.I[env$nI + idx, 1L] <- env$stateSupport[env$Sdx[next_i, env$I_]]
    S.J[env$nJ + jdx, 1L] <- env$extendedStateSupport[env$Sdx[next_i, env$nI + env$J_]]

    X.I[idx, 1L] <- S.I[idx, 1L] + env$Q$vals[env$scndx[k, env$I_]] - S.I[idx + env$nI, 1L]
    X.J[jdx, 1L] <- S.J[jdx + env$nJ, 1L] - S.J[jdx, 1L] + env$D$vals[env$scndx[k, env$nI + env$J_]]
  }
  idx <- t * env$nI + env$I_
  jdx <- t * env$nJ + env$J_
  cost[t+1L] <- c(cbind(S.I[idx, 1L], pmax(S.J[jdx, 1L], 0L), - pmin(S.J[jdx, 1L], 0L)) %*% env$alpha)

  list(
    "S.I" = S.I,
    "S.J" = S.J,
    "X.I" = X.I,
    "X.J" = X.J,
    "cost" = sum(cost)
  )
}

#' Compute Statistics
#' 
#' Computes statistics for a given set of data.
#' 
#' @param v The data.
#' @param N The number of observations.
#' @export
compute_stats <- function(v, N) {
  cumv <- apply(v, 2L, cumsum)
  cumv.avg <- apply(cumv, 1L, mean)
  cumv.cnf <- cumv.avg + (apply(cumv, 1L, sd) / sqrt(N)) %*% matrix(qnorm(0.975) * c(-1.0, 1.0), ncol = 2L)
  unname(cbind(cumv.avg, cumv.cnf))
}
