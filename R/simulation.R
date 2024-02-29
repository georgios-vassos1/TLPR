#' Get From Routes
#' 
#' Retrieves the routes originating from different sources.
#' 
#' @param env The environment variable.
#' @export
get_from_routes <- function(env, ...) {
  env$from_i <- vector(mode = 'list', length = env$nI)
  for (i in seq(env$nI)) {
    # All possible routes from origin i to the destinations
    idx <- (i-1L)*env$nJ + seq(env$nJ)
    # All indices in the bids that start from origin i
    msk <- which(apply(outer(env$L_, idx, '=='), 1L, any))
    # Adding spot indices that start from origin i
    idx <- env$nL_ + c(outer(idx, (seq(env$nCO)-1L)*env$nL, '+'))
    # Store result
    env$from_i[[i]] <- c(msk, idx)
  }
}

#' Get To Routes
#' 
#' Retrieves the routes leading to different destinations.
#' 
#' @param env The environment variable.
#' @export
get_to_routes <- function(env, ...) {
  env$to_j <- vector(mode = 'list', length = env$nJ)
  for (j in seq(env$nJ)) {
    # All possible routes from all origins to destination j
    idx <- (seq(env$nJ)-1L)*env$nI
    # All indices in the bids that go to destination j
    msk <- which(apply(outer(env$L_, idx + j, '=='), 1L, any))
    # Adding spot indices that go to destination j
    idx <- env$nL_ + c(outer(idx, (seq(env$nCO) - 1L)*env$nL, '+')) + j
    # Store result
    env$to_j[[j]] <- c(msk, idx)
  }
}

#' Sample Exogenous State
#' 
#' Samples the exogenous state based on a given function and parameters.
#' 
#' @param func The function for sampling.
#' @param params The parameters for the sampling function.
sample_exogenous_state <- function(func, params, ...) {
  do.call(func, params)
}

#' Simulate System
#' 
#' Simulates the system based on policies, exogenous states, and environment parameters.
#' 
#' @param env The environment variable.
#' @param policy The policy for simulation.
#' @param args Arguments for policy functions.
#' @param exog Exogenous states for simulation.
#' @export
simulate_system <- function(env, policy, args, exog, ...) {
  npi <- length(policy)
  # Simulation metrics
  S.I  <- matrix(0.0, nrow = (env$tau+1L)*env$nI, ncol = npi)
  S.J  <- matrix(0.0, nrow = (env$tau+1L)*env$nJ, ncol = npi)
  X.I  <- matrix(0.0, nrow = env$tau*env$nI, ncol = npi)
  X.J  <- matrix(0.0, nrow = env$tau*env$nJ, ncol = npi)
  cost <- matrix(NA, nrow = env$tau, ncol = npi)
  q    <- numeric(env$tau)
  # Exogenous states
  Q <- sample_exogenous_state(exog$Q$func, exog$Q$params)
  D <- sample_exogenous_state(exog$D$func, exog$D$params)
  # Simulation
  for (t in seq(env$tau)) {
    idx <- (t-1L)*env$nI+env$I_
    jdx <- (t-1L)*env$nJ+env$J_

    q[t] <- round(mean(S.I[idx,] + Q[idx]))

    args[["obj_"]] <- c(env$CTb, env$CTo[t,])
    args[["n"]]    <- q[t]

    for (pdx in seq_along(policy)) {
      args[["rhs_"]] <- c(env$Cb[t,], env$Co[t,], env$R - S.J[idx,pdx], S.I[idx,pdx] + Q[idx], q[t])

      optx <- do.call(policy[[pdx]], args)
      if (!is.null(optx[["status"]])) {
        if (optx[["status"]] != "OPTIMAL") {
          print(optx[["status"]])
          return(list("status" = FALSE))
        }
      }
      a.t  <- optx$x
      cost[t,pdx] <- h.t(env, S.I[idx,pdx], S.J[jdx,pdx], env$alpha) + optx$objval

      X.I[idx,pdx] <- unlist(lapply(env$from_i, function(k) sum(a.t[k])))
      X.J[jdx,pdx] <- unlist(lapply(env$to_j,   function(k) sum(a.t[k])))

      S.I[idx+env$nI,pdx] <- pmin(pmax(S.I[idx,pdx] + Q[idx] - X.I[idx,pdx], 0L), env$R)
      S.J[jdx+env$nJ,pdx] <- pmin(S.J[jdx,pdx] - D[jdx] + X.J[jdx,pdx], env$R)
    }
  }
  list(
    "cost" = cost, "S.I" = S.I, "S.J" = S.J, "X.I" = X.I, "X.J" = X.J, "Q" = Q, "D" = D, "status" = TRUE
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
