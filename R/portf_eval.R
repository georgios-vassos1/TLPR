#' Evaluating a portfolio contract without state-dependent constraints
#' 
#' This function evaluates a portfolio contract without state-dependent constraints.
#' 
#' @param env The environment containing the parameters of the portfolio contract.
#' @param A_ A vector containing the values of the assets.
#' @param plot Logical indicating whether to plot the cost.
#' 
#' @return A list containing the cost and allocation of the portfolio contract.
#' 
#' @details This function evaluates a portfolio contract without considering state-dependent constraints. 
#' It computes the cost and allocation for each asset value provided in the vector \code{A_}. If \code{plot} is \code{TRUE}, it plots the cost.
#' 
#' @export
eval_stateless_portfolio <- function(env, A_, plot = FALSE) {

  cost <- rep(NaN, env$tau * nA)
  allocation <- rep(NaN, env$tau * nA * env$nvars)
  nA <- length(A_)

  offdx <- cumsum(c(env$nCS+env$nCO, 1L))

  rhs_ <- c(rep(NA, offdx[1L]), NA)
  for (t in seq(env$tau)) {
    obj_ <- c(env$CTb, env$CTo[t,])
    rhs_[seq(offdx[1L])] <- c(env$Cb[t,], env$Co[t,])
    for (j in seq(nA)) {
      rhs_[offdx[2L]] <- A_[j]
      optx <- optimal_assignment(model, obj_, rhs_)
      if (optx$status == "OPTIMAL") {
        cost[(t-1L)*nA+j] <- optx$objval
        allocation[((t-1L)*nA+(j-1L))*env$nvars+seq(env$nvars)] <- optx$x
      } else {
        print(optx$status)
      }
    }
  }

  if (plot) {
    matplot(A_, matrix(cost, nrow = nA), type = 'l')
  }

  list(
    "cost" = cost,
    "assignment" = allocation
  )
}

#' Evaluating a portfolio contract with state-dependent constraints
#' 
#' This function evaluates a portfolio contract considering state-dependent constraints.
#' 
#' @param env The environment containing the parameters of the portfolio contract.
#' @param A_ A vector containing the values of the assets.
#' @param S0 Initial state vector.
#' @param Q Vector representing inflow constraints.
#' @param D Vector representing outflow constraints.
#' @param ... Additional arguments to be passed.
#' 
#' @return A list containing the states, cost, and allocation of the portfolio contract.
#' 
#' @details This function evaluates a portfolio contract considering state-dependent constraints. 
#' It computes the cost and allocation for each asset value provided in the vector \code{A_} and based on the given state constraints. 
#' It also updates the state vector according to the transactions made during each period.
#' 
#' @export
eval_portfolio <- function(env, A_, S0, Q, D, ...) {

  nA <- length(A_)

  offdx <- cumsum(c(env$nCS+env$nCO, env$nI, env$nJ, 1L))
  rhs_  <- rep(NA, offdx[4L])

  cost <- vector(mode = "list", length = env$tau)
  allocation <- vector(mode = "list", length = env$tau)
  S_ <- vector(mode = "list", length = env$tau+1L)
  S_[[1L]] <- S0

  for (t in seq(env$tau)) {
    m <- nrow(S_[[t]])
    S_[[t+1L]] <- matrix(0L, nrow = m*nA, ncol = env$nL)
    cost[[t]]  <- matrix(NA, nrow = m, ncol = nA)
    allocation[[t]] <- numeric(m*nA*env$nvars)
    obj_ <- c(env$CTb, env$CTo[t,])
    rhs_[seq(offdx[1L])] <- c(env$Cb[t,], env$Co[t,])
    for (i in seq(m)) {
      rhs_[offdx[1L]+env$I_] <- env$R - S_[[t]][i,env$nI+env$J_]
      rhs_[offdx[2L]+env$J_] <- S_[[t]][i,env$I_] + Q[(t-1L)*env$nI+env$I_]
      for (j in seq(nA)) {
        rhs_[offdx[4L]] <- A_[j]
        optx <- optimal_assignment(model, obj_, rhs_)
        if (optx$status == "OPTIMAL") {
          cost[[t]][i,j] <- optx$objval
          allocation[[t]][((i-1L)*nA+(j-1L))*env$nvars+seq(env$nvars)] <- optx$x
        } else {
          print(optx$status)
        }
        # Allocation
        X.I <- unlist(lapply(env$to_j, function(k) sum(optx$x[k])))
        X.J <- unlist(lapply(env$from_i, function(k) sum(optx$x[k])))
        # Transition
        S_[[t+1L]][(i-1L)*nA+j, env$I_] <- S_[[t]][i, env$I_] + Q[(t-1L)*env$nI+env$I_] - X.I
        S_[[t+1L]][(i-1L)*nA+j, env$nI+env$J_] <- S_[[t]][i,env$nI+env$J_] - D[(t-1L)*env$nJ+env$J_] + X.J
      }
    }
  }

  list(
    "states" = S_,
    "cost"   = cost,
    "assignment" = allocation
  )
}

#' Single-stage portfolio contract evaluation with state-dependent constraints
#' 
#' This function evaluates a single-stage portfolio contract considering state-dependent constraints.
#' 
#' @param env The environment containing the parameters of the portfolio contract.
#' @param A_ A vector containing the values of the assets.
#' @param t The time period for which the evaluation is performed.
#' @param S.t The state vector at time t.
#' @param Q.t The inflow constraint at time t.
#' @param ... Additional arguments to be passed.
#' 
#' @return A list containing the cost and allocation of the portfolio contract for the specified time period.
#' 
#' @details This function evaluates a single-stage portfolio contract considering state-dependent constraints for the specified time period. 
#' It computes the cost and allocation for each asset value provided in the vector \code{A_} based on the given state at time \code{t} and 
#' the associated constraints.
#' 
#' @export
eval_portfolio_t <- function(env, A_, t, S.t, Q.t, ...) {

  nA <- length(A_)
  offdx <- cumsum(c(env$nCS+env$nCO, env$nI, env$nJ, 1L))
  rhs_  <- rep(NA, offdx[4L])

  cost <- rep(NaN, nA)
  allocation <- numeric(nA*env$nvars)

  obj_ <- c(env$CTb, env$CTo[t,])
  rhs_[seq(offdx[1L])]   <- c(env$Cb[t,], env$Co[t,])
  rhs_[offdx[1L]+env$I_] <- env$R - S.t[env$nI+env$J_]
  rhs_[offdx[2L]+env$J_] <- S.t[env$I_] + Q.t
  
  for (j in seq(nA)) {
    rhs_[offdx[4L]] <- A_[j]
    optx <- optimal_assignment(model, obj_, rhs_)
    if (optx$status == "OPTIMAL") {
      cost[j] <- optx$objval
      allocation[(j-1L)*env$nvars+seq(env$nvars)] <- optx$x
    } else {
      print(optx$status)
    }
  }

  list(
    "cost" = cost,
    "assignment" = allocation
  )
}
