#' Reverse the columns of a matrix
#'
#' This function reverses the order of columns in a matrix.
#'
#' @param X A matrix.
#' @return A matrix with columns in reversed order.
#' @examples
#' mat <- matrix(1:9, ncol = 3)
#' rev_cols(mat)
rev_cols <- function(X) {
  X[, rev(seq(ncol(X)))]
}

#' Consolidate indices from lists
#'
#' This function takes a list of indices and consolidates them into a matrix where each row contains a unique combination of indices.
#'
#' @param lists A list of indices.
#' @return A matrix where each row represents a unique combination of indices.
#' @examples
#' lists <- list(c(1, 2), c(3, 4), c(5, 6))
#' consolidate_idx__(lists)
consolidate_idx__ <- function(lists) {
  lists |>
    expand.grid() |>
    unname() |>
    as.matrix()
}

#' Consolidate indices from lists using data.table::CJ
#'
#' This function takes a list of indices and consolidates them into a matrix where each row contains a unique combination of indices using data.table::CJ function.
#'
#' @param lists A list of indices.
#' @return A matrix where each row represents a unique combination of indices.
#' @examples
#' lists <- list(c(1, 2), c(3, 4), c(5, 6))
#' consolidate_idx(lists)
#' @export
consolidate_idx <- function(lists) {
  do.call(data.table::CJ, rev(lists)) |>
    unname() |>
    as.matrix() |>
    rev_cols()
}

#' Endogenous state vector
#'
#' @param env Environment to host dynamic programming data
#' @param nI Number of origin hubs
#' @param nJ Number of destination hubs
#' @param max_ Maximum value
#' @param incr_ Increment value
#' @export
get_state_indices <- function(env, nI, nJ, max_, incr_) {
  env$SI_ <- seq(0L, max_)
  env$SJ_ <- seq(-max_, max_)
  env$nSI <- length(SI_)
  env$nSJ <- length(SJ_)

  env$nSdx <- (nSI^nI)*(nSJ^nJ)
  env$Sdx  <- consolidate_idx(c(replicate(nI, seq(nSI), simplify = FALSE), replicate(nJ, seq(nSJ), simplify = FALSE)))
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
          max.S, weights)

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
