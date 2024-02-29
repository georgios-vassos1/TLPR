#' Endogenous state vector
#'
#' @param env Environment parameters
#' @param max_ Maximum value
#' @param incr_ Increment value
#' @export
get_state_indices <- function(env, max_, incr_) {
  SI_ <<- seq(0L, max_)
  SJ_ <<- seq(-max_, max_)
  nSI <<- length(SI_)
  nSJ <<- length(SJ_)

  nSdx <<- (nSI^env$nI)*(nSJ^env$nJ)
  Sdx  <<- c(
    replicate(env$nI, seq(nSI), simplify = FALSE),
    replicate(env$nJ, seq(nSJ), simplify = FALSE)) |>
    expand.grid() |>
    unname() |>
    as.matrix()
}

#' Get scenario space
#'
#' @param Q Q parameter
#' @param D D parameter
#' @param W W parameter
#' @param hmg Boolean indicating whether to use hmg or not (default: TRUE)
#' @export
get_scenario_space <- function(Q, D, W, hmg=TRUE) {
  nQ <- length(Q$vals)
  nD <- length(D$vals)
  nW <- length(W$vals)
  
  len_spot <- 1L
  if (!hmg) {
    len_spot <- env$nL
  }
  
  c(replicate(env$nI, seq(nQ), simplify = FALSE), 
    replicate(env$nJ, seq(nD), simplify = FALSE), 
    replicate(len_spot, seq(nW), simplify = FALSE)) |>
    expand.grid() |>
    unname() |>
    as.matrix() -> idx
  
  nOmega <<- nrow(idx)
  
  Phi.t <<- cbind(
    matrix(Q$vals[idx[,env$I_]], ncol = env$nI), 
    matrix(D$vals[idx[,env$nI+env$J_]], ncol = env$nJ), 
    matrix(W$vals[idx[,env$nI+env$nJ+seq(len_spot)]], ncol = len_spot))
  
  PPhi.t <<- apply(cbind(
    matrix(Q$prob[idx[,env$I_]], ncol = env$nI), 
    matrix(D$prob[idx[,env$nI+env$J_]], ncol = env$nJ), 
    matrix(W$prob[idx[,env$nI+env$nJ+seq(len_spot)]], ncol = len_spot)), 1L, prod)
}

#' All possible allocations of n containers to k routes
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
#' @param t Value representing time
#' @param start Start index
#' @param end End index
#' @param scenaria A vector of scenarios
#' @param weights Numeric vector of weights
#' @return A matrix representing transit
#' @export
run_scenarios <- function(env, t, start, end, scenaria, weights) {
  nScen   <- length(scenaria)
  transit <- matrix(NA, nrow = (end-start+1L)*nA*nScen, ncol = 5L)
  ldx <- 1L
  for (i in seq(start, end)) {
    if(!(i%%10L)) print(i)
    for (j in seq(nA)) {
      for (k in seq(nScen)) {
        scndx <- scenaria[k]
        q <- Phi.t[scndx,env$I_]
        d <- Phi.t[scndx,env$nI+env$J_]
        w <- rep(Phi.t[scndx,env$nI+env$nJ+1L], env$nL)

        optx <- optimal_assignment(
          model, obj_ = c(env$CTb, w), 
          rhs_ = c(env$Cb[t,], env$Co[t,], env$R - SJ_[Sdx[i,env$nI+env$J_]], SI_[Sdx[i,env$I_]] + q, A_[j]))
        if (optx$status == "INFEASIBLE") next

        xI <- unlist(lapply(env$from_i, function(l) sum(optx$x[l])))
        xJ <- unlist(lapply(env$to_j,   function(l) sum(optx$x[l])))

        transit[(ldx-1L)*nA*nScen+(j-1L)*nScen+k,] <- c(stateIdx(
          env,
          pmin(pmax(SI_[Sdx[i,env$I_]] + q - xI, 0L), env$R),
          pmin(pmax(SJ_[Sdx[i,env$nI+env$J_]] - d + xJ, -env$R), env$R),
          max.S, weights),
          h.t(env, SI_[Sdx[i,env$I_]], SJ_[Sdx[i,env$nI+env$J_]]) + optx$objval, i, j, scndx)
      }
    }
    ldx <- ldx + 1L
  }
  transit
}
