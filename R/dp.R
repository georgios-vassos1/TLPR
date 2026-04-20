#' Configure DP Environment
#'
#' Populates \code{env} with all index arrays, cardinalities, and scenario
#' structures required by \code{computeEnvironmentRx}, \code{computeEnvironmentCx},
#' \code{dynamic_programming}, and \code{system_transition}.
#'
#' Fixes applied vs. the inline version previously found in demo scripts:
#' \itemize{
#'   \item \code{env$I_} and \code{env$J_} are set (\code{seq(nI)}, \code{seq(nJ)}).
#'   \item \code{env$nScen} is assigned directly (not via \code{with()}).
#' }
#'
#' @param env An R environment containing at minimum \code{R}, \code{nI},
#'   \code{nJ}, \code{nCO}, and \code{Q}, \code{D}, \code{W} lists with
#'   \code{$vals} and \code{$prob} entries.
#' @return \code{env}, invisibly, modified in place.
#' @export
dp_config <- function(env) {
  env$I_ <- seq(env$nI)
  env$J_ <- seq(env$nJ)

  env$stateSupport         <- 0L:env$R
  env$extendedStateSupport <- -env$R:env$R

  env$Sdx  <- do.call(TLPR::CartesianProductX, c(
    replicate(env$nI, seq_along(env$stateSupport),         simplify = FALSE),
    replicate(env$nJ, seq_along(env$extendedStateSupport), simplify = FALSE)))
  env$nSdx <- nrow(env$Sdx)

  env$actionSupport <- 0L:env$R
  env$Adx           <- seq_along(env$actionSupport)
  env$nAdx          <- length(env$Adx)

  env$nQ <- length(env$Q$vals)
  env$nD <- length(env$D$vals)
  env$nW <- length(env$W$vals)

  env$scndx <- do.call(TLPR::CartesianProductX, c(
    replicate(env$nI,  seq(env$nQ), simplify = FALSE),
    replicate(env$nJ,  seq(env$nD), simplify = FALSE),
    replicate(env$nCO, seq(env$nW), simplify = FALSE)))
  env$scnpb <- apply(env$scndx, 1L, function(x)
    prod(env$Q$prob[x[seq(env$nI)]],
         env$D$prob[x[env$nI + seq(env$nJ)]],
         env$W$prob[x[env$nI + env$nJ + seq(env$nCO)]]))

  env$Qdx  <- do.call(TLPR::CartesianProductX, replicate(env$nI,  seq(env$nQ), simplify = FALSE))
  env$Ddx  <- do.call(TLPR::CartesianProductX, replicate(env$nJ,  seq(env$nD), simplify = FALSE))
  env$Wdx  <- do.call(TLPR::CartesianProductX, replicate(env$nCO, seq(env$nW), simplify = FALSE))
  env$nQdx <- nrow(env$Qdx)
  env$nDdx <- nrow(env$Ddx)
  env$nWdx <- nrow(env$Wdx)
  env$nScen <- env$nQdx * env$nDdx * env$nWdx

  # Mixed-radix keys for encoding (Qdx, Ddx, Wdx) sub-indices into a flat
  # scenario index.  Used by computeEnvironmentRx to compute kdx.
  env$flowKeys <- c(
    env$nQ  ^ (seq(env$nI)  - 1L),
    env$nQdx * env$nD  ^ (seq(env$nJ)  - 1L),
    env$nQdx * env$nDdx * env$nW ^ (seq(env$nCO) - 1L))

  invisible(env)
}

#' Stochastic Backward Induction
#'
#' Runs backward induction over the pre-computed transit table, integrating
#' over all feasible scenarios using \code{env$scnpb} as weights.
#'
#' @param env Environment produced by \code{dp_config}.
#' @param transit Transit matrix (\code{tau*nSdx*nAdx*nScen x 6}) with columns
#'   \code{(next_i, objval, i, j, kdx, t)}.
#' @param ... Unused; reserved for future extension.
#' @return Named list: \code{V}, \code{Q}, \code{pi_star}, \code{pi_rand}.
#' @export
dynamic_programming <- function(env, transit, ...) {
  V       <- matrix(NA, nrow = env$tau + 1L, ncol = env$nSdx)
  Q       <- matrix(NA, nrow = env$tau,      ncol = env$nSdx * env$nAdx)
  pi_rand <- matrix(NA, nrow = env$tau,      ncol = env$nSdx * env$nAdx)
  pi_star <- matrix(NA, nrow = env$tau,      ncol = env$nSdx * env$nAdx)

  V[env$tau + 1L, ] <- -c(
    cbind(
      apply(env$Sdx[, env$I_,          drop = FALSE], 2L, function(sdx) env$stateSupport[sdx]),
      apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L, function(sdx) pmax(env$extendedStateSupport[sdx], 0L)),
     -apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L, function(sdx) pmin(env$extendedStateSupport[sdx], 0L))
    ) %*% c(env$alpha))

  for (t in seq(env$tau, 1L)) {
    for (i in seq(env$nSdx)) {
      hold.t <- h.t(env,
                    env$stateSupport[env$Sdx[i, env$I_]],
                    env$extendedStateSupport[env$Sdx[i, env$nI + env$J_]],
                    env$alpha)
      for (j in seq(env$nAdx)) {
        p      <- (((t - 1L) * env$nSdx + (i - 1L)) * env$nAdx + (j - 1L)) * env$nScen
        next_p <- p + seq(env$nScen)
        kdx    <- which(!is.na(transit[next_p, 1L]))
        if (length(kdx) == 0L) next

        next_i <- transit[next_p, 1L][kdx]
        costs  <- transit[next_p, 2L][kdx]
        prob   <- env$scnpb[kdx]

        Q[t, (i - 1L) * env$nAdx + j] <- -hold.t - sum(prob * costs) + sum(V[t + 1L, next_i] * prob)
      }
      idx    <- (i - 1L) * env$nAdx + env$Adx
      V[t, i] <- max(Q[t, idx], na.rm = TRUE)
      Qxs     <- Q[t, idx] - min(Q[t, idx], na.rm = TRUE) + 1.0
      Qxs[is.nan(Qxs)] <- NA_real_
      pi_rand[t, idx] <- data.table::fcoalesce(Qxs / sum(Qxs, na.rm = TRUE), 0.0)
      pi_star[t, idx] <- data.table::fcoalesce(
        as.numeric(Q[t, idx] == max(Q[t, idx], na.rm = TRUE)), 0.0)
    }
  }
  list("V" = V, "Q" = Q, "pi_star" = pi_star, "pi_rand" = pi_rand)
}

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

            # Uncertainty key encoding
            # kdx <- ((k3 - 1L) * env$nQdx + (k2 - 1L)) * env$nDdx + k1 # (Only works for 1x1 instance)
            kdx <- sum(c(env$Qdx[k1,] - 1L, env$Ddx[k2,] - 1L, env$Wdx[k3,] - 1L) * env$flowKeys) + 1L

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
          pmax(pmin(env.dp$SI_[env.dp$Sdx[i,env$I_]] + q - xI, env$R), 0L),
          pmin(pmax(env.dp$SJ_[env.dp$Sdx[i,env$nI+env$J_]] - d + xJ, -env$R), env$R),
          env.dp$max.S, weights)

        # Store into the transition matrix (objval only; h.t is added by dynamic_programming)
        transit[((i-start)*env.dp$nA+(j-1L))*nScen+k,] <- c(
          next_i,
          optx$objval,
          i, j, scndx)
      }
    }
  }
  transit
}

#' Rolling DP via C++ Bellman Updates
#'
#' Memory-efficient backward induction using \code{bellmanUpdateCx} for each
#' period. Avoids materialising the full \code{tau * nSdx * nAdx * nScen x 6}
#' transit table; working memory is O(nSdx * nAdx) per period.
#'
#' @param env Environment produced by \code{dp_config} (must have \code{I_},
#'   \code{J_}, \code{Sdx}, \code{stateSupport}, \code{extendedStateSupport},
#'   \code{scnpb}, \code{alpha}, \code{nSdx}, \code{nAdx}, \code{tau}).
#' @param jsonFile Path to the instance JSON file (passed to \code{bellmanUpdateCx}).
#' @param numThreads Number of OMP threads for \code{bellmanUpdateCx}. Default 8.
#' @return A named list with elements \code{V}, \code{Q}, \code{pi_star},
#'   \code{pi_rand} — identical structure to \code{dynamic_programming}.
#' @export
rolling_dp_cx <- function(env, jsonFile, numThreads = 8L,
                          traversalOrder = "lexicographic", chunkSize = 32L) {

  V       <- matrix(NA_real_, nrow = env$tau + 1L, ncol = env$nSdx)
  Q       <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  pi_star <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  pi_rand <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)

  # Terminal value — identical to dynamic_programming
  V[env$tau + 1L, ] <- -c(
    cbind(
      apply(env$Sdx[, env$I_,          drop = FALSE], 2L,
            function(sdx) env$stateSupport[sdx]),
      apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
            function(sdx) pmax(env$extendedStateSupport[sdx], 0L)),
     -apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
            function(sdx) pmin(env$extendedStateSupport[sdx], 0L))
    ) %*% c(env$alpha))

  for (t in seq(env$tau, 1L)) {
    step <- bellmanUpdateCx(
      jsonFile     = jsonFile,
      t            = t - 1L,              # 0-based for C++
      stateSupport = as.double(env$stateSupport),
      flowSupport  = as.double(env$Q$vals),
      scnpb        = env$scnpb,
      alpha        = env$alpha,
      V_next         = V[t + 1L, ],
      numThreads     = numThreads,
      traversalOrder = traversalOrder,
      chunkSize      = chunkSize
    )
    V      [t, ] <- step$V_t
    Q      [t, ] <- step$Q_t
    pi_star[t, ] <- step$pi_star_t
    pi_rand[t, ] <- step$pi_rand_t
  }

  list("V" = V, "Q" = Q, "pi_star" = pi_star, "pi_rand" = pi_rand)
}

#' Rolling DP via XPtr Bellman Updates
#'
#' XPtr-based variant of \code{rolling_dp_cx}. Loads the instance JSON once
#' into a \code{ProblemData} external pointer and passes it to
#' \code{bellmanUpdatePtr} each period, eliminating repeated file I/O and
#' JSON parsing.
#'
#' @param env Environment produced by \code{dp_config}.
#' @param jsonFile Path to the instance JSON file.
#' @param numThreads Number of OMP threads. Default 8.
#' @return Same structure as \code{rolling_dp_cx}: a named list with
#'   \code{V}, \code{Q}, \code{pi_star}, \code{pi_rand}.
#' @export
rolling_dp_ptr <- function(env, jsonFile, numThreads = 8L,
                           traversalOrder = "lexicographic", chunkSize = 32L) {

  prob <- loadProblemDataCx(jsonFile)

  V       <- matrix(NA_real_, nrow = env$tau + 1L, ncol = env$nSdx)
  Q       <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  pi_star <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)
  pi_rand <- matrix(NA_real_, nrow = env$tau,       ncol = env$nSdx * env$nAdx)

  V[env$tau + 1L, ] <- -c(
    cbind(
      apply(env$Sdx[, env$I_,          drop = FALSE], 2L,
            function(sdx) env$stateSupport[sdx]),
      apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
            function(sdx) pmax(env$extendedStateSupport[sdx], 0L)),
     -apply(env$Sdx[, env$nI + env$J_, drop = FALSE], 2L,
            function(sdx) pmin(env$extendedStateSupport[sdx], 0L))
    ) %*% c(env$alpha))

  for (t in seq(env$tau, 1L)) {
    step <- bellmanUpdatePtr(
      problem_ptr  = prob,
      t            = t - 1L,
      stateSupport = as.double(env$stateSupport),
      flowSupport  = as.double(env$Q$vals),
      scnpb        = env$scnpb,
      alpha        = env$alpha,
      V_next         = V[t + 1L, ],
      numThreads     = numThreads,
      traversalOrder = traversalOrder,
      chunkSize      = chunkSize
    )
    V      [t, ] <- step$V_t
    Q      [t, ] <- step$Q_t
    pi_star[t, ] <- step$pi_star_t
    pi_rand[t, ] <- step$pi_rand_t
  }

  list("V" = V, "Q" = Q, "pi_star" = pi_star, "pi_rand" = pi_rand)
}
