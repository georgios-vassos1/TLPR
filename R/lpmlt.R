#' Carrier Capacity with Padding
#'
#' This function calculates the carrier capacity with added padding for the given environment.
#'
#' @param env The environment list containing various parameters such as nI, nJ, nCS, nCO, and others.
#'
#' @return A list containing the matrix 'A' with added padding, right-hand side 'rhs', and 'sense'.
#' @export
carrier_capacity_padded <- function(env) {
  # Compute standard carrier capacity
  ccstd <- carrier_capacity(env)

  # Create padding matrices
  padding <- Matrix::spMatrix(ncol = env$nI + 2L * env$nJ, nrow = env$nCS + env$nCO)

  # Return the padded carrier capacity components
  list(
    'A'     = cbind(padding, ccstd$A, padding),
    'rhs'   = c(env$Cb[1L, ], env$Co[1L, ]),
    'sense' = rep("<", env$nCS + env$nCO)
  )
}

#' Transition Logic
#'
#' This function constructs the transition logic matrix and its components for the given environment.
#'
#' @param env The environment list containing various parameters such as nI, nJ, nvars, and others.
#' @param q A numeric vector representing the initial conditions.
#' @param d A numeric vector representing the demands.
#'
#' @return A list containing the matrix 'A', right-hand side 'rhs', and 'sense'.
#' @export
transition_logic <- function(env, q, d) {
  # Initialize the transition logic matrix
  A3 <- Matrix::spMatrix(ncol = env$nI + 2L * env$nJ + env$nvars + env$nI + 2L * env$nJ, nrow = env$nI + env$nJ)

  # Fill the matrix for items
  for (i in seq(env$nI)) {
    A3[i, i] <- -1L
    idx <- env$from_i[[i]]
    A3[i, env$nI + 2L * env$nJ + idx] <- 1L
    A3[i, env$nI + 2L * env$nJ + env$nvars + i] <- 1L
  }

  # Fill the matrix for jobs
  for (j in seq(env$nJ)) {
    A3[env$nI + j, env$nI + (j - 1L) * 2L + seq(2L)] <- c(1L, -1L)
    jdx <- env$to_j[[j]]
    A3[env$nI + j, env$nI + 2L * env$nJ + jdx] <- 1L
    A3[env$nI + j, env$nI + 2L * env$nJ + env$nvars + env$nI + (j - 1L) * 2L + seq(2L)] <- c(-1L, 1L)
  }

  # Return the transition logic components
  list(
    'A'     = A3,
    'rhs'   = c(q, d),
    'sense' = rep("=", env$nI + env$nJ)
  )
}

#' Storage Limits
#'
#' This function constructs the storage limit constraints for the given environment.
#'
#' @param env The environment list containing various parameters such as nI, nJ, nvars, and others.
#' @param q A numeric vector representing the initial conditions.
#'
#' @return A list containing the matrix 'A', right-hand side 'rhs', and 'sense'.
#' @export
storage_limits <- function(env, q) {
  # Initialize the storage limit matrix for items
  A4 <- Matrix::spMatrix(ncol = env$nI, nrow = env$nI)
  for (i in seq(env$nI)) {
    A4[i, i] <- 1L
  }

  # Create padding for the A4 matrix
  padding <- Matrix::spMatrix(ncol = 2L * env$nJ + env$nvars + env$nI + 2L * env$nJ, nrow = env$nI)
  A4 <- cbind(A4, padding)

  # Initialize the storage limit matrix for jobs
  A5 <- Matrix::spMatrix(ncol = 2L * env$nJ + env$nvars, nrow = env$nJ)
  for (j in seq(env$nJ)) {
    A5[j, (j - 1L) * 2L + 1L] <- 1L
    jdx <- env$to_j[[j]]
    A5[j, 2L * env$nJ + jdx] <- 1L
  }
  
  # Create left and right padding for the A5 matrix
  lpadx <- Matrix::spMatrix(ncol = env$nI, nrow = env$nJ)
  rpadx <- Matrix::spMatrix(ncol = env$nI + 2L * env$nJ, nrow = env$nJ)
  A5 <- cbind(lpadx, A5, rpadx)
  
  # Return the storage limit components
  list(
    'A'     = rbind(A4, A5),
    'rhs'   = c(env$R - q, rep(env$R, env$nJ)),
    'sense' = rep("<", env$nI + env$nJ)
  )
}

#' Multiperiod Expansion
#'
#' This function expands the given constraints and objective function for multiple periods.
#'
#' @param env The environment list containing various parameters such as nI, nJ, nvars, tau, alpha, and others.
#' @param Q A numeric vector representing the exogenous arrivals.
#' @param D A numeric vector representing the exogenous demands.
#' @param A The constraint matrix.
#' @param obj The objective function coefficients.
#' @param rhs The right-hand side vector.
#' @param sns The sense vector.
#'
#' @return A list containing the expanded objective function 'obj', constraint matrix 'A', right-hand side 'rhs', and 'sense'.
#' @export
multiperiod_expansion <- function(env, Q, D, A, obj, rhs, sns) {
  offset <- env$nI + 2L * env$nJ
  A.tau  <- Matrix::spMatrix(nrow = env$tau * nrow(A), ncol = env$tau * (ncol(A) - offset) + offset)
  A.tau[seq(nrow(A)), seq(ncol(A))] <- A

  obj.tau <- numeric(ncol(A.tau))
  rhs.tau <- numeric(nrow(A.tau))

  obj.tau[seq(ncol(A) - offset)] <- obj[seq(ncol(A) - offset)]
  rhs.tau[seq(nrow(A))] <- rhs
  for (t in seq(2L, env$tau)) {
    rdx <- (t - 1L) * nrow(A) + seq(nrow(A))
    cdx <- (t - 1L) * (ncol(A) - offset) + seq(ncol(A))
    A.tau[rdx, cdx] <- A
    obj.tau[(t - 1L) * (ncol(A) - offset) + seq(ncol(A) - offset)] <- c(env$alpha, env$CTb, env$CTo[t, ])
    rhs.tau[(t - 1L) * nrow(A) + seq(nrow(A))] <- c(
      c(env$Cb[t, ], env$Co[t, ]),
      c(Q[(t - 1L) * env$nI + env$I_], D[(t - 1L) * env$nJ + env$J_]),
      env$R - Q[(t - 1L) * env$nI + env$I_],
      rep(env$R, env$nJ)
    )
  }
  obj.tau[seq(ncol(A.tau) - offset + 1L, ncol(A.tau))] <- env$alpha

  # Return the expanded components for multiple periods
  list(
    'obj'   = obj.tau,
    'A'     = A.tau,
    'rhs'   = rhs.tau,
    'sense' = rep(sns, env$tau)
  )
}

#' Post-Hoc Simulation
#'
#' This function simulates the results post hoc for the given environment and decision variables.
#'
#' @param env The environment list containing various parameters such as tau, nI, nJ, nvars, alpha, R, I_, J_, and others.
#' @param x A numeric vector representing the decision variables.
#'
#' @return A list containing the simulated inventory levels 'S.I', job status 'S.J', and 'allocation'.
#' @export
post_hoc_simulation <- function(env, x) {
  S.I <- matrix(NA, nrow = (env$tau + 1L) * env$nI, ncol = 1L)
  S.J <- matrix(NA, nrow = (env$tau + 1L) * env$nJ, ncol = 1L)
  allocation <- matrix(NA, env$tau * env$nvars, ncol = 1L)

  offset <- env$nI + 2L * env$nJ
  blk <- offset + env$nvars
  for (t in seq(env$tau + 1L)) {
    i <- (t - 1L) * env$nI + env$I_
    j <- (t - 1L) * env$nJ + env$J_
    idx <- (t - 1L) * blk + seq(blk)
    S.I[i, 1L] <- x[idx][env$I_]
    S.J[j, 1L] <- c(Reduce('-', x[idx][env$nI + env$J_]), Reduce('-', x[idx][env$nI + env$nJ + env$J_]))
    
    if (t > env$tau) break
    allocation[(t - 1L) * env$nvars + seq(env$nvars), 1L] <- x[idx][offset + seq(env$nvars)]
  }

  # Return the post hoc simulation results
  list(
    'S.I' = S.I,
    'S.J' = S.J,
    'allocation' = allocation
  )
}
