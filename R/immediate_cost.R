#' h.t Function
#' 
#' Calculate some cost based on inputs.
#' 
#' @param env The environment variable.
#' @param SIt Some variable.
#' @param SJt Another variable.
#' @param alpha A coefficient vector. Defaults to NULL.
#' @return The calculated cost.
#' @export
h.t <- function(env, SIt, SJt, alpha=NULL) {
  if (is.null(alpha)) {
    alpha <- rep(1.0, (env$nI + 2L * env$nJ))
  }
  # Inventory holding cost (origin)
  c(SIt %*% alpha[seq(env$nI)]) + 
    # Inventory holding cost (destination)
    c(pmax(SJt, 0.0) %*% alpha[env$nI+seq(env$nJ)]) - 
    # Backorder cost (destination)
    c(pmin(SJt, 0.0) %*% alpha[env$nI+env$nJ+seq(env$nJ)])
}

#' carrier_capacity Function
#' 
#' Calculate carrier capacity constraints.
#' 
#' @param env The environment variable.
#' @return A list with constraints.
#' @import Matrix
carrier_capacity <- function(env) {
  ## Capacity constraints

  # Strategic carriers
  A1 <- spMatrix(ncol = env$nvars, nrow = env$nCS)
  for (k in env$CS) {
    A1[k, sum(env$nLc[1L:k])+seq(env$nLc[k+1L])] <- 1L
  }
  rhs1 <- NULL # env$Cb
  sns1 <- rep("<", env$nCS)

  # Spot carriers
  A2 <- spMatrix(ncol = env$nvars, nrow = env$nCO)
  for (k in seq(env$nCO)) {
    A2[k, env$nL_+(k-1L)*env$nL+seq(env$nL)] <- 1L
  }
  rhs2 <- NULL # env$Co
  sns2 <- rep("<", env$nCO)

  list('A' = rbind(A1, A2), 'b' = c(rhs1, rhs2), 's' = c(sns1, sns2))
}

#' storage_limit Function
#' 
#' Calculate storage constraints.
#' 
#' @param env The environment variable.
#' @return A list with constraints.
#' @import Matrix
storage_limit <- function(env, ...) {
  # Storage constrains
  A <- spMatrix(ncol = env$nvars, nrow = env$nJ)
  for (j in seq(env$nJ)) {
    if (is.null(to_j <- env$to_j[[j]])) {
      # All possible routes to destination j from all origins
      idx <- (j-1L)*env$nI + seq(env$nI)
      # All indices in the bids that do to destination j
      msk <- which(apply(outer(env$L_, idx, '=='), 1L, any))
      # Adding spot indices that go to destination j
      idx <- env$nL_ + c(outer(idx, (seq(env$nCO)-1L)*env$nL, '+'))
      ##
      ## Reverse indexing
      # idx <- (seq(env$nJ) - 1L)*env$nI
      # msk <- which(apply(outer(env$L_, idx + j, '=='), 1L, any))
      # idx <- env$nL_ + c(outer(idx, (seq(env$nCO) - 1L)*env$nL, '+')) + j
      ##
      # A[j, c(msk, idx)] <- 1L
      to_j <- c(msk, idx)
    }
    A[j, to_j] <- 1L
  }
  rhs <- NULL
  sns <- rep("<", env$nJ)

  list('A' = A, 'b' = rhs, 's' = sns)
}

#' positivity_ Function
#' 
#' Calculate positivity constraints.
#' 
#' @param env The environment variable.
#' @return A list with constraints.
#' @import Matrix
positivity_ <- function(env, ...) {
  # Positivity constrains
  A <- spMatrix(ncol = env$nvars, nrow = env$nI)
  for (i in seq(env$nI)) {
    if (is.null(from_i <- env$from_i[[i]])) {
      # All possible routes from origin i to the destinations
      idx <- (seq(env$nJ) - 1L) * env$nI
      # All indices in the bids that start from origin i
      msk <- which(apply(outer(env$L_, idx + i, '=='), 1L, any))
      # Adding spot indices that start from origin i
      idx <- env$nL_ + c(outer(idx, (seq(env$nCO) - 1L)*env$nL, '+')) + i
      ##
      ## Reverse indexing
      # idx <- (i - 1L)*env$nJ + seq(env$nJ)
      # msk <- which(apply(outer(env$L_, idx, '=='), 1L, any))
      # idx <- env$nL_ + c(outer(idx, (seq(env$nCO) - 1L)*env$nL, '+'))
      ##
      # A[i, c(msk, idx)] <- 1L
      from_i <- c(msk, idx)
    }
    A[i, from_i] <- 1L
  }
  rhs <- NULL
  sns <- rep("<", env$nI)

  list('A' = A, 'b' = rhs, 's' = sns)
}

#' order_quantity Function
#' 
#' Calculate order quantity constraints.
#' 
#' @param env The environment variable.
#' @return A list with constraints.
#' @import Matrix
order_quantity <- function(env, ...) {
  # Order quantity constraint
  A   <- spMatrix(ncol = env$nvars, nrow = 1L, i = rep(1L, env$nvars), j = seq(env$nvars), x = rep(1L, env$nvars))
  rhs <- NULL
  sns <- "="

  list('A' = A, 'b' = rhs, 's' = sns)
}

#' create_model Function
#' 
#' Create a model with specified constraints.
#' 
#' @param env The environment variable.
#' @param ... Additional arguments.
#' @return A model.
#' @export
create_model <- function(env, Idx = seq(4L), ...) {
  args <- list(...)
  if (constraints <- is.null(args[["constraints"]])) {
    constraints <- c(carrier_capacity, storage_limit, positivity_, order_quantity)[Idx]
  }

  constr_list <- lapply(constraints, do.call, args = list(env))
  A   <- do.call(rbind, sapply(constr_list, function(item) item$A, simplify = FALSE))
  rhs <- do.call(c, sapply(constr_list, function(item) item$b))
  sns <- do.call(c, sapply(constr_list, function(item) item$s))

  model <- list()

  model$A          <- A
  model$obj        <- NULL
  model$modelsense <- 'min'
  model$rhs        <- rhs
  model$sense      <- sns
  model$lb         <- rep(0L, env$nvars)
  model$vtype      <- rep('I', env$nvars)

  model
}

#' Random Single Stage Shipment Assignment
#' 
#' Generates a random assignment of shipments.
#' 
#' @param n The number of shipments.
#' @param k The number of carriers.
#' @param obj_ The objective function coefficients.
#' @param ... Additional arguments.
#' @return A list containing the assignment and objective value.
#' @export
random_assignment <- function(n, k, obj_, ...) {
  if (n == 0L) {
    return(list(
      "x" = numeric(k),
      "objval" = 0.0
    ))
  }
  # Generate n-1 random integers between 1 and n-1
  allocation <- sort(sample(n-1L, k-1L, replace = TRUE))
  # Calculate the differences between consecutive elements
  x <- c(allocation[1L], diff(allocation), n - allocation[k-1L])
  # Return structure
  list(
    "x" = x,
    "objval" = c(obj_ %*% x)
  )
}

#' Capacitated Random Assignment
#' 
#' Randomly assigns items to bins while respecting capacity constraints.
#' 
#' @param n Integer, number of items.
#' @param k Integer, number of bins.
#' @param obj_ Numeric vector, objective coefficients.
#' @param x Numeric vector, bin capacities.
#' @return List with 'x' (allocation) and 'objval' (objective value).
#' @examples
#' # Assign 10 items to 3 bins with capacities 5, 5, and 10 respectively
#' cap_assignment <- capacitated_random_assignment(10, 3, c(1, 2, 3), c(5, 5, 10))
#' cap_assignment$x
#' cap_assignment$objval
#' @export
capacitated_random_assignment <- function(n, k, obj_, x, ...) {
  if (n == 0L) {
    return(list(
      "x" = numeric(k),
      "objval" = 0.0
    ))
  }

  allocation <- numeric(k)
  residue <- min(n, sum(x))

  while (residue > 0L) {
    for (i in seq(k)) {
      a_i <- sample(min(x[i], residue), 1L)
      allocation[i] <- allocation[i] + a_i
      x[i] <- max(x[i] - a_i, 0L)
      residue <- residue - a_i
    }
  }

  # Return structure
  list(
    "x" = allocation,
    "objval" = c(obj_ %*% allocation)
  )
}

#' Optimal Single Stage Shipment Assignment
#' 
#' Solves for the optimal assignment of shipments.
#' 
#' @param model The optimization model.
#' @param obj_ The objective function coefficients.
#' @param rhs_ The right-hand side of constraints.
#' @param params A list of parameters for the optimization solver. Defaults to list(OutputFlag = 0L).
#' @param ... Additional arguments.
#' @return The result of the optimization.
#' @export
#' @import gurobi
optimal_assignment <- function(model, obj_, rhs_, params = list(OutputFlag = 0L), ...) {
  # Update dynamic parameters
  model$obj <- obj_
  model$rhs <- c(model$rhs, rhs_)
  # solve
  gurobi(model, params = params)
}

#' Example Demo
#' 
#' Demonstrates an example scenario.
#' 
#' @param env The environment variable.
#' @param ... Additional arguments.
#' @return A list containing the status, assignment, and objective value.
#' @export
example_demo <- function(env, ...) {
  # Set the state of the system
  S.I <- c(70L, 90L)
  S.J <- c(10L, 10L)
  # Exogenous parameters
  Q  <- rpois(env$nI, 10L)
  D  <- rpois(env$nJ, 10L)
  # Set the total order quantity
  A.t <- 100L
  # Set the model
  model <- create_model(env)
  # Optimization
  obj_ <- c(env$CTb, env$CTo[1L,])
  rhs_ <- c(rep(10L, env$nCS), rep(50L, env$nCO), env$R - S.J, S.I + Q, A.t)
  result <- optimal_assignment(model, obj_, rhs_)
  list(
    "status"     = result$status,
    "assignment" = result$x,
    "objective"  = result$objval
  )
}
