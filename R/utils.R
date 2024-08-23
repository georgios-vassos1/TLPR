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
  X[, seq(ncol(X), 1L)]
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

## Parallel run
#' Chunk Up a Sequence
#'
#' This function divides a sequence of length `n` into `k` approximately equal-sized chunks.
#'
#' @param n Integer. The total length of the sequence to be divided.
#' @param k Integer. The number of chunks to divide the sequence into.
#'
#' @return A vector of integers representing the start of each chunk. The last element in the vector will always be `n + 1`, which represents the end of the last chunk.
#'
#' @examples
#' # Divide a sequence of 10 elements into 3 chunks
#' chunkup(10, 3)
#' # This will return something like: c(1, 4, 7, 11)
#'
#' # Divide a sequence of 15 elements into 4 chunks
#' chunkup(15, 4)
#' # This will return something like: c(1, 4, 8, 12, 16)
#'
#' @export
chunkup <- function(n, k) {
  chunk_size <- n %/% k
  chunks <- seq(1L, n, by = chunk_size)
  if (length(chunks) < k + 1) {
    chunks <- c(chunks, n + 1L)
  } else {
    chunks[k + 1] <- n + 1L
  }
  chunks
}

#' Compute Graph Data Table
#'
#' This function computes a data.table representing the graph for a given environment,
#' source nodes (S.I), destination nodes (S.J), and allocation matrix.
#'
#' @param env A list containing the environment variables including `L_`, `nL`, `tau`, `nvars`, and `L`.
#' @param S.I A vector representing the source nodes.
#' @param S.J A vector representing the destination nodes.
#' @param allocation A matrix representing the allocation data.
#' 
#' @return A data.table containing the columns: `t` (time), `origin`, `destination`, and `assignment`.
#' The `assignment` column is the sum of allocations for each combination of time, origin, and destination.
#' 
#' @import data.table
#' @examples
#' # Assuming `env`, `S.I`, `S.J`, and `allocation` are defined
#' graph_dt <- compute_graph_dt(env, S.I, S.J, allocation)
#'
#' @export
compute_graph_dt <- function(env, S.I, S.J, allocation) {
  # Create a vector of lanes
  lanes <- c(env$L_, seq(env$nL))
  
  # Construct a data table with columns t, lane, origin, destination, and assignment
  dt <- as.data.table(cbind(
    t           = rep(seq(env$tau), each = env$nvars),
    lane        = lanes,
    origin      = env$L[lanes, 1L],
    destination = env$L[lanes, 2L] + env$nI,
    assignment  = allocation[, 1L]
  ))
  
  # Aggregate the assignments by time, origin, and destination
  dt[, .(assignment = sum(assignment, na.rm = TRUE)), by = .(t, origin, destination)]
}

#' Generate a Positive Definite Matrix
#'
#' This function generates a random symmetric positive definite matrix of a given size.
#'
#' @param p An integer specifying the dimension of the matrix (number of rows and columns).
#' @return A p x p symmetric positive definite matrix.
#' @examples
#' # Generate a 3x3 positive definite matrix
#' matrix <- generatePositiveDefiniteMatrix(3)
#' print(matrix)
#' @export
generatePositiveDefiniteMatrix <- function(p) {
  # Generate a random p x p matrix
  A <- matrix(rnorm(p * p), nrow = p, ncol = p)
  
  # Create a symmetric positive definite matrix by multiplying A by its transpose
  symmetricMatrix <- t(A) %*% A
  
  # Adding a small value to the diagonal to ensure positive definiteness
  symmetricMatrix <- symmetricMatrix + diag(p) * 1e-3
  
  return(symmetricMatrix)
}

#' Create a Gradient Color Palette
#'
#' This function generates a gradient color palette between two colors, specified by `low` and `high`.
#' If `norm_` is provided, it maps the normalized z-values to corresponding colors in the gradient palette.
#'
#' @param n Integer. The number of colors to generate in the gradient palette.
#' @param norm_ Numeric vector (optional). A vector of normalized values to map to the gradient palette.
#' @param low Character. The hex code for the color representing the low end of the gradient. Default is `#80FFFF`.
#' @param high Character. The hex code for the color representing the high end of the gradient. Default is `#FF80FF`.
#'
#' @return A vector of colors. If `norm_` is not provided, it returns a vector of `n` colors representing the gradient. 
#'         If `norm_` is provided, it returns a vector of colors corresponding to the mapped values.
#'
#' @examples
#' create_gradient_palette(5)
#' create_gradient_palette(10, norm_ = runif(10))
#'
#' @export
create_gradient_palette <- function(n, norm_=NULL, low = "#80FFFF", high = "#FF80FF") {
  # Create a gradient color palette
  gradient_colors <- colorRampPalette(c(low, high))(n)
  if (is.null(norm_)) return(gradient_colors)

  # Map normalized z-values to colors in the gradient palette
  gradient_colors[cut(norm_, breaks = n, include.lowest = TRUE)]
}

#' Action-Value Surface for Fixed Exit State
#'
#' This function computes the action-value surface given a fixed exit state in a specified environment.
#' It calculates the state-action index values based on the state support and the environment configuration.
#'
#' @param env List. The environment object containing the state support, action indices, and other relevant parameters.
#' @param exit_state Numeric vector. The fixed exit state used to calculate the action-value surface.
#' @param ... Additional arguments passed to the function.
#'
#' @return An integer vector representing the action-value surface indices.
#'
#' @examples
#' env <- list(nI = 2, stateSupport = c(1, 2), R = 1, nAdx = 3, stateKeys = c(1, 10))
#' exit_state <- c(1, 2)
#' get_Qdx_fixed_exit(env, exit_state)
#'
#' @export
get_Qdx_fixed_exit <- function(env, exit_state, ...) {
  Sdx  <- do.call(TLPR::CartesianProductX, c(replicate(env$nI, seq_along(env$stateSupport), simplify = FALSE)))
  nSdx <- nrow(Sdx)

  rng  <- seq(env$nAdx)

  Qdx <- integer(nSdx * env$nAdx)
  for (i in seq(nSdx)) {
    idx <- (sum(c(env$stateSupport[Sdx[i,]], env$R + exit_state) * env$stateKeys) + 1L)
    Qdx[(i - 1L) * env$nAdx + rng] <- (idx - 1L) * env$nAdx + rng
  }

  Qdx
}

#' Action-Value Surface for Fixed Entry State
#'
#' This function computes the action-value surface given a fixed entry state in a specified environment.
#' It calculates the state-action index values based on the extended state support and the environment configuration.
#'
#' @param env List. The environment object containing the extended state support, action indices, and other relevant parameters.
#' @param entry_state Numeric vector. The fixed entry state used to calculate the action-value surface.
#' @param ... Additional arguments passed to the function.
#'
#' @return An integer vector representing the action-value surface indices.
#'
#' @examples
#' env <- list(nJ = 2, extendedStateSupport = c(1, 2), R = 1, nAdx = 3, stateKeys = c(1, 10))
#' entry_state <- c(1, 2)
#' get_Qdx_fixed_entry(env, entry_state)
#'
#' @export
get_Qdx_fixed_entry <- function(env, entry_state, ...) {
  Sdx  <- do.call(TLPR::CartesianProductX, c(replicate(env$nJ, seq_along(env$extendedStateSupport), simplify = FALSE)))
  nSdx <- nrow(Sdx)

  rng  <- seq(env$nAdx)

  Qdx <- integer(nSdx * env$nAdx)
  for (j in seq(nSdx)) {
    jdx <- (sum(c(entry_state, env$R + env$extendedStateSupport[Sdx[j,]]) * env$stateKeys) + 1L)
    Qdx[(j - 1L) * env$nAdx + rng] <- (jdx - 1L) * env$nAdx + rng
  }

  Qdx
}
