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

#' Plot 2x2 Instance
#'
#' This function plots a series of graphs representing the state of the environment at different times.
#'
#' @param env A list containing the environment variables including `tau`, `nI`, and `nJ`.
#' @param graph.dt A data.table containing the graph data with columns `t`, `origin`, `destination`, and `assignment`.
#' @param S.I A vector representing the source nodes.
#' @param S.J A vector representing the destination nodes.
#' @param Q A vector representing some quantities related to the nodes.
#' @param D A vector representing some other quantities related to the nodes.
#' 
#' @import data.table igraph
#' @examples
#' # Assuming `env`, `graph.dt`, `S.I`, `S.J`, `Q`, and `D` are defined
#' plot2x2instance(env, graph.dt, S.I, S.J, Q, D)
#'
#' @export
plot2x2instance <- function(env, graph.dt, S.I, S.J, Q, D) {
  n <- 2L
  reflex <- diag(n)[, seq(n, 1L)]
  
  # Set up the plotting area to display multiple plots
  par(mfrow = c(3L, 4L))
  for (t_ in seq(env$tau)) {
    # Create an igraph object for the current time step
    g <- graph_from_data_frame(graph.dt[t == t_, -1L], directed = TRUE, vertices = seq(env$nI + env$nJ))
    
    # Set vertex attributes
    V(g)$color <- "lightblue"
    V(g)$size <- 50L  # Size of the vertices
    
    # Set edge attributes
    E(g)$arrow.size <- 0.25  # Size of the arrows
    E(g)$weight <- graph.dt$assignment[(t_ - 1L) * env$nL + seq(env$nL)]
    
    # Plot the graph with edge labels
    plot(g, edge.label = E(g)$weight, layout = igraph::layout_on_grid(g) %*% reflex, main = paste0("t = ", t_))
    
    # Add text labels to the plot
    text(x = c(-1.34, -1.34, 1.34, 1.34, -1.00, -1.00, 1.00, 1.00), 
         y = c(-1.00, 1.00, -1.00, 1.00, -0.60, 0.60, -0.60, 0.60), 
         labels = c(
           Q[(t_ - 1L) * env$nI + seq(env$nI)],
           D[(t_ - 1L) * env$nI + seq(env$nI)],
           S.I[(t_ - 1L) * env$nI + env$I_, 1L],
           S.J[(t_ - 1L) * env$nJ + env$J_, 1L]
         ),
         adj = c(.5, .5))
  }
  # Reset plotting area to single plot
  par(mfrow = c(1L, 1L))
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
