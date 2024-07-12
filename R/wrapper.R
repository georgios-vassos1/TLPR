#' Compute Cartesian Indices
#'
#' This function takes a list of integer vectors and computes the Cartesian indices
#' using a C++ function.
#'
#' @param vectors A list of integer vectors.
#' @return A matrix of integers representing the Cartesian indices.
#' @examples
#' result <- computeCartesianIndices(c(1, 2), c(3, 4), c(5, 6))
#' print(result)
#' @useDynLib TLPR
#' @importFrom Rcpp sourceCpp
#' @export
CartesianProduct <- function(...) {
    vectors <- list(...)
    # Call the C++ function to compute Cartesian indices
    .Call('_TLPR_CartesianProductRcpp', vectors)
}

#' Compute Cartesian Indices
#'
#' This function takes a list of integer vectors and computes the Cartesian indices
#' using a C++ function.
#'
#' @param ... Number of integer vectors.
#' @param numThreads Number of threads to use for parallel computation.
#' @return A matrix of integers representing the Cartesian indices.
#' @examples
#' result <- CartesianProductX(c(1, 2), c(3, 4), c(5, 6), numThreads = 4)
#' print(result)
#' @useDynLib TLPR 
#' @importFrom Rcpp sourceCpp
#' @export
CartesianProductX <- function(..., numThreads = 8) {
  vectors <- list(...)
  # Call the C++ function to compute Cartesian indices
  .Call('_TLPR_CartesianProductRcppParallel', vectors, numThreads)
}

#' Compute Cartesian Indices
#'
#' This function takes a list of integer vectors and computes the Cartesian indices
#' using a C++ function.
#'
#' @param ... Number of integer vectors.
#' @param numThreads Number of threads to use for parallel computation.
#' @return A matrix of integers representing the Cartesian indices.
#' @examples
#' result <- CartesianProductX(c(1, 2), c(3, 4), c(5, 6), numThreads = 4)
#' print(result)
#' @useDynLib TLPR 
#' @importFrom Rcpp sourceCpp
#' @export
CartesianProductLB <- function(..., numThreads = 8) {
  vectors <- list(...)
  # Call the C++ function to compute Cartesian indices
  .Call('_TLPR_CartesianProductRcppParallelxLB', vectors, numThreads)
}
