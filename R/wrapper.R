#' Process List of Vectors
#'
#' This function creates a list of numeric vectors and passes it to a C++ function
#' to be processed.
#'
#' @examples
#' createAndPassList()
#' @useDynLib TLPR
#' @importFrom Rcpp sourceCpp
#' @export
createAndPassList <- function(list_of_vectors, is_map = FALSE, key_type = "int") {
    # Call the C++ function to process the list
    invisible(.Call('_TLPR_processListSEXP', list_of_vectors, is_map, key_type))
}

#' Process R matrix object
#'
#' This function takes an R matrix object and passes it to a C++ function to be
#' processed.
#'
#' @examples
#' passRmatrixToCPP()
#' @useDynLib TLPR
#' @importFrom Rcpp sourceCpp
#' @export
passRmatrixToCPP <- function(R_matrix) {
    # Call the C++ function to process the list
    .Call('_TLPR_r_to_eigen', R_matrix)
}

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

# #' Optimize Model from JSON file
# #'
# #' This function takes a JSON file and optimizes the model using a C++ function.
# #'
# #' @param json_file Path to the JSON file.
# #' @examples
# #' result <- optimizeModelFromJSON("path/to/json/file.json")
# #' print(result$objval)
# #' print(result$x)
# #' @useDynLib TLPR
# #' @importFrom Rcpp sourceCpp
# #' @export
# optimizeModelFromJSON <- function(jsonFile) {
#     # Call the C++ function to optimize the model
#     .Call('_TLPR_optimizeModelFromJSON', jsonFile)
# }
