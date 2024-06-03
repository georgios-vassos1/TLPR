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
