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
createAndPassList <- function(list_of_vectors) {
    # Call the C++ function to process the list
    invisible(.Call('_TLPR_processListSEXP', list_of_vectors))
}
