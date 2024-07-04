# #' Initialize the package
# #'
# #' This function initializes the package by loading required Rcpp modules. It is 
# #' automatically called when the package is loaded, ensuring all necessary components 
# #' are set up properly.
# #' @useDynLib TLPR
# #' @importFrom Rcpp loadModule
# .onLoad <- function(libname, pkgname) {
#   Rcpp::loadModule("echo_module", TRUE)
#   message("Package loaded and Rcpp module loaded.")
# }
