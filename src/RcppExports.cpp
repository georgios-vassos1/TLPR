// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// r_to_eigen
Eigen::MatrixXd r_to_eigen(const Eigen::MatrixXd& mat);
RcppExport SEXP _TLPR_r_to_eigen(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(r_to_eigen(mat));
    return rcpp_result_gen;
END_RCPP
}
// processListSEXP
void processListSEXP(SEXP list_sexp, bool is_map, const std::string KeyTypeArg);
RcppExport SEXP _TLPR_processListSEXP(SEXP list_sexpSEXP, SEXP is_mapSEXP, SEXP KeyTypeArgSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type list_sexp(list_sexpSEXP);
    Rcpp::traits::input_parameter< bool >::type is_map(is_mapSEXP);
    Rcpp::traits::input_parameter< const std::string >::type KeyTypeArg(KeyTypeArgSEXP);
    processListSEXP(list_sexp, is_map, KeyTypeArg);
    return R_NilValue;
END_RCPP
}
// CartesianProductRcpp
Eigen::MatrixXi CartesianProductRcpp(List vectors);
RcppExport SEXP _TLPR_CartesianProductRcpp(SEXP vectorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type vectors(vectorsSEXP);
    rcpp_result_gen = Rcpp::wrap(CartesianProductRcpp(vectors));
    return rcpp_result_gen;
END_RCPP
}
// CartesianProductRcppParallel
Eigen::MatrixXi CartesianProductRcppParallel(List vectors, int numThreads);
RcppExport SEXP _TLPR_CartesianProductRcppParallel(SEXP vectorsSEXP, SEXP numThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type vectors(vectorsSEXP);
    Rcpp::traits::input_parameter< int >::type numThreads(numThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(CartesianProductRcppParallel(vectors, numThreads));
    return rcpp_result_gen;
END_RCPP
}
// CartesianProductRcppParallelxLB
Eigen::MatrixXi CartesianProductRcppParallelxLB(List vectors, int numThreads);
RcppExport SEXP _TLPR_CartesianProductRcppParallelxLB(SEXP vectorsSEXP, SEXP numThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type vectors(vectorsSEXP);
    Rcpp::traits::input_parameter< int >::type numThreads(numThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(CartesianProductRcppParallelxLB(vectors, numThreads));
    return rcpp_result_gen;
END_RCPP
}
// computeEnvironmentSTL
std::vector<std::vector<int>> computeEnvironmentSTL(const std::string jsonFile, const std::vector<double>& stateSupport, int numThreads);
RcppExport SEXP _TLPR_computeEnvironmentSTL(SEXP jsonFileSEXP, SEXP stateSupportSEXP, SEXP numThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type jsonFile(jsonFileSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type stateSupport(stateSupportSEXP);
    Rcpp::traits::input_parameter< int >::type numThreads(numThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(computeEnvironmentSTL(jsonFile, stateSupport, numThreads));
    return rcpp_result_gen;
END_RCPP
}
// optimizeModelFromJSON
Rcpp::List optimizeModelFromJSON(std::string jsonFile);
RcppExport SEXP _TLPR_optimizeModelFromJSON(SEXP jsonFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type jsonFile(jsonFileSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizeModelFromJSON(jsonFile));
    return rcpp_result_gen;
END_RCPP
}
// rmvnorm
Eigen::MatrixXd rmvnorm(int n, int p, const Eigen::VectorXd& mean, const Eigen::MatrixXd& covar, int nThreads);
RcppExport SEXP _TLPR_rmvnorm(SEXP nSEXP, SEXP pSEXP, SEXP meanSEXP, SEXP covarSEXP, SEXP nThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type covar(covarSEXP);
    Rcpp::traits::input_parameter< int >::type nThreads(nThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvnorm(n, p, mean, covar, nThreads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TLPR_r_to_eigen", (DL_FUNC) &_TLPR_r_to_eigen, 1},
    {"_TLPR_processListSEXP", (DL_FUNC) &_TLPR_processListSEXP, 3},
    {"_TLPR_CartesianProductRcpp", (DL_FUNC) &_TLPR_CartesianProductRcpp, 1},
    {"_TLPR_CartesianProductRcppParallel", (DL_FUNC) &_TLPR_CartesianProductRcppParallel, 2},
    {"_TLPR_CartesianProductRcppParallelxLB", (DL_FUNC) &_TLPR_CartesianProductRcppParallelxLB, 2},
    {"_TLPR_computeEnvironmentSTL", (DL_FUNC) &_TLPR_computeEnvironmentSTL, 3},
    {"_TLPR_optimizeModelFromJSON", (DL_FUNC) &_TLPR_optimizeModelFromJSON, 1},
    {"_TLPR_rmvnorm", (DL_FUNC) &_TLPR_rmvnorm, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_TLPR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
