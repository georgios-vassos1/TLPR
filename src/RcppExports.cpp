// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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
// updateStateIdx
std::vector<int> updateStateIdx(const std::vector<int>& stateIdx, const std::vector<int>& inflowIdx, const std::vector<std::vector<int>>& outflowIndices, const std::vector<double>& stateSupport, const std::vector<double>& extendedStateSupport, const std::vector<double>& flowSupport, const std::vector<double>& xI, const std::vector<double>& xJ, const double& storageLimit, const std::vector<int>& stateKeys, const int& nOrigins, const int& nDestinations);
RcppExport SEXP _TLPR_updateStateIdx(SEXP stateIdxSEXP, SEXP inflowIdxSEXP, SEXP outflowIndicesSEXP, SEXP stateSupportSEXP, SEXP extendedStateSupportSEXP, SEXP flowSupportSEXP, SEXP xISEXP, SEXP xJSEXP, SEXP storageLimitSEXP, SEXP stateKeysSEXP, SEXP nOriginsSEXP, SEXP nDestinationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type stateIdx(stateIdxSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type inflowIdx(inflowIdxSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<int>>& >::type outflowIndices(outflowIndicesSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type stateSupport(stateSupportSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type extendedStateSupport(extendedStateSupportSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type flowSupport(flowSupportSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type xI(xISEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type xJ(xJSEXP);
    Rcpp::traits::input_parameter< const double& >::type storageLimit(storageLimitSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type stateKeys(stateKeysSEXP);
    Rcpp::traits::input_parameter< const int& >::type nOrigins(nOriginsSEXP);
    Rcpp::traits::input_parameter< const int& >::type nDestinations(nDestinationsSEXP);
    rcpp_result_gen = Rcpp::wrap(updateStateIdx(stateIdx, inflowIdx, outflowIndices, stateSupport, extendedStateSupport, flowSupport, xI, xJ, storageLimit, stateKeys, nOrigins, nDestinations));
    return rcpp_result_gen;
END_RCPP
}
// computeEnvironmentCx
Eigen::MatrixXd computeEnvironmentCx(const std::string jsonFile, const std::vector<double>& stateSupport, const std::vector<double>& flowSupport, int numThreads);
RcppExport SEXP _TLPR_computeEnvironmentCx(SEXP jsonFileSEXP, SEXP stateSupportSEXP, SEXP flowSupportSEXP, SEXP numThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type jsonFile(jsonFileSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type stateSupport(stateSupportSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type flowSupport(flowSupportSEXP);
    Rcpp::traits::input_parameter< int >::type numThreads(numThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(computeEnvironmentCx(jsonFile, stateSupport, flowSupport, numThreads));
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
Eigen::MatrixXd rmvnorm(int n, const Eigen::VectorXd& mean, const Eigen::MatrixXd& covar, int nThreads);
RcppExport SEXP _TLPR_rmvnorm(SEXP nSEXP, SEXP meanSEXP, SEXP covarSEXP, SEXP nThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type covar(covarSEXP);
    Rcpp::traits::input_parameter< int >::type nThreads(nThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvnorm(n, mean, covar, nThreads));
    return rcpp_result_gen;
END_RCPP
}
// convertListToMapTest
void convertListToMapTest(SEXP rlist, const std::string& KeyTypeArg);
RcppExport SEXP _TLPR_convertListToMapTest(SEXP rlistSEXP, SEXP KeyTypeArgSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rlist(rlistSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type KeyTypeArg(KeyTypeArgSEXP);
    convertListToMapTest(rlist, KeyTypeArg);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TLPR_CartesianProductRcpp", (DL_FUNC) &_TLPR_CartesianProductRcpp, 1},
    {"_TLPR_CartesianProductRcppParallel", (DL_FUNC) &_TLPR_CartesianProductRcppParallel, 2},
    {"_TLPR_CartesianProductRcppParallelxLB", (DL_FUNC) &_TLPR_CartesianProductRcppParallelxLB, 2},
    {"_TLPR_updateStateIdx", (DL_FUNC) &_TLPR_updateStateIdx, 12},
    {"_TLPR_computeEnvironmentCx", (DL_FUNC) &_TLPR_computeEnvironmentCx, 4},
    {"_TLPR_optimizeModelFromJSON", (DL_FUNC) &_TLPR_optimizeModelFromJSON, 1},
    {"_TLPR_rmvnorm", (DL_FUNC) &_TLPR_rmvnorm, 4},
    {"_TLPR_convertListToMapTest", (DL_FUNC) &_TLPR_convertListToMapTest, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_TLPR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
