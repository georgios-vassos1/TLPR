#include <Eigen/Dense>
#include <random>
#include <omp.h>
#include <stdexcept>
#include <Rcpp.h>

// Function to generate samples from a multivariate normal distribution.
// Accepts the pre-computed lower Cholesky factor L instead of the covariance
// matrix to avoid redundant decompositions when called in a loop.
Eigen::VectorXd multivariateNormalSample(const Eigen::VectorXd& mean, const Eigen::MatrixXd& L,
                                         std::mt19937& gen) {

  std::normal_distribution<> dist(0.0, 1.0);

  Eigen::VectorXd randomVec(mean.size());
  for (int i = 0; i < mean.size(); ++i) {
    randomVec(i) = dist(gen);
  }

  return mean + L * randomVec;
}

// Function to fill the matrix
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd rmvnorm(int n, const Eigen::VectorXd& mean, const Eigen::MatrixXd& covar,
                        int nThreads = 8) {

  int p = mean.size();

  Eigen::MatrixXd outputMatrix(n, p);

  // Pre-compute the Cholesky factor once; reused for every sample.
  Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);
  if (cholSolver.info() != Eigen::Success)
    Rcpp::stop("Cholesky decomposition failed");
  const Eigen::MatrixXd L = cholSolver.matrixL();

  // Pre-generate random number generators for each thread
  int maxT = omp_get_max_threads();
  std::vector<std::mt19937> generators(maxT);
  for (int t = 0; t < maxT; ++t) {
    std::random_device rd;
    generators[t] = std::mt19937(rd());
  }

// Parallelize the outer loops
#pragma omp parallel for num_threads(nThreads)
  for (int idx = 0; idx < n; ++idx) {
    int thread_id = omp_get_thread_num();
    std::mt19937& gen = generators[thread_id];

    Eigen::VectorXd sample = multivariateNormalSample(mean, L, gen);

    // Fill the output matrix
    outputMatrix.row(idx) = sample.transpose();
  }

  return outputMatrix;
}
