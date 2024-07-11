#include <Eigen/Dense>
#include <random>
#include <omp.h>
#include <stdexcept>

// Function to generate samples from a multivariate normal distribution
Eigen::VectorXd multivariateNormalSample(const Eigen::VectorXd& mean, const Eigen::MatrixXd& covar, std::mt19937& gen) {

  std::normal_distribution<> dist(0.0, 1.0);

  Eigen::VectorXd randomVec(mean.size());
  for(int i = 0; i < mean.size(); ++i) {
    randomVec(i) = dist(gen);
  }

  Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);
  if (cholSolver.info() != Eigen::Success) {
    throw std::runtime_error("Decomposition failed");
  }
  Eigen::MatrixXd L = cholSolver.matrixL();

  return mean + L * randomVec;
}

// Function to fill the matrix
//' @useDynLib TLPR
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd rmvnorm(int n, const Eigen::VectorXd& mean, const Eigen::MatrixXd& covar, int nThreads = 8) {

  int p = mean.size();

  Eigen::MatrixXd outputMatrix(n, p);

  // Set the number of threads
  omp_set_num_threads(nThreads);

  // Pre-generate random number generators for each thread
  std::vector<std::mt19937> generators(nThreads);
  for (int t = 0; t < nThreads; ++t) {
    std::random_device rd;
    generators[t] = std::mt19937(rd());
  }

  // Parallelize the outer loops
  #pragma omp parallel for
  for (int idx = 0; idx < n; ++idx) {
    int thread_id = omp_get_thread_num();
    std::mt19937& gen = generators[thread_id];

    Eigen::VectorXd sample = multivariateNormalSample(mean, covar, gen);

    // Fill the output matrix
    outputMatrix.row(idx) = sample.transpose();
  }

  return outputMatrix;
}

// Illustration Only 
//
// Function to fill the matrix
Eigen::MatrixXd fillMatrixWithSamples(int nSdx, int nAdx, int nUdx, int p, const Eigen::VectorXd& mean, const Eigen::MatrixXd& covar) {

  if (p != mean.size()) {
    throw std::invalid_argument("The parameter p must be equal to the size of the mean vector.");
  }

  int n = nSdx * nAdx * nUdx;
  Eigen::MatrixXd outputMatrix(n, p);

  std::random_device rd;
  std::mt19937 gen(rd());

  for (int i = 0; i < nSdx; ++i) {
    for (int j = 0; j < nAdx; ++j) {
      for (int k = 0; k < nUdx; ++k) {

        Eigen::VectorXd sample = multivariateNormalSample(mean, covar, gen);

        // Calculate the row index
        int idx = i * nAdx * nUdx + j * nUdx + k;

        // Fill the output matrix
        outputMatrix.row(idx) = sample.transpose();
      }
    }
  }

  return outputMatrix;
}

// Function to fill the matrix
Eigen::MatrixXd fillMatrixWithSamplesOMP3(int nSdx, int nAdx, int nUdx, int p, const Eigen::VectorXd& mean, const Eigen::MatrixXd& covar, int nThreads) {

  if (p != mean.size()) {
    throw std::invalid_argument("The parameter p must be equal to the size of the mean vector.");
  }

  int n = nSdx * nAdx * nUdx;
  Eigen::MatrixXd outputMatrix(n, p);

  // Set the number of threads
  omp_set_num_threads(nThreads);

  // Parallelize the outer loops
  #pragma omp parallel for collapse(3)
  for (int i = 0; i < nSdx; ++i) {
    for (int j = 0; j < nAdx; ++j) {
      for (int k = 0; k < nUdx; ++k) {
        // Each thread will have its own random number generator
        std::random_device rd;
        static thread_local std::mt19937 gen(rd());  // thread_local for thread safety

        Eigen::VectorXd sample = multivariateNormalSample(mean, covar, gen);

        // Calculate the row index
        int idx = i * nAdx * nUdx + j * nUdx + k;

        // Directly write to the output matrix without a critical section
        outputMatrix.row(idx) = sample.transpose();
      }
    }
  }

  return outputMatrix;
}
