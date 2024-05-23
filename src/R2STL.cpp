#include <Rcpp.h>
#include <vector>
#include <iostream>

using namespace Rcpp;

// Function to convert R list of vectors to std::vector<std::vector<int>>
std::vector< std::vector<int> > convertListToVector(SEXP list_sexp) {
    List list_of_vectors(list_sexp);
    std::vector< std::vector<int> > vec;

    for (int i = 0; i < list_of_vectors.size(); i++) {
        NumericVector temp = list_of_vectors[i];
        std::vector<int> temp_vec(temp.begin(), temp.end());
        vec.push_back(temp_vec);
    }

    return vec;
}

// Function to print the std::vector<std::vector<int>>
void printVector(const std::vector< std::vector<int> >& vec) {
    for (const auto& sub_vec : vec) {
        for (const auto& val : sub_vec) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

// Wrapper function to be called from R
// [[Rcpp::export]]
void processListSEXP(SEXP list_sexp) {
    std::vector< std::vector<int> > vec = convertListToVector(list_sexp);
    printVector(vec);
}
