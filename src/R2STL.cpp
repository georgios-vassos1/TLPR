#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Rcpp.h>

using namespace Rcpp;

// Function to convert R list of vectors to std::map<int, std::vector<int>>
template <typename KeyType, typename ValueType>
std::map< KeyType, std::vector<ValueType> > convertListToMap(SEXP inputList) {
    std::map< KeyType, std::vector<ValueType> > result;

    // Ensure the input is a list
    if (!Rf_isNewList(inputList)) {
        Rcpp::stop("Input must be a list");
    }

    Rcpp::List list(inputList);
    bool useIndices = Rf_isNull(list.names());
    // std::cout << "useIndices: " << useIndices << std::endl;
    Rcpp::CharacterVector names;
    if (!useIndices) {
        names = list.names();
    }

    for (int i = 0; i < list.size(); i++) {
        KeyType key;
        if (std::is_same<KeyType, int>::value) {
            key = i;
        } else {
            key = Rcpp::as<KeyType>(names[i]);
        }

        SEXP valuesSEXP = list[i];
        Rcpp::NumericVector values(valuesSEXP);
        std::vector<ValueType> vec(values.begin(), values.end());

        result[key] = vec;
    }

    return result;
}

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

// This function reads a matrix from R and converts it to an Eigen matrix
// [[Rcpp::export]]
Eigen::MatrixXd r_to_eigen(const Eigen::MatrixXd& mat) {
    // Here, mat is already an Eigen::MatrixXd object because RcppEigen handles the conversion
    // You can now perform any Eigen operations on mat
    Eigen::MatrixXd result = mat.transpose(); // Example operation: transpose
    return result;
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

// Function to print the std::map<int, std::vector<int>>
template <typename KeyType, typename ValueType>
void printMap(const std::map< KeyType, std::vector<ValueType> >& map) {
    for (const auto& pair : map) {
        std::cout << pair.first << ": ";
        for (const auto& val : pair.second) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

// Wrapper function to be called from R
// [[Rcpp::export]]
void processListSEXP(SEXP list_sexp, bool is_map, const std::string KeyTypeArg) {
    if (is_map) {
        if (KeyTypeArg == "int") {
            std::map< int, std::vector<float> > map = convertListToMap<int, float>(list_sexp);
            printMap<int, float>(map);
        } else if (KeyTypeArg == "std::string") {
            std::map< std::string, std::vector<float> > map = convertListToMap<std::string, float>(list_sexp);
            printMap<std::string, float>(map);
        } else {
            Rcpp::stop("KeyTypeArg must be either 'int' or 'std::string'");
        }
    } else {
        std::vector< std::vector<int> > vec = convertListToVector(list_sexp);
        printVector(vec);
    }
}
