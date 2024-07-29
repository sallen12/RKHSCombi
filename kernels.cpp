// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;

// Euclidean norm
// [[Rcpp::export]]
double euclnormC(arma::colvec x){
  double out = sqrt(sum(square(x)));
  return(out);
}

// energy kernel
// [[Rcpp::export]]
double enerk(double x, double y){
  double out = abs(x) + abs(y) - abs(x - y);
  return(out);
}

// mv energy kernel
// [[Rcpp::export]]
double enerk_mv(arma::colvec x, arma::colvec y){
  double out = euclnormC(x) + euclnormC(y) - euclnormC(x - y);
  return(out);
}

// gaussian kernel
// [[Rcpp::export]]
double gaussk(double x, double y){
  double out = exp(-pow(abs(x - y), 2));
  return (out);
}

// mv gaussian kernel
// [[Rcpp::export]]
double gaussk_mv(arma::colvec x, arma::colvec y){
  double out = exp(-pow(euclnormC(x - y), 2));
  return (out);
}

// laplace kernel
// [[Rcpp::export]]
double laplk(double x, double y){
  double out = exp(-abs(x - y));
  return (out);
}

// mv laplace kernel
// [[Rcpp::export]]
double laplk_mv(arma::colvec x, arma::colvec y){
  double out = exp(-euclnormC(x - y));
  return (out);
}



// gaussian kernel matrix
// [[Rcpp::export]]
arma::mat gaussk_mat(const arma::mat& x, const arma::mat& y) {
  int m = x.n_cols;
  int n = x.n_rows;
  arma::mat kmat(m, m, arma::fill::zeros);
  
  // Precompute constants
  double factor = 1.0 / n;
  
  for (int i = 0; i < m; ++i) {
    for (int j = i; j < m; ++j) {
      // Compute the sum of Gaussian kernel values
      double sum = arma::accu(exp(-arma::square(x.col(i) - y.col(j))));
      
      // Assign to the upper triangle of kmat
      kmat(i, j) = sum * factor;
      
      // If i != j, assign the symmetric element
      if (i != j)
        kmat(j, i) = kmat(i, j);
    }
  }
  return kmat;
}

// mv gaussian kernel matrix
// [[Rcpp::export]]
arma::mat gaussk_mv_mat(const arma::mat& x, const arma::mat& y) {
  int m = x.n_cols;
  arma::mat kmat(m, m, arma::fill::zeros);

  for (int i = 0; i < m; ++i) {
    for (int j = i; j < m; ++j) {
      // Compute the sum of Gaussian kernel values
      kmat(i, j) = gaussk_mv(x.col(i), y.col(j));
      
      // If i != j, assign the symmetric element
      if (i != j)
        kmat(j, i) = kmat(i, j);
    }
  }
  return kmat;
}


// laplace kernel matrix
// [[Rcpp::export]]
arma::mat laplk_mat(const arma::mat& x, const arma::mat& y, double sigma) {
  int m = x.n_cols;
  int n = x.n_rows;
  arma::mat kmat(m, m, arma::fill::zeros);
  
  // Precompute constants
  double factor = 1.0 / n;
  
  for (int i = 0; i < m; ++i) {
    for (int j = i; j < m; ++j) {
      // Compute the sum of Laplace kernel values
      double sum = arma::accu(exp(abs(x.col(i) - y.col(j)) / sigma));
      
      // Assign to the upper triangle of kmat
      kmat(i, j) = sum * factor;
      
      // If i != j, assign the symmetric element
      if (i != j)
        kmat(j, i) = kmat(i, j);
    }
  }
  return kmat;
}

// mv laplace kernel matrix
// [[Rcpp::export]]
arma::mat laplk_mv_mat(const arma::mat& x, const arma::mat& y) {
  int m = x.n_cols;
  arma::mat kmat(m, m, arma::fill::zeros);
  
  for (int i = 0; i < m; ++i) {
    for (int j = i; j < m; ++j) {
      // Compute the sum of Laplace kernel values
      kmat(i, j) = laplk_mv(x.col(i), y.col(j));
      
      // If i != j, assign the symmetric element
      if (i != j)
        kmat(j, i) = kmat(i, j);
    }
  }
  return kmat;
}


// energy kernel matrix
// [[Rcpp::export]]
arma::mat enerk_mat(const arma::mat& x, const arma::mat& y) {
  int m = x.n_cols;
  int n = x.n_rows;
  arma::mat kmat(m, m, arma::fill::zeros);
  
  // Precompute constants
  double factor = 1.0 / n;
  
  for (int i = 0; i < m; ++i) {
    for (int j = i; j < m; ++j) {
      // Compute the sum of Energy kernel values
      double sum = arma::accu(abs(x.col(i)) + abs(y.col(j)) - abs(x.col(i) - y.col(j)));
      
      // Assign to the upper triangle of kmat
      kmat(i, j) = sum * factor;
      
      // If i != j, assign the symmetric element
      if (i != j)
        kmat(j, i) = kmat(i, j);
    }
  }
  return kmat;
}

// mv energy kernel matrix
// [[Rcpp::export]]
arma::mat enerk_mv_mat(const arma::mat& x, const arma::mat& y) {
  int m = x.n_cols;
  arma::mat kmat(m, m, arma::fill::zeros);
  
  for (int i = 0; i < m; ++i) {
    for (int j = i; j < m; ++j) {
      // Compute the sum of Energy kernel values
      kmat(i, j) = enerk_mv(x.col(i), y.col(j));
      
      // If i != j, assign the symmetric element
      if (i != j)
        kmat(j, i) = kmat(i, j);
    }
  }
  return kmat;
}
