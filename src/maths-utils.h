#ifndef MATHS_UTILS_H
#define MATHS_UTILS_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Fonctions exportées vers R
arma::mat mat_inv(const arma::mat& A);
arma::mat solve_linear_syst(const arma::mat& A, const arma::mat& B);
arma::mat kronecker_gamma_diag_mult(const arma::mat& gamma, const arma::vec& d_diag, int p);

// Fonctions internes (utilisées par d'autres cpp)
arma::vec rmvnorm_cpp(const arma::vec& mean, const arma::mat& cov);
double logdmvnorm_cpp(const arma::vec& x, const arma::vec& mean, const arma::mat& sigma);

#endif
