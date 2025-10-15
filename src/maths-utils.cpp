#include "maths-utils.h"

// // [[Rcpp::export]]
// arma::mat mat_inv(const arma::mat& A) {
//   return inv(A);
// }

// [[Rcpp::export]]
arma::mat mat_inv(const arma::mat& A) {
  try {
    // Check that the matrix contains only finite values
    // (no NA, NaN, or Inf). Armadillo's is_finite() does this efficiently.
    if (!A.is_finite()) {
      stop("Error: the estimated covariance matrix contains non-finite values (NA, NaN, or Inf).");
    }

    // Try to invert the matrix.
    // arma::inv() throws a std::runtime_error if the matrix is singular or not invertible.
    arma::mat invA = arma::inv(A);
    return invA;

  } catch (std::runtime_error &e) {
    // Typical case: the matrix is singular or ill-conditioned.
    stop("Error: the estimated covariance matrix is not invertible (probably singular or ill-conditioned).");
  } catch (...) {
    // Catch-all for any other unexpected errors.
    stop("Unknown error occurred while trying to invert the estimated covariance matrix.");
  }
}


// [[Rcpp::export]]
arma::mat solve_linear_syst(const arma::mat& A, const arma::mat& B) {
  return arma::solve(A, B, solve_opts::fast);  // résout A * X = B, solve_opts::fast pour une résolution rapide si A symétrique définie positive
}

// [[Rcpp::export]]
arma::mat kronecker_gamma_diag_mult(const arma::mat& gamma, const arma::vec& d_diag, int p) {
  int q = gamma.n_rows;  // ici 2
  int dim = p * q;       // ici 2 * 501 = 1002

  arma::mat res(dim, dim, fill::zeros);

  for (int i = 0; i < q; ++i) {
    for (int j = 0; j < q; ++j) {
      double gij = gamma(i, j);
      if (gij == 0.0) continue;

      int row_offset = i * p;
      int col_offset = j * p;

      for (int k = 0; k < p; ++k) {
        int row = row_offset + k;
        int col = col_offset + k;
        res(row, col) = gij * d_diag(row);
      }
    }
  }

  return res;
}


// arma::vec rmvnorm_cpp(const arma::vec& mean, const arma::mat& cov) {
//   int d = mean.n_elem;
//   arma::vec z = arma::randn(d);
//   return mean + arma::chol(cov, "lower") * z;
// }


// double logdmvnorm_cpp(
//     const arma::vec& x,
//     const arma::vec& mean,
//     const arma::mat& sigma
// ) {
//   int d = x.n_elem;
//   arma::vec diff = x - mean;
//   double quadform = arma::as_scalar(diff.t() * arma::inv(sigma) * diff);
//   double logdet_sigma;
//   double sign;
//   arma::log_det(logdet_sigma, sign, sigma);
//   return -0.5 * (d * std::log(2 * M_PI) + logdet_sigma + quadform);
// }
