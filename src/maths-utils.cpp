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
    arma::mat invA;
    bool ok = arma::inv(invA, A);  // retourne false si A non inversible
    if (!ok) {
      Rcpp::stop("Covariance matrix is singular or ill-conditioned.");
    }
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
  arma::mat X;  // matrice resultat
  bool ok = arma::solve(X, A, B);  // version securisee

  if (!ok) {
    Rcpp::stop("Error: linear system cannot be solved (matrix singular or ill-conditioned).");
  }

  return X;
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
