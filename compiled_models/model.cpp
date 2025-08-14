
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
inline double g_scalar_cpp(const arma::vec& phi, double t) {
  return {    phi[3-1] + phi[1-1]/(1 + exp(-(t - phi[2-1])))};
}

// [[Rcpp::export]]
arma::vec g_vector_cpp(const arma::vec& phi, const arma::vec& t) {
  arma::vec out(t.n_elem);
  for (int i = 0; i < t.n_elem; i++) {
      out[i] = g_scalar_cpp(phi, t[i]);
  }
  return out;
}

// [[Rcpp::export]]
arma::vec rmvnorm_cpp(const arma::vec& mean, const arma::mat& sigma) {
  arma::vec z = arma::randn(mean.n_elem);
  return mean + arma::chol(sigma) * z;
}

// [[Rcpp::export]]
double logdmvnorm_cpp(const arma::vec& x, const arma::vec& mean, const arma::mat& sigma) {
  int k = x.n_elem;
  double det_sigma = arma::det(sigma);
  arma::vec diff = x - mean;
  double quadform = arma::as_scalar(diff.t() * arma::inv(sigma) * diff);
  return -0.5 * (k * std::log(2 * M_PI) + std::log(det_sigma) + quadform);
}

// [[Rcpp::export]]
List metropolis_vector_cpp(
    const List& y,
    const List& t,
    const List& phi_current,
    const List& mean_prop,
    const arma::mat& var_prop_mat,
    const double& sigma2,
    int niter_mh,
    const double& kappa,
    const std::string& kernel = "pop"
) {
    int n = y.size();
    List phi_chains(n);

    for (int i = 0; i < n; ++i) {
        NumericVector y_i = y[i];
        NumericVector t_i = t[i];
        int ni = y_i.size();

        arma::vec phi0 = as<arma::vec>(phi_current[i]);
        arma::vec mean_i = as<arma::vec>(mean_prop[i]);

        int d = phi0.n_elem;
        arma::mat chain(d, niter_mh + 1);
        chain.col(0) = phi0;

        for (int r = 0; r < niter_mh; ++r) {
            arma::vec phi_old = chain.col(r);
            arma::vec phi_prop;

            if (kernel == "pop") {
                phi_prop = rmvnorm_cpp(mean_i, var_prop_mat);
            } else {
                phi_prop = rmvnorm_cpp(phi_old, kappa * var_prop_mat);
            }

            double logratio = 0.0;
            double sd = std::sqrt(sigma2);

            for (int j = 0; j < ni; ++j) {
                double mean_new = g_scalar_cpp(phi_prop, t_i[j]);
                double mean_old = g_scalar_cpp(phi_old, t_i[j]);
                logratio += R::dnorm(y_i[j], mean_new, sd, true)
                          - R::dnorm(y_i[j], mean_old, sd, true);
            }

            if (kernel == "random_walk") {
                logratio += logdmvnorm_cpp(phi_prop, mean_i, var_prop_mat)
                          - logdmvnorm_cpp(phi_old, mean_i, var_prop_mat);
            }

            if (std::log(R::runif(0.0, 1.0)) <= logratio) {
                chain.col(r + 1) = phi_prop;
            } else {
                chain.col(r + 1) = phi_old;
            }
        }

        phi_chains[i] = chain.col(niter_mh);
    }

    return phi_chains;
}

