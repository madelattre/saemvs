#' Compile a user-defined model into C++ with RcppArmadillo
#'
#' Internal function that compiles a user-specified R function into optimized C++
#' code using \pkg{RcppArmadillo}.  
#'
#' The function extracts the body of an R function \code{g_fun}, adapts the
#' indexing from R (1-based) to C++ (0-based), and generates C++ code containing:
#' \itemize{
#'   \item A scalar evaluation function \code{g_scalar_cpp}.
#'   \item A vectorized evaluation function \code{g_vector_cpp}.
#'   \item Utility functions for multivariate normal simulation and density.
#'   \item A Metropolis-Hastings sampler for parameter chains.
#' }
#' The generated code is compiled on the fly and loaded into the current R session.
#'
#' @param g_fun A user-defined R function of the form \code{g(phi, t)} representing
#'   the model function, where \code{phi} is a parameter vector and \code{t} a scalar
#'   or vector of inputs. This function body will be translated into C++.
#' @param build_dir Character string indicating the directory where compiled models
#'   should be stored. Defaults to \code{"compiled_models"}. Currently, the function
#'   writes the generated C++ code to a temporary file rather than persisting it in
#'   \code{build_dir}.
#'
#' @return
#' Invisibly returns \code{TRUE} if compilation succeeds. The compiled functions
#' \code{g_scalar_cpp}, \code{g_vector_cpp}, \code{rmvnorm_cpp},
#' \code{logdmvnorm_cpp}, and \code{metropolis_vector_cpp} are loaded and callable
#' within the R session.
#'
#' @details
#' The following steps are performed:
#' \enumerate{
#'   \item Extract the body of \code{g_fun} and collapse it into a single string.
#'   \item Replace R-style 1-based indexing \code{phi[i]} with C++-compatible
#'         0-based indexing.
#'   \item Embed the function body in a templated C++ source code.
#'   \item Write the code to a temporary \code{.cpp} file.
#'   \item Compile and load the code with \code{Rcpp::sourceCpp()}.
#' }
#'
#' This function is intended for internal use. Error handling is minimal and assumes
#' that the input function \code{g_fun} is well-formed and compatible with the
#' generated C++ code.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' g_fun <- function(phi, t) { phi[1] + phi[2] * t }
#' compile_model(g_fun)
#'
#' # Then call the compiled function
#' g_scalar_cpp(c(1, 2), 0.5)
#' }
compile_model <- function(g_fun, build_dir = "compiled_models") {
  
  # Extract the body of the R function
  body_txt <- deparse(body(g_fun))
  body_txt <- paste(body_txt, collapse = "")
  
  # Adjustment of R indexing (1-based) to C++ (0-based)
  body_txt <- gsub("phi\\[(\\d+)\\]", "phi[\\1-1]", body_txt)

  # Generate the complete C++ code
  code_cpp <- sprintf('
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
inline double g_scalar_cpp(const arma::vec& phi, double t) {
  return %s;
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
', body_txt)

  # Creates a temporary file for C++ code
  tmp_cpp <- tempfile(fileext = ".cpp")
  writeLines(code_cpp, tmp_cpp)

  # Compiles and loads the code into memory
  Rcpp::sourceCpp(tmp_cpp, rebuild = TRUE, verbose = FALSE)

  invisible(TRUE)
}



#' Load a compiled shared object model
#'
#' Internal function that dynamically loads a compiled shared object (\code{.so})
#' file containing a previously compiled model. This provides direct access to
#' the C++ functions generated by \code{\link{compile_model}}.
#'
#' @param so_path Character string giving the path to the shared object
#'   (\code{.so}) file to be loaded.
#'
#' @return
#' Returns \code{TRUE} invisibly if the shared object is successfully loaded.
#'
#' @details
#' This function is a thin wrapper around \code{dyn.load()}, and is intended for
#' internal use.  
#' It does not check whether the symbols are already loaded or whether the file
#' was compiled with \code{\link{compile_model}}.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' so_file <- "path/to/compiled_model.so"
#' load_compiled_model(so_file)
#' }
load_compiled_model <- function(so_path) {
  dyn.load(so_path)
}
