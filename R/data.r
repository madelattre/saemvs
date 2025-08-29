#' Small example dataset for saemvs
#'
#' A small simulated dataset to demonstrate the usage of `saemvs`.
#' This dataset is intended to illustrate the usage of the `saemvs` package for variable selection in nonlinear mixed-effects models.
#'
#' @docType data
#' @name small_example_data
#' @usage data(small_example_data)
#' @format A data frame with 100 individuals and 4 covariates:
#' \describe{
#'   \item{y_list}{Responses. List of numeric vectors, each one containing observed responses for one individual.}
#'   \item{t_list}{Time points. List of numeric vectors, each one containing observation time points for one individual.}
#'   \item{x}{Matrix of covariates with dimensions 100 x 4 (rows = individuals, columns = covariates).}
#' }
NULL