#' Example datasets for saemvs
#'
#' Simulated datasets used to demonstrate the usage of the `saemvs` package
#' for variable selection in nonlinear mixed-effects models.
#'
#' The simulated data correspond to 200 individuals, each measured over
#' time, with 200 candidate covariates for variable selection.
#' Only a subset of these covariates truly influences the model parameters.
#'
#' The data are provided in data-frame format:
#' - `df_long`: longitudinal measurements in long format
#' - `df_cov`: covariates in wide format (one row per individual)
#'
#' @docType data
#' @name example_data
#'
#' @usage data(df_long)
#' data(df_cov)
#'
#' @format
#' **`df_long`**
#' A data frame with columns:
#' \describe{
#'   \item{id}{Individual identifier}
#'   \item{time}{Measurement time}
#'   \item{y}{Observed response}
#' }
#'
#' **`df_cov`**
#' A data frame in wide format with:
#' \describe{
#'   \item{id}{Individual identifier}
#'   \item{X1, ..., X200}{Numeric covariates used for variable selection}
#' }
#'
NULL


#' Example dataset: df_long
#'
#' Long-format simulated longitudinal data for `saemvs`.
#'
#' Each row corresponds to one observation for one individual at a given
#' time point.
#'
#' @docType data
#' @name df_long
#' @usage data(df_long)
#' @format A data frame with columns \code{id}, \code{time}, and \code{y}.
NULL


#' Example dataset: df_cov
#'
#' Simulated covariate dataset for `saemvs`.
#'
#' Contains 200 candidate covariates in wide format, with one row per
#' individual.
#'
#' @docType data
#' @name df_cov
#' @usage data(df_cov)
#' @format A data frame with one row per individual and 200 covariate columns.
NULL
