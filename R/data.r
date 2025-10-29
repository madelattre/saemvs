#' Example datasets for saemvs
#'
#' Simulated datasets used to demonstrate the usage of the `saemvs` package
#' for variable selection in nonlinear mixed-effects models.
#'
#' These datasets represent the same simulated data in two different formats:
#' - `example_data_list`: data stored as lists of individual responses and
#' time points.
#' - `example_data_df`: data stored as data frames, one for the
#' longitudinal measurements and one for the covariates.
#'
#' The simulated data correspond to 200 individuals, each measured over
#' time, with 200 candidate covariates for variable selection.
#' Only a subset of these covariates truly influences the model parameters.
#'
#' @docType data
#' @name example_data
#'
#' @usage data(example_data_list)
#' data(example_data_df)
#'
#' @format
#' **`example_data_list`**
#' A list with the following elements:
#' \describe{
#'   \item{y_list}{List of numeric vectors — observed responses for
#' each individual.}
#'   \item{t_list}{List of numeric vectors — observation time points for
#' each individual.}
#'   \item{x}{Numeric matrix (200 × 200) — candidate covariates,
#' rows = individuals, columns = covariates.}
#' }
#'
#' **`example_data_df`**
#' A list of two data frames:
#' \describe{
#'   \item{df_long}{Data frame with columns:
#'     \code{id} (individual ID), \code{time} (measurement time),
#' and \code{y} (observed response).}
#'   \item{covar_df}{Data frame in wide format, one row per individual,
#' containing 200 candidate covariates.}
#' }
#'

NULL


#' Example dataset: df_long
#'
#' Long-format simulated data for saemvs.
#' @docType data
#' @name df_long
NULL

#' Example dataset: df_cov
#'
#' Covariate matrix for saemvs.
#' @docType data
#' @name df_cov
NULL

#' Example dataset: example_data_list
#'
#' List-format simulated data for saemvs.
#' @docType data
#' @name example_data_list
NULL
