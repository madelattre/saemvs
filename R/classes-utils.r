#' @title Extract variables and special declarations from a model formula
#'
#' @description
#' Internal utility that parses a model formula used in \code{saemvsData}
#' construction. The formula is expected to follow an extended syntax of the
#' form:
#' \preformatted{
#'   y ~ . + repeated(time) + group(id)
#' }
#'
#' In addition to the standard response and covariates, this function extracts:
#' \itemize{
#'   \item the response variable,
#'   \item the grouping (subject) variable via \code{group()},
#'   \item the repeated-measure (time) variable via \code{repeated()},
#'   \item covariates explicitly included or excluded in the formula.
#' }
#'
#' The dot (\code{.}) is detected but not expanded at this stage, as its
#' resolution requires access to the data. Expansion is deferred to the data
#' preprocessing stage.
#'
#' @param f A model formula specifying the response, covariates, and special
#'   terms \code{group()} and \code{repeated()}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{y}{Character scalar giving the response variable name.}
#'   \item{id}{Character scalar giving the subject/group identifier, or
#'     \code{NULL} if not specified.}
#'   \item{t}{Character scalar giving the time variable, or \code{NULL} if not
#'     specified.}
#'   \item{forced}{Character vector of covariates explicitly included in the
#'     formula, or \code{NULL}.}
#'   \item{excluded}{Character vector of covariates explicitly excluded via
#'     subtraction, or \code{NULL}.}
#' }
#'
#' @keywords internal
get_variables_from_formula <- function(f) {
  # Parse formula of form: y ~ . + repeated(time) + group(id)
  # [+ forced_vars] [- excluded_vars]

  # Extract LHS (response variable)
  y <- deparse(f[[2]])

  # Extract RHS as string for parsing
  rhs_expr <- f[[3]]

  # Initialize variables
  id <- NULL
  t <- NULL
  forced <- character(0)
  excluded <- character(0)
  has_dot <- FALSE

  # Helper function to extract function arguments
  extract_function_arg <- function(expr_str, func_name) {
    pattern <- paste0(func_name, "\\(([^)]+)\\)")
    match <- regexpr(pattern, expr_str)
    if (match > 0) {
      matched_str <- regmatches(expr_str, match)
      arg <- sub(paste0(func_name, "\\("), "", matched_str)
      arg <- sub("\\)", "", arg)
      return(trimws(arg))
    }
    return(NULL) # nolint: return_linter
  }

  # Remove spaces and convert to string
  rhs_string <- deparse(rhs_expr, width.cutoff = 500)
  rhs_string <- paste(rhs_string, collapse = " ")

  # Extract group(id)
  id <- extract_function_arg(rhs_string, "group")

  # Extract repeated(time)
  t <- extract_function_arg(rhs_string, "repeated")

  # Remove special function calls from RHS for further parsing
  rhs_cleaned <- gsub("group\\([^)]+\\)", "", rhs_string)
  rhs_cleaned <- gsub("repeated\\([^)]+\\)", "", rhs_cleaned)
  rhs_cleaned <- trimws(rhs_cleaned)

  # Clean up dangling operators (+ or - at the end or multiple consecutive)
  rhs_cleaned <- gsub("\\s*\\+\\s*\\+", "+", rhs_cleaned) # + + -> +
  rhs_cleaned <- gsub("\\s*-\\s*-", "-", rhs_cleaned) # - - -> -
  rhs_cleaned <- gsub("^\\s*[+]\\s*", "", rhs_cleaned) # Remove leading +
  rhs_cleaned <- gsub("\\s*[+]\\s*$", "", rhs_cleaned) # Remove trailing +
  rhs_cleaned <- gsub("\\s*[-]\\s*$", "", rhs_cleaned) # Remove trailing -
  # Clean up . followed by operators
  rhs_cleaned <- gsub("\\.\\s*\\+", ".", rhs_cleaned) # . + -> .
  rhs_cleaned <- gsub("\\.\\s*-", ".", rhs_cleaned) # . - -> .
  rhs_cleaned <- trimws(rhs_cleaned)

  # Check if there's a dot in the cleaned RHS
  has_dot <- grepl("\\.", rhs_cleaned)

  # Parse the remaining terms
  if (nchar(rhs_cleaned) > 0) {
    # If there's a . in the formula, we need to handle it specially
    # because terms() can't evaluate . without data context
    if (has_dot) {
      # Remove the . and parse the rest
      rhs_without_dot <- gsub("\\.", "", rhs_cleaned)
      # Clean up + +
      rhs_without_dot <- gsub("\\s*\\+\\s*\\+", "+", rhs_without_dot)
      # Remove leading +
      rhs_without_dot <- gsub("^\\s*[+]\\s*", "", rhs_without_dot)
      # Remove trailing +
      rhs_without_dot <- gsub("\\s*[+]\\s*$", "", rhs_without_dot)
      rhs_without_dot <- trimws(rhs_without_dot)

      if (nchar(rhs_without_dot) > 0) {
        temp_formula <- as.formula(paste("dummy ~", rhs_without_dot))
        trm <- terms(temp_formula)
        term_labels <- attr(trm, "term.labels")
      } else {
        term_labels <- character(0)
      }
    } else {
      # No dot, we can safely use terms()
      temp_formula <- as.formula(paste("dummy ~", rhs_cleaned))
      trm <- terms(temp_formula)
      term_labels <- attr(trm, "term.labels")
    }

    # Process each term
    for (term in term_labels) {
      if (startsWith(term, "-")) {
        # Excluded variable
        excluded <- c(excluded, sub("^-\\s*", "", term))
      } else {
        # Forced (explicit in formula)
        forced <- c(forced, term)
      }
    }
  }

  # If . is used, we need to expand it
  if (has_dot) {
    # Note: all.vars will include the special function names (group, repeated)
    # We need to ignore them and extract only actual variables
    # The expansion of . will be done by saemvsDataFromDFs which has access
    # to data
    # For now, we mark that we have a dot and set forced to indicate
    # "all except excluded"
    # This will be handled in saemvsDataFromDFs
    forced <- character(0)
    # Reset to be empty - actual variables will be determined from data
  }

  # Remove empty strings and convert to NULL if empty
  forced <- forced[forced != ""]
  if (length(forced) == 0) forced <- NULL

  excluded <- excluded[excluded != ""]
  if (length(excluded) == 0) excluded <- NULL

  return(list( # nolint: return_linter
    y = y,
    id = id,
    t = t,
    forced = forced,
    excluded = excluded
  ))
}

#' @title Standardize a user-defined model function to \code{function(t, phi)}
#'
#' @description
#' Internal helper that converts a user-supplied model function into a
#' standardized form with signature \code{function(t, phi)}.
#'
#' The function body is rewritten so that:
#' \itemize{
#'   \item the user-defined time variable is renamed to \code{t},
#'   \item individual model parameters are replaced by indexed elements of
#'     \code{phi}.
#' }
#'
#' This transformation allows the internal SAEMVS algorithms to operate on a
#' uniform model representation, independently of how the user specified the
#' model.
#'
#' @param g_user A user-defined model function, either of the form
#'   \code{function(t, p1, p2, ...)}.
#'
#' @return A function with signature \code{function(t, phi)} equivalent to the
#'   original model.
#'
#' @keywords internal

make_phi_fn <- function(g_user) {
  fmls <- names(formals(g_user))
  if (length(fmls) < 1) {
    stop("Function must have at least one argument (time)")
  }

  time_name <- fmls[1]
  param_names <- if (length(fmls) > 1) fmls[-1] else character(0)

  # Recursive AST transformer
  transform_expr <- function(expr) {

    # NULL / empty expressions
    if (is.null(expr) || length(expr) == 0) {
      return(expr)
    }

    # Atomic literals: numeric, logical, character, etc.
    if (is.atomic(expr)) {
      return(expr)
    }

    # Symbol (variable name)
    if (is.name(expr)) {
      name <- as.character(expr)

      if (name == time_name) {
        return(as.symbol("t"))
      }

      idx <- match(name, param_names)
      if (!is.na(idx)) {
        # Replace parameter by phi[idx]
        return(substitute(phi[I], list(I = as.numeric(idx))))
      }

      # Any other symbol: leave untouched
      return(expr)
    }

    # Call
    if (is.call(expr)) {
      op <- expr[[1]]

      # Special case: assignment
      if (as.character(op) %in% c("<-", "=", "<<-")) {
        # LHS must NOT be transformed
        lhs <- expr[[2]]
        rhs <- expr[[3]]
        return(as.call(list(
          op,
          lhs,
          transform_expr(rhs)
        )))
      }

      # General call: transform arguments only
      new_args <- lapply(as.list(expr)[-1], transform_expr)
      return(as.call(c(list(op), new_args)))
    }

    # Fallback (should not happen often)
    expr
  }

  # Transform body
  new_body <- transform_expr(body(g_user))

  # Build new function: function(t, phi)
  new_fn <- function(t, phi) {}
  body(new_fn) <- new_body
  environment(new_fn) <- .GlobalEnv
  new_fn

}


#' @title Validate forced covariate declarations for model parameters
#'
#' @description
#' Internal validation helper for the \code{x_forced_support} argument of
#' \code{saemvsModel}. This function checks that forced covariate declarations
#' are provided as a named list mapping model parameter names to covariate
#'  names.
#'
#' The function ensures that:
#' \itemize{
#'   \item the object is either \code{NULL} or a named list,
#'   \item all names correspond to valid model parameters,
#'   \item each list element is a character vector of covariate names.
#' }
#'
#' An empty list is allowed and indicates that no covariates are forced in the
#' model.
#'
#' @param x A named list or \code{NULL}. Each element corresponds to a model
#'   parameter and contains the names of covariates that must be forced for that
#'   parameter.
#' @param phi_names Character vector giving the valid model parameter names.
#'
#' @return The validated \code{x_forced_support} object (possibly an empty
#'  list).
#'
#' @keywords internal
validate_x_forced_support <- function(x, phi_names) {

  if (is.null(x)) {
    return(list())
  }

  if (!is.list(x)) {
    stop("'x_forced_support' must be NULL or a named list.")
  }

  if (length(x) == 0) {
    return(x)
  }

  if (is.null(names(x)) || any(names(x) == "")) {
    stop("'x_forced_support' must be a named list of model parameters.")
  }

  bad_phi <- setdiff(names(x), phi_names)
  if (length(bad_phi) > 0) {
    stop(
      "Unknown parameter(s) in 'x_forced_support': ",
      paste(bad_phi, collapse = ", ")
    )
  }

  for (nm in names(x)) {
    if (!is.character(x[[nm]])) {
      stop(
        "Element '", nm,
        "' of 'x_forced_support' must be a character vector of covariate names."
      )
    }
  }

  x
}
