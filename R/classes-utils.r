#' Parse Right-Hand Side Terms of a Formula
#'
#' Internal utility that parses the right-hand side (RHS) of a model
#' formula and extracts included terms, excluded terms, and special
#' directives such as grouping or repeated-measurement variables.
#'
#' @param expr An R expression corresponding to the RHS of a formula.
#'
#' @return A named list with components:
#' \describe{
#'   \item{forced}{A character vector of terms explicitly included in the RHS.}
#'   \item{excluded}{A character vector of terms explicitly excluded
#'  from the RHS.}
#'   \item{id}{A character scalar identifying the grouping variable specified
#'   via \code{group()}, or \code{NULL} if absent.}
#'   \item{t}{A character scalar identifying the repeated-measurement variable
#'   specified via \code{repeated()}, or \code{NULL} if absent.}
#'   \item{has_dot}{A logical value indicating whether the symbol \code{.}
#'   was present in the RHS.}
#' }
#'
#' @details
#' The function recursively traverses the expression Abstract Syntax Tree
#' and classifies terms according to their sign:
#'
#' \itemize{
#'   \item Terms combined using \code{+} are treated as included.
#'   \item Terms combined using \code{-} are treated as excluded.
#'   \item Calls to \code{group()} identify a grouping variable.
#'   \item Calls to \code{repeated()} identify a repeated-measurement variable.
#'   \item The symbol \code{.} indicates inclusion of all candidate covariates.
#' }
#'
#' Duplicate terms are removed from the returned vectors.
#'
#' @keywords internal
#' @noRd
parse_rhs_terms <- function(expr) {
  state <- new.env(parent = emptyenv())

  state$forced <- character()
  state$excluded <- character()
  state$id <- NULL
  state$t <- NULL
  state$has_dot <- FALSE

  walk <- function(e, sign = "+") {
    if (is.call(e)) {
      op <- as.character(e[[1]])

      if (op == "+") {
        walk(e[[2]], sign)
        walk(e[[3]], sign)
        return()
      }

      if (op == "-") {
        walk(e[[2]], sign)
        walk(e[[3]], ifelse(sign == "+", "-", "+"))
        return()
      }

      if (op == "group") {
        state$id <- deparse(e[[2]])
        return()
      }

      if (op == "repeated") {
        state$t <- deparse(e[[2]])
        return()
      }

      term <- deparse(e)
    } else {
      term <- deparse(e)
    }

    if (term == ".") {
      state$has_dot <- TRUE
      return()
    }

    if (sign == "+") {
      state$forced <- c(state$forced, term)
    } else {
      state$excluded <- c(state$excluded, term)
    }
  }

  walk(expr)

  list(
    forced = unique(state$forced),
    excluded = unique(state$excluded),
    id = state$id,
    t = state$t,
    has_dot = state$has_dot
  )
}



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
#' @noRd
get_variables_from_formula <- function(f) {
  # ----- Extract LHS -----
  y <- deparse(f[[2]])

  # ----- Parse RHS via AST -----
  rhs_info <- parse_rhs_terms(f[[3]])

  forced <- rhs_info$forced
  excluded <- rhs_info$excluded
  id <- rhs_info$id
  t <- rhs_info$t

  # Convert empty -> NULL
  if (length(forced) == 0) forced <- NULL
  if (length(excluded) == 0) excluded <- NULL

  list(
    y = y,
    id = id,
    t = t,
    forced = forced,
    excluded = excluded
  )
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
#' @noRd

make_phi_fn <- function(g_user) {
  fmls <- names(formals(g_user))

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
#' @noRd
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
