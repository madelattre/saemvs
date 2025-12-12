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
