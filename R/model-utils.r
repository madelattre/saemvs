# R to C++ transpiler using AST
# Converts an R function to C++ code compilable with RcppArmadillo


#' Transpile R Function to C++ Code
#'
#' Converts an R function into equivalent C++ code that can be used with Rcpp.
#' Supports mathematical functions, control flow structures, and recursive
#' transpilation of nested functions.
#'
#' @param g_fun A function object to be transpiled to C++
#' @param function_name Character string specifying the name of the C++
#'  function.
#'   Defaults to "g_scalar_cpp"
#' @param is_main Logical. If TRUE, adds Rcpp::export attribute to the function.
#'   Defaults to TRUE
#' @param hierarchy_path Character string representing the current position in
#'   the function hierarchy for nested function naming. Used internally for
#'   recursive calls. Defaults to NULL
#'
#' @return A character string containing the complete transpiled C++ function
#'   code, including function signature, variable declarations, translated
#'   expressions, and any auxiliary functions
#'
#' @details
#' The function supports:
#' \itemize{
#'   \item Mathematical functions (exp, log, sin, cos, etc.)
#'   \item Binary operators (+, -, *, /, ^, ==, <, >, &&, ||, etc.)
#'   \item Unary operators (-, +, !)
#'   \item Control flow (if/else, for loops, while loops)
#'   \item Variable declarations and assignments
#'   \item Function calls and nested function transpilation
#'   \item Vector indexing with automatic 1-to-0-based conversion
#'   \item Type inference for environment variables
#' }
#'
#' Warnings are accumulated during transpilation and reported to the user
#' for unsupported constructs or automatic conversions.
#'
#' @examples
#' \dontrun{
#' # Simple scalar function
#' f <- function(x, y) {
#'   z <- x + y
#'   return(exp(z))
#' }
#' cpp_code <- transpile_to_cpp(f, "my_function")
#' cat(cpp_code)
#' }
#'
#' @keywords internal
transpile_to_cpp <- function(
    g_fun, function_name = "g_scalar_cpp",
    is_main = TRUE,
    hierarchy_path = NULL) {
  # Context to track variables and code generation
  context <- new.env()
  context$variables <- list() # Declared variables
  context$param_names <- names(formals(g_fun)) # Parameter names
  context$code_lines <- c() # Generated C++ code lines
  context$declarations <- c() # Variable declarations
  context$temp_counter <- 0 # Counter for temporary variables
  context$warnings <- c() # Accumulated warnings
  context$indent_level <- 1 # Current indentation level (1: function body level)
  context$aux_functions_code <- list() # Map from C++ funct. names to their code
  context$aux_func_counter <- 0 # Counter for aux. function names at this level
  context$func_env <- environment(g_fun) # Environment of the original function
  context$hierarchy_path <- hierarchy_path # Current position in hierarchy (e.g., "1", "1_2", "2_1") # nolint : line_length_linter

  # Supported mathematical functions (R -> C++)
  math_functions <- list(
    "exp" = "std::exp",
    "log" = "std::log",
    "log10" = "std::log10",
    "sqrt" = "std::sqrt",
    "abs" = "std::abs",
    "fabs" = "std::fabs",
    "sin" = "std::sin",
    "cos" = "std::cos",
    "tan" = "std::tan",
    "asin" = "std::asin",
    "acos" = "std::acos",
    "atan" = "std::atan",
    "atan2" = "std::atan2",
    "sinh" = "std::sinh",
    "cosh" = "std::cosh",
    "tanh" = "std::tanh",
    "floor" = "std::floor",
    "ceiling" = "std::ceil",
    "round" = "std::round",
    "trunc" = "std::trunc",
    "max" = "std::max",
    "min" = "std::min",
    "pow" = "std::pow",
    "fmod" = "std::fmod",
    "gamma" = "R::gammafn",
    "lgamma" = "R::lgammafn",
    "digamma" = "R::digamma",
    "trigamma" = "R::trigamma",
    "is.nan" = "std::isnan",
    "is.na" = "R_IsNA",
    "sum" = "arma::sum",
    "c" = "arma::vec",
    "is.infinite" = "std::isinf"
  )

  # Supported binary operators
  binary_ops <- list(
    "+" = "+",
    "-" = "-",
    "*" = "*",
    "/" = "/",
    "^" = "std::pow",
    "%%" = "std::fmod",
    "%/%" = "std::floor", # Integer division approximation
    "==" = "==",
    "!=" = "!=",
    "<" = "<",
    "<=" = "<=",
    ">" = ">",
    ">=" = ">=",
    "&" = "&&",
    "|" = "||"
  )

  # Unary operators
  unary_ops <- list(
    "-" = "-",
    "+" = "+",
    "!" = "!"
  )

  # Generate a temporary variable name
  get_temp_var <- function() {
    context$temp_counter <- context$temp_counter + 1
    return(paste0("temp_", context$temp_counter)) # nolint: return_linter
  }

  # Infer type and generate initialization code for environment variables
  # infer_and_init_env_var <- function(var_name) {
  #   # Check if variable exists in function environment
  #   if (!is.null(context$func_env)) {
  #     if (exists(var_name, envir = context$func_env, inherits = TRUE)) { # modif TRUE/FALSE
  #       var_value <- get(var_name, envir = context$func_env, inherits = TRUE) # modif TRUE/FALSE

  #       # Infer type based on value
  #       if (is.numeric(var_value)) {
  #         if (length(var_value) == 1) {
  #           # Scalar
  #           return(list(
  #             type = "double",
  #             declaration = paste0(
  #               "double ",
  #               var_name,
  #               " = ",
  #               as.character(var_value),
  #               ";"
  #             )
  #           ))
  #         } else {
  #           # Vector
  #           values_str <- paste(var_value, collapse = ", ")
  #           return(list(
  #             type = "arma::vec",
  #             declaration = paste0(
  #               "arma::vec ",
  #               var_name,
  #               " = arma::vec({",
  #               values_str, "});"
  #             )
  #           ))
  #         }
  #       } else if (is.logical(var_value)) {
  #         if (length(var_value) == 1) {
  #           return(list(
  #             type = "bool",
  #             declaration = paste0(
  #               "bool ",
  #               var_name,
  #               " = ",
  #               ifelse(var_value, "true", "false"),
  #               ";"
  #             )
  #           ))
  #         } else {
  #           # Logical vector - convert to numeric
  #           values_str <- paste(as.numeric(var_value), collapse = ", ")
  #           return(list(
  #             type = "arma::vec",
  #             declaration = paste0(
  #               "arma::vec ",
  #               var_name,
  #               " = arma::vec({",
  #               values_str, "});"
  #             )
  #           ))
  #         }
  #       } else if (is.character(var_value)) {
  #         # Character - not fully supported, but we'll try
  #         if (length(var_value) == 1) {
  #           return(list(
  #             type = "std::string",
  #             declaration = paste0(
  #               "std::string ",
  #               var_name,
  #               " = \"",
  #               var_value,
  #               "\";"
  #             )
  #           ))
  #         } else {
  #           # Character vector - not supported,
  #           # use first element
  #           add_warning(paste(
  #             "Character vector",
  #             var_name,
  #             "truncated to first element"
  #           ))
  #           return(list(
  #             type = "std::string",
  #             declaration = paste0(
  #               "std::string ",
  #               var_name,
  #               " = \"",
  #               var_value[1],
  #               "\";"
  #             )
  #           ))
  #         }
  #       } else {
  #         # Unknown type - default to double
  #         add_warning(paste(
  #           "Unknown type for variable",
  #           var_name,
  #           "- defaulting to double"
  #         ))
  #         return(list(
  #           type = "double",
  #           declaration = paste0("double ", var_name, " = 0.0;")
  #         ))
  #       }
  #     }
  #   }
  #   # Variable not found in environment - default to double
  #   return(list( # nolint: return_linter
  #     type = "double",
  #     declaration = paste0("double ", var_name, ";")
  #   ))
  # }



  declare_variable <- function(var_name, init_value = NULL) {
    if (!var_name %in% context$variables &&
          !var_name %in% context$param_names) {
      context$variables <- c(context$variables, var_name)

      # Déterminer type et déclaration selon init_value
      if (!is.null(init_value)) {
        if (is.numeric(init_value) && length(init_value) == 1) {
          type <- "double"
          decl <- paste0("double ", var_name, " = ", init_value, ";")
        } else if (is.logical(init_value) && length(init_value) == 1) {
          type <- "bool"
          decl <- paste0(
            "bool ", var_name, " = ", ifelse(init_value,
              "true", "false"
            ),
            ";"
          )
        } else {
          type <- "double"
          decl <- paste0("double ", var_name, ";")
          add_warning(paste("Unknown type for variable", 
                            var_name, "- defaulting to double"))
        }
      } else {
        # Pas d'init_value : fallback double
        type <- "double"
        decl <- paste0("double ", var_name, ";")
      }

      context$declarations <- c(context$declarations, decl)
      return(type)
    }
    return(NULL)
  }

  # Add a code line with proper indentation
  add_code_line <- function(line) {
    indent <- paste(rep("  ", context$indent_level), collapse = "")
    context$code_lines <- c(context$code_lines, paste0(indent, line))
  }

  # Check if an expression needs parentheses (is a complex expression)
  needs_parens <- function(expr_str) {
    expr_str <- trimws(expr_str)
    # If already parenthesized, don't add more
    if (grepl("^\\s*\\(.*\\)\\s*$", expr_str)) {
      return(FALSE)
    }
    # Simple expressions don't need parentheses:
    # variables, numbers, function calls, indexing
    # Match: variable name, number, function call (name(...)),
    # or indexing (name[...])
    if (grepl("^[a-zA-Z_][a-zA-Z0-9_]*(\\[|\\()", expr_str) ||
          grepl("^[0-9]+\\.?[0-9]*([eE][+-]?[0-9]+)?$", expr_str) ||
          grepl("^[a-zA-Z_][a-zA-Z0-9_]*$", expr_str)) {
      return(FALSE)
    }
    # Otherwise, it's a complex expression that needs parentheses
    return(TRUE) # nolint: return_linter
  }

  # Main recursive translation function
  translate_expr <- function(expr) {
    # Case: NULL or empty expression
    if (is.null(expr) || length(expr) == 0) {
      return("")
    }

    # Case: NA of any type (character, logical, etc.) - treat as numeric NA
    # Only check for NA on atomic values,
    # not on symbols (variable names) or calls
    if (length(expr) == 1 && !is.name(expr) &&
          !is.call(expr) && is.atomic(expr)) {
      # Check if it's NA (but not NaN) -
      # use suppressWarnings to avoid errors on symbols
      if (!is.nan(expr)) {
        na_result <- suppressWarnings(tryCatch(is.na(expr),
          error = function(e) FALSE
        ))
        if (isTRUE(na_result)) {
          return("R_NaReal")
        }
      }
    }

    # Case: numeric literal
    if (is.numeric(expr)) {
      if (length(expr) == 1) {
        if (is.infinite(expr)) {
          if (expr > 0) {
            return("std::numeric_limits<double>::infinity()")
          } else {
            return("-std::numeric_limits<double>::infinity()")
          }
        }
        # Check for NA first (is.na returns TRUE for both NA and NaN)
        if (is.na(expr)) {
          # If it's specifically a NaN (not NA), use NaN
          if (is.nan(expr)) {
            return("std::numeric_limits<double>::quiet_NaN()")
          } else {
            # It's an NA, use R_NaReal so R recognizes it as NA
            return("R_NaReal")
          }
        }
        return(as.character(expr))
      } else {
        # Vector - create a temporary variable
        var_name <- get_temp_var()
        declare_variable(var_name)
        values <- paste(expr, collapse = ", ")
        add_code_line(paste0(var_name, " = arma::vec({", values, "});"))
        return(var_name)
      }
    }

    # Case: character string
    if (is.character(expr)) {
      # For now, we don't really support strings in numeric calculations
      add_warning("Character string detected, conversion not supported")
      return('""')
    }

    # Case: variable name
    if (is.name(expr)) {
      var_name <- as.character(expr)
      # If it's a parameter, return it as is
      if (var_name %in% context$param_names) {
        return(var_name)
      }
      # Otherwise, declare the variable if necessary
      declare_variable(var_name)
      return(var_name)
    }

    # Case: function call or operator
    if (is.call(expr)) {
      operator <- as.character(expr[[1]])
      args <- if (length(expr) > 1) expr[-1] else list()

      # Assignment (<- or =)
      if (operator %in% c("<-", "=", "<<-")) {
        if (length(args) < 2) {
          add_warning(paste("Invalid assignment:", deparse(expr)))
          stop("Syntax error in assignment")
        }
        var_name <- as.character(args[[1]])
        value <- translate_expr(args[[2]])
        declare_variable(var_name)
        add_code_line(paste0(var_name, " = ", value, ";"))
        return(var_name)
      }

      # Conditional if structure
      if (operator == "if") {
        if (length(args) < 1) {
          add_warning("Invalid if structure")
          stop("Error in if structure")
        }
        condition <- translate_expr(args[[1]])

        add_code_line(paste0("if (", condition, ") {"))

        # Increase indentation for if block
        context$indent_level <- context$indent_level + 1

        # Process if block (only translate once)
        if (length(args) >= 2) {
          translate_expr(args[[2]])
        }

        # Decrease indentation
        context$indent_level <- context$indent_level - 1

        # Process else block if present
        if (length(args) >= 3) {
          add_code_line("} else {")
          context$indent_level <- context$indent_level + 1
          translate_expr(args[[3]])
          context$indent_level <- context$indent_level - 1
        }

        add_code_line("}")
        return("")
      }

      # For loop
      if (operator == "for") {
        if (length(args) < 2) {
          add_warning("Invalid for loop")
          stop("Error in for loop")
        }
        var_name <- as.character(args[[1]])
        range_expr <- args[[2]]
        body_expr <- if (length(args) >= 3) args[[3]] else NULL

        # Process range (usually 1:n or seq(...))
        range_code <- translate_range(range_expr)

        add_code_line(paste0(
          "for (int ", var_name, " = ", range_code$start, "; ",
          var_name, " <= ", range_code$end, "; ", var_name, "++) {"
        ))

        # Increase indentation for loop body
        context$indent_level <- context$indent_level + 1

        # Process body
        if (!is.null(body_expr)) {
          translate_expr(body_expr)
        }

        # Decrease indentation
        context$indent_level <- context$indent_level - 1
        add_code_line("}")

        return("")
      }

      # While loop
      if (operator == "while") {
        if (length(args) < 1) {
          add_warning("Invalid while loop")
          stop("Error in while loop")
        }
        condition <- translate_expr(args[[1]])
        body_expr <- if (length(args) >= 2) args[[2]] else NULL

        add_code_line(paste0("while (", condition, ") {"))

        # Increase indentation for loop body
        context$indent_level <- context$indent_level + 1

        if (!is.null(body_expr)) {
          translate_expr(body_expr)
        }

        # Decrease indentation
        context$indent_level <- context$indent_level - 1
        add_code_line("}")

        return("")
      }

      # Block {}
      if (operator == "{") {
        if (length(args) > 0 && !is.null(args)) {
          last_result <- ""
          for (i in seq_along(args)) {
            # Translate expression (may add code lines or return a value)
            expr_result <- translate_expr(args[[i]])

            # Keep track of last non-empty result for implicit return
            # Only expressions that return a value
            # (not assignments, if, etc.) count
            if (expr_result != "") {
              last_result <- expr_result
            }
          }
          # Return the last expression result for implicit return
          return(last_result)
        }
        return("")
      }

      # Return
      if (operator == "return") {
        if (length(args) >= 1) {
          value <- translate_expr(args[[1]])
          add_code_line(paste0("return ", value, ";"))
        } else {
          add_code_line("return;")
        }
        return("")
      }

      # Explicit parentheses - just return the inner expression
      if (operator == "(") {
        if (length(args) >= 1) {
          return(translate_expr(args[[1]]))
        } else {
          add_warning("Empty parentheses")
          return("")
        }
      }


      # Indexing [
      if (operator == "[" || operator == "[[" || operator == "$") {
        if (length(args) < 2) {
          add_warning("Invalid indexing")
          stop("Error in indexing")
        }
        obj <- translate_expr(args[[1]])
        index <- args[[2]]

        # Process index (can be an expression)
        if (is.numeric(index)) {
          # 0-based adjustment
          index_cpp <- as.character(as.integer(index) - 1)
        } else {
          # Computed index - create temporary variable
          index_var <- get_temp_var()
          declare_variable(index_var)
          index_expr <- translate_expr(index)
          add_code_line(paste0(
            "int ",
            index_var,
            " = (int)(",
            index_expr,
            ") - 1;"
          ))
          index_cpp <- index_var
        }

        return(paste0(obj, "[", index_cpp, "]"))
      }

      # Binary operators
      # if (operator %in% names(binary_ops)) {
      if (operator %in% names(binary_ops) && length(args) == 2) {
        if (length(args) < 2) {
          add_warning(paste(
            "Binary operator",
            operator,
            "with less than 2 arguments"
          ))
          stop(paste("Error with operator", operator))
        }
        left <- translate_expr(args[[1]])
        right <- translate_expr(args[[2]])
        op_cpp <- binary_ops[[operator]]

        if (operator == "^") {
          # Power requires std::pow
          return(paste0(op_cpp, "(", left, ", ", right, ")"))
        } else if (operator == "%/%") {
          # Integer division
          return(paste0(op_cpp, "(", left, " / ", right, ")"))
        } else {
          # Add parentheses only if needed
          left_paren <- if (needs_parens(left)) paste0("(", left, ")") else left
          right_paren <- if (needs_parens(right)) paste0("(", right, ")") else right # nolint: line_length_linter
          return(paste0(left_paren, " ", op_cpp, " ", right_paren))
        }
      }

      # Unary operators
      if (operator %in% names(unary_ops)) {
        if (length(args) < 1) {
          add_warning(paste("Unary operator", operator, "without argument"))
          stop(paste("Error with unary operator", operator))
        }
        operand <- translate_expr(args[[1]])
        op_cpp <- unary_ops[[operator]]
        return(paste0(op_cpp, "(", operand, ")"))
      }

      # Mathematical functions
      if (operator %in% names(math_functions)) {
        if (length(args) < 1) {
          add_warning(paste("Function", operator, "without argument"))
          stop(paste("Error with function", operator))
        }
        func_cpp <- math_functions[[operator]]
        translated_args <- sapply(args, translate_expr)
        args_str <- paste(translated_args, collapse = ", ")

        # Special handling for is.nan: test both NaN and NA
        if (operator == "is.nan") {
          # Test both std::isnan and R_IsNA
          arg_expr <- translated_args[1]
          return(paste0(
            "(std::isnan(",
            arg_expr,
            ") || R_IsNA(",
            arg_expr,
            "))"
          ))
        }

        # Special handling for certain functions
        if (operator == "max" || operator == "min") {
          if (length(args) == 2) {
            return(paste0(func_cpp, "(", args_str, ")"))
          } else {
            # More than 2 arguments - create temporary variable
            var_name <- get_temp_var()
            declare_variable(var_name)
            add_code_line(paste0(var_name, " = ", translated_args[1], ";"))
            for (i in 2:length(translated_args)) {
              add_code_line(paste0(
                var_name, " = ", func_cpp, "(", var_name, ", ",
                translated_args[i], ");"
              ))
            }
            return(var_name)
          }
        }

        return(paste0(func_cpp, "(", args_str, ")"))
      }

      # Unrecognized generic function
      # Try to find it in the function's environment
      func_found <- FALSE
      aux_func_name <- NULL

      if (!is.null(context$func_env)) {
        # Try to find the function in the environment
        if (exists(operator, envir = context$func_env, inherits = FALSE)) {
          func_obj <- get(operator, envir = context$func_env, inherits = FALSE)
          # Check if it's a function
          if (is.function(func_obj)) {
            func_found <- TRUE
            # Always create a new function - no deduplication
            # Generate hierarchical name based on current position
            context$aux_func_counter <- context$aux_func_counter + 1

            # Build hierarchical path: current_path + new_counter
            if (is.null(context$hierarchy_path)) {
              # Top level: just the counter
              new_hierarchy_path <- as.character(context$aux_func_counter)
            } else {
              # Nested: append counter to existing path
              new_hierarchy_path <- paste0(
                context$hierarchy_path,
                "_",
                context$aux_func_counter
              )
            }

            # Generate name: aux_<hierarchy_path>
            aux_func_name <- paste0("aux_", new_hierarchy_path)

            # Transpile the auxiliary function recursively (not main function)
            # Pass the new hierarchy path
            aux_cpp_code <- transpile_to_cpp(func_obj, aux_func_name,
              is_main = FALSE,
              hierarchy_path = new_hierarchy_path
            )

            # Store the transpiled function
            context$aux_functions_code[[aux_func_name]] <- aux_cpp_code
          }
        }
      }

      if (func_found && !is.null(aux_func_name)) {
        # Use the transpiled auxiliary function
        translated_args <- sapply(args, translate_expr)
        args_str <- paste(translated_args, collapse = ", ")
        return(paste0(aux_func_name, "(", args_str, ")"))
      } else {
        # Function not found, attempt direct translation
        add_warning(paste(
          "Unrecognized function:",
          operator,
          "- attempting direct translation"
        ))
        translated_args <- sapply(args, translate_expr)
        args_str <- paste(translated_args, collapse = ", ")
        return(paste0(operator, "(", args_str, ")"))
      }
    }

    # Unhandled case
    add_warning(paste("Unhandled expression:", deparse(expr)))
    stop(paste("Unable to transpile expression:", deparse(expr)))
  }

  # Translate a range (for for loops)
  translate_range <- function(range_expr) {
    if (is.call(range_expr)) {
      op <- as.character(range_expr[[1]])
      if (op == ":") {
        # 1:n
        start <- translate_expr(range_expr[[2]])
        end <- translate_expr(range_expr[[3]])
        return(list(start = start, end = end))
      } else if (op == "seq" || op == "seq_len" || op == "seq_along") {
        add_warning("seq() in for loop - simplified conversion")
        if (length(range_expr) >= 2) {
          end <- translate_expr(range_expr[[2]])
          return(list(start = "1", end = end))
        }
      }
    }
    # Default
    return(list(start = "1", end = translate_expr(range_expr))) # nolint: return_linter
  }

  # Add a warning
  add_warning <- function(msg) {
    context$warnings <- c(context$warnings, msg)
    warning(msg)
  }

  # Main body
  body_expr <- body(g_fun)

  # Display accumulated warnings before starting
  if (length(context$warnings) > 0) {
    for (w in context$warnings) {
      warning(w)
    }
  }

  # Translate function body
  result <- translate_expr(body_expr)

  # If the last element is not a return, add the return
  if (length(context$code_lines) == 0 ||
        !grepl("^return", context$code_lines[length(context$code_lines)])) {
    if (result != "") {
      context$code_lines <- c(
        context$code_lines,
        paste0("return ", result, ";")
      )
    } else {
      # Last expression of the block
      context$code_lines <- c(context$code_lines, "return 0.0;")
    }
  }

  # Build C++ function signature
  param_decls <- c()
  for (param in context$param_names) {
    # By default, assume double for scalars
    # We could detect vectors, but for now assume arma::vec
    if (param == "phi") {
      param_decls <- c(param_decls, "const arma::vec& phi")
    } else {
      param_decls <- c(param_decls, paste0("double ", param))
    }
  }
  signature <- paste0(
    "inline double ", function_name, "(",
    paste(param_decls, collapse = ", "), ")"
  )

  # Assemble C++ code
  # Note: code_lines already have proper indentation from add_code_line
  declarations_str <- if (length(context$declarations) > 0) {
    paste("  ", context$declarations, collapse = "\n")
  } else {
    ""
  }

  code_str <- if (length(context$code_lines) > 0) {
    paste(context$code_lines, collapse = "\n")
  } else {
    "  return 0.0;"
  }

  # Add Rcpp::export only for main function
  export_attr <- if (is_main) "// [[Rcpp::export]]\n" else ""

  cpp_code <- paste0(
    export_attr, signature, " {\n",
    if (declarations_str != "") paste0(declarations_str, "\n") else "",
    code_str, "\n}"
  )

  # Prepend auxiliary functions if any
  if (length(context$aux_functions_code) > 0) {
    aux_functions_str <- paste(context$aux_functions_code, collapse = "\n\n")
    cpp_code <- paste0(aux_functions_str, "\n\n", cpp_code)
  }

  # Check if there were fatal errors
  if (length(context$warnings) > 0) {
    warning("Warnings were emitted during transpilation")
  }

  return(cpp_code)
}



#' Build R Backend for Model Functions
#'
#' Creates a list of backend functions for model evaluation and
#' Metropolis-Hastings sampling in R, providing scalar and vectorized
#' versions of a model function along with multivariate normal sampling
#' and density evaluation utilities.
#'
#' @param g_fun A function that takes parameters `phi` and time `t` and returns
#'   a scalar numeric value representing the model prediction.
#'
#' @return A list containing the following functions:
#'   \item{g_scalar}{Wrapper around `g_fun` that evaluates the model for a
#'  single time point.}
#'   \item{g_vector}{Vectorized version of the model function that evaluates
#'     `g_fun` across multiple time points using `vapply`.}
#'   \item{rmvnorm}{Generates random samples from a multivariate normal
#'  distribution using `MASS::mvrnorm`. Requires the MASS package.}
#'   \item{rmvnorm_mat(mean, sigma, n)}{Generates `n` draws from a
#' multivariate normal distribution.}
#'   \item{logdmvnorm}{Computes the log-density of a multivariate normal
#'  distribution given a mean vector and covariance matrix.}
#'   \item{metropolis_vector}{Implements a Metropolis-Hastings algorithm for
#'  sampling from the posterior distribution of parameters. Supports
#' population and random walk proposal kernels. Returns a list of sampled
#'  parameter chains.}
#'   \item{ll(yi_list, ti_list, phi_samples, sigma2)}{Computes the
#'  log-likelihood for a list of observed responses `yi_list` over times
#'  `ti_list`, given sampled parameters `phi_samples` and variance `sigma2`.}
#'
#' @details
#' The `metropolis_vector` function performs adaptive Metropolis-Hastings
#'  sampling for multiple individuals in parallel, accepting or rejecting
#'  proposed parameter values based on the log-likelihood ratio and proposal
#'  density ratio (for random walk kernel).
#'
#' @examples
#' \dontrun{
#' g_model <- function(phi, t) phi[1] * exp(-phi[2] * t)
#' backend <- build_backend_r(g_model)
#' samples <- backend$g_vector(c(1, 0.5), c(0, 1, 2, 3))
#' }
#'
#' @keywords internal
build_backend_r <- function(g_fun) {
  g_scalar <- function(t, phi) {
    if (is.matrix(phi) && nrow(phi) > 1) phi <- as.vector(phi)
    g_fun(t, phi)
  }


  g_vector <- function(t, phi) {
    vapply(t, function(ti) g_fun(ti, phi), numeric(1))
  }

  rmvnorm <- function(mean, sigma) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package MASS is required for R fallback.")
    }
    as.numeric(MASS::mvrnorm(1, mu = mean, Sigma = sigma))
  }

  rmvnorm_mat <- function(mean, sigma, n) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package MASS is required for R fallback.")
    }
    MASS::mvrnorm(n = n, mu = mean, Sigma = sigma)
  }

  logdmvnorm <- function(x, mean, sigma) {
    k <- length(x)
    diff <- x - mean
    quadform <- drop(t(diff) %*% solve(sigma) %*% diff)
    -0.5 * (k * log(2 * pi) + log(det(sigma)) + quadform)
  }

  loglik <- function(yi_list, ti_list, phi_samples, sigma2) {
    n <- length(yi_list)
    loglike <- 0
    for (i in seq_len(n)) {
      yi <- yi_list[[i]]
      ti <- ti_list[[i]]
      ni <- length(yi)
      n_samples <- nrow(phi_samples)
      contribution <- 0
      for (s in seq_len(n_samples)) {
        phi_i <- phi_samples[s, ]
        sumsq <- sum((yi -
                        vapply(ti, function(tj) g_scalar(tj, phi_i),
                               numeric(1)))^2)
        contribution <- contribution + exp(-sumsq / (2 * sigma2))
      }
      loglike <- loglike +
        log((2 * pi * sigma2)^(-ni/2) * contribution / n_samples)
    }
    loglike
  }

  metropolis_vector <- function(y, t, phi_current, mean_prop, var_prop_mat,
                                sigma2, niter_mh, kappa, kernel = "pop") {
    n <- length(y)
    phi_chains <- vector("list", n)
    sd <- sqrt(sigma2)

    for (i in seq_len(n)) {
      y_i <- y[[i]]
      t_i <- t[[i]]
      ni <- length(y_i)

      phi0 <- phi_current[[i]]
      mean_i <- mean_prop[[i]]

      d <- length(phi0)
      chain <- matrix(NA_real_, d, niter_mh + 1)
      chain[, 1] <- phi0

      for (r in seq_len(niter_mh)) {
        phi_old <- chain[, r]

        phi_prop <- if (kernel == "pop") {
          rmvnorm(mean_i, var_prop_mat)
        } else {
          rmvnorm(phi_old, kappa * var_prop_mat)
        }

        logratio <- 0

        for (j in seq_len(ni)) {
          mean_new <- g_scalar(t_i[j], phi_prop)
          mean_old <- g_scalar(t_i[j], phi_old)
          logratio <- logratio +
            dnorm(y_i[j], mean_new, sd, log = TRUE) -
            dnorm(y_i[j], mean_old, sd, log = TRUE)
        }

        if (kernel == "random_walk") {
          logratio <- logratio +
            logdmvnorm(phi_prop, mean_i, var_prop_mat) -
            logdmvnorm(phi_old, mean_i, var_prop_mat)
        }

        if (log(runif(1)) <= logratio) {
          chain[, r + 1] <- phi_prop
        } else {
          chain[, r + 1] <- phi_old
        }
      }

      phi_chains[[i]] <- chain[, niter_mh + 1]
    }

    phi_chains
  }

  list(
    g_scalar = g_scalar,
    g_vector = g_vector,
    rmvnorm = rmvnorm,
    rmvnorm_mat = rmvnorm_mat,
    logdmvnorm = logdmvnorm,
    metropolis_vector = metropolis_vector,
    ll = loglik
  )
}


#' Build C++ Backend for Model functions
#'
#' Generates and compiles a C++ backend for efficient computation of a
#' model provided by the user. This function transpiles an R
#' function to C++ and creates compiled functions for model evaluation,
#' sampling, and likelihood computation.
#'
#' @param g_fun A function that defines the structural model.
#'   This function should take parameters (phi) and time (t) as inputs.
#'
#' @param debug Logical. If TRUE, returns the generated C++ code as a string
#'   without compilation. Useful for inspecting generated code. Default is
#'  FALSE.
#'
#' @return A list containing five compiled C++ functions:
#'   \describe{
#'     \item{g_scalar}{Evaluates the structural model at a single time point.
#'       Takes phi (parameters vector) and t (time point).}
#'     \item{g_vector}{Vectorized version of g_scalar for multiple time points.
#'       Takes phi (parameters vector) and t (times vector).}
#'     \item{rmvnorm}{Generates random samples from a multivariate normal
#'       distribution. Takes mean (vector) and sigma (covariance matrix).}
#'     \item{rmvnorm_mat(mean, sigma, n)}{Generates `n` draws from a
#'  multivariate normal distribution. Returns an `n x d` matrix.}
#'     \item{logdmvnorm}{Computes log-probability density of multivariate
#'  normal.
#'       Takes x (observation vector), mean (vector), and sigma (covariance
#'  matrix).}
#'     \item{metropolis_vector}{Implements Metropolis-Hastings MCMC sampler
#'       for population parameters. Takes y (observations), t (times),
#'       phi_current (current parameters), mean_prop (proposal mean),
#'       var_prop_mat (proposal covariance), sigma2 (residual variance),
#'       niter_mh (number of iterations), kappa (scaling factor), and
#'       kernel (proposal type: "pop" or "random_walk").}
#'    \item{ll(yi_list, ti_list, phi_samples, sigma2)}{Computes the
#'  log-likelihood of observed data given parameter samples.}
#' }
#'
#'
#' @details
#' When debug = FALSE, this function:
#' \enumerate{
#'   \item Transpiles the input R function to C++
#'   \item Embeds it in a C++ template with supporting functions
#'   \item Writes the code to a temporary file
#'   \item Compiles using Rcpp::sourceCpp
#'   \item Returns compiled functions in the global environment
#' }
#'
#' @examples
#' \dontrun{
#' # Define a simple PK model
#' pk_model <- function(phi, t) {
#'   phi[1] * exp(-phi[2] * t)
#' }
#'
#' # Build the C++ backend
#' backend <- build_backend_cpp(pk_model)
#'
#' # Use compiled functions for inference
#' samples <- backend$metropolis_vector(
#'   y, t, phi_init, mean_prop, var_prop, sigma2, niter
#' )
#' }
#'
#' @keywords internal
build_backend_cpp <- function(g_fun, debug = FALSE) {
  cpp_g_function <- transpile_to_cpp(g_fun, "g_scalar_cpp")

  cpp_code <- sprintf('
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

%s

// [[Rcpp::export]]
arma::vec g_vector_cpp(const arma::vec& t, const arma::vec& phi) {
  arma::vec out(t.n_elem);
  for (int i = 0; i < t.n_elem; i++)
    out[i] = g_scalar_cpp(t[i], phi);
  return out;
}

// [[Rcpp::export]]
arma::vec rmvnorm_cpp(const arma::vec& mean, const arma::mat& sigma) {
  return mean + arma::chol(sigma) * arma::randn(mean.n_elem);
}

// [[Rcpp::export]]
arma::mat rmvnorm_mat_cpp(const arma::vec& mean, const arma::mat& sigma, int n) {
  int d = mean.n_elem;
  arma::mat out(n, d);
  arma::mat L = arma::chol(sigma); 
  for (int i = 0; i < n; i++) {
    out.row(i) = (mean + L * arma::randn(d)).t();
  }
return out;
}

// [[Rcpp::export]]
double loglik_cpp(List yi_list, List ti_list, List phi_samples_list, double sigma2) {
  int n = yi_list.size();
  double loglike = 0.0;
  for(int i = 0; i < n; i++) {
    NumericVector yi = yi_list[i];
    NumericVector ti = ti_list[i];
    int ni = yi.size();
    arma::mat phi_samples = phi_samples_list[i];
    int n_samples = phi_samples.n_rows;
    double contribution = 0.0;
    for(int s = 0; s < n_samples; s++) {
      double sumsq = 0.0;
        for(int j = 0; j < ni; j++) {
          arma::vec phi_i = phi_samples.row(s).t().eval();
          double g = g_scalar_cpp(ti[j], phi_i);
          sumsq += pow(yi[j] - g, 2);
        } 
      contribution += exp(-sumsq / (2.0 * sigma2));
    }
    loglike += log(pow(2.0 * M_PI * sigma2, -ni/2.0) * contribution / n_samples);
  }
  return loglike;
}

// [[Rcpp::export]]
double logdmvnorm_cpp(const arma::vec& x,
                      const arma::vec& mean,
                      const arma::mat& sigma) {
  arma::vec d = x - mean;
  return -0.5 * (x.n_elem * log(2 * M_PI) +
                 log(arma::det(sigma)) +
                 arma::as_scalar(d.t() * arma::inv(sigma) * d));
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
        double mean_new = g_scalar_cpp(t_i[j], phi_prop);
        double mean_old = g_scalar_cpp(t_i[j], phi_old);
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
', cpp_g_function)

  if (debug) {
    return(cpp_code)
  }

  # Create temporary file for C++ code
  tmp_cpp <- tempfile(fileext = ".cpp")
  writeLines(cpp_code, tmp_cpp)

  # Compile from C++ file and force globalenv
  Rcpp::sourceCpp(tmp_cpp, env = .GlobalEnv)

  # retourner la liste backend
  list(
    g_scalar = g_scalar_cpp,
    g_vector = g_vector_cpp,
    rmvnorm = rmvnorm_cpp,
    rmvnorm_mat = rmvnorm_mat_cpp,
    logdmvnorm = logdmvnorm_cpp,
    metropolis_vector = metropolis_vector_cpp,
    ll = loglik_cpp
  )
}


#' Compile Model with Optional C++ Backend
#'
#' Compiles a model by attempting to build a C++ backend if requested,
#' with automatic fallback to a pure R implementation.
#'
#' @param g_fun A function representing the model to be compiled.
#' @param use_cpp Logical. If TRUE, attempts to compile using C++ backend;
#'   if FALSE, uses pure R backend directly.
#' @param debug Logical. If TRUE, enables debug mode for C++ compilation.
#'   Default is FALSE.
#' @param silent Logical. If TRUE, suppresses informational messages.
#'   Default is FALSE.
#'
#' @return A compiled backend object. Returns a C++ backend if compilation
#'   succeeds and use_cpp is TRUE, otherwise returns a pure R backend.
#'
#' @details
#' The function attempts C++ compilation only if use_cpp is TRUE.
#' If C++ compilation fails, it automatically falls back to the R backend.
#' Messages are displayed during the process unless silent is TRUE.
#'
#' @examples
#' \dontrun{
#' backend <- compile_model(my_model_func, use_cpp = TRUE)
#' backend_r <- compile_model(my_model_func, use_cpp = FALSE, silent = TRUE)
#' }
#'
#' @seealso
#' \code{\link{build_backend_cpp}} for C++ backend compilation.
#' \code{\link{build_backend_r}} for R backend compilation.
#'
#' @keywords internal
compile_model <- function(g_fun, use_cpp, debug = FALSE, silent = FALSE) {
  if (use_cpp) {
    backend <- tryCatch(
      build_backend_cpp(g_fun, debug),
      error = function(e) NULL
    )

    if (!is.null(backend)) {
      if (!silent) {
        message("C++ backend successfully loaded.")
      }

      return(backend)
    } else {
      if (!silent) {
        message("C++ compilation failed. Falling back to pure R backend.")
      }
      backend <- build_backend_r(g_fun)
      return(backend)
    }
  } else {
    if (!silent) {
      message("Using pure R backend (user choice).")
    }
    backend <- build_backend_r(g_fun)
    return(backend)
  }
}
