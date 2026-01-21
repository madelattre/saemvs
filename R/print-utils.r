############################################################
# Utility functions (internal, not exported)
############################################################

#' @keywords internal
#' @noRd
.format_matrix <- function(mat, max_rows = 5, max_cols = 5, digits = 3,
                           rownames = TRUE, colnames = TRUE,
                           phi_colnames = FALSE, phi_prefix = "phi") {
  if (is.null(mat)) {
    return("NULL")
  }

  nr <- nrow(mat)
  nc <- ncol(mat)

  mat_disp <- mat
  if (nr > max_rows) mat_disp <- mat_disp[1:max_rows, , drop = FALSE]
  if (nc > max_cols) mat_disp <- mat_disp[, 1:max_cols, drop = FALSE]

  mat_disp <- as.data.frame(mat_disp)
  for (j in seq_along(mat_disp)) {
    x <- mat_disp[[j]]
    x_fmt <- formatC(x, digits = digits, format = "f")
    w <- max(nchar(x_fmt), na.rm = TRUE)
    mat_disp[[j]] <- formatC(x, digits = digits, format = "f", width = w)
  }

  if (colnames) {
    cn <- colnames(mat_disp)

    if (phi_colnames) {
      default_cn <- paste0("V", seq_len(ncol(mat_disp)))

      if (is.null(cn) || identical(cn, default_cn)) {
        cn <- paste0(phi_prefix, seq_len(ncol(mat_disp)))
      } else {
        cn <- paste0(phi_prefix, "_", cn)
      }
    }

    colnames(mat_disp) <- cn
  }
  mat_disp <- as.matrix(mat_disp)

  if (!rownames) rownames(mat_disp) <- NULL
  if (!colnames) colnames(mat_disp) <- NULL

  txt <- capture.output(print(mat_disp, quote = FALSE, right = TRUE))

  if (nr > max_rows) txt <- c(txt, "...")

  header <- paste0("(", nr, " rows x ", nc, " columns)")
  paste(c(header, txt), collapse = "\n")
}

#' @keywords internal
#' @noRd
.format_vector <- function(
    x, max_length = 5, digits = 3, indent = "",
    full = FALSE) {
  n <- length(x)
  if (full || n <= max_length) {
    x_disp <- x
  } else {
    x_disp <- x[1:max_length]
  }

  x_fmt <- formatC(x_disp, digits = digits, format = "f")
  w <- max(nchar(x_fmt), na.rm = TRUE)
  x_fmt <- formatC(x_disp, digits = digits, format = "f", width = w)

  txt <- paste(x_fmt, collapse = ", ")
  if (!full && n > max_length) {
    txt <- paste0(txt, ", ... [length=", n, "]")
  }

  txt <- gsub("\n", paste0("\n", indent), txt)
  paste0("c(", txt, ")")
}

#' @keywords internal
#' @noRd
.format_list_of_vectors <- function(lst,
                                   max_items = 3,
                                   max_vec_length = 5,
                                   digits = 3,
                                   indent = "  ") {
  if (length(lst) == 0) {
    return("list()")
  }

  out <- vapply(
    seq_len(min(length(lst), max_items)),
    function(i) {
      vec_str <- .format_vector(
        lst[[i]],
        max_length = max_vec_length,
        digits = digits,
        indent = paste0(indent, "  ")
      )
      paste0(indent, "[[", i, "]] ", vec_str)
    },
    character(1)
  )

  if (length(lst) > max_items) {
    out <- c(
      out,
      paste0(indent, "... (", length(lst) - max_items, " more elements)")
    )
  }

  paste(out, collapse = "\n")
}

