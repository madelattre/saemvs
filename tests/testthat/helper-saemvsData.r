# Helper to build valid arguments for saemvsData
valid_args <- list(
  y = list(rnorm(5), rnorm(3)),
  t = list(1:5, 1:3),
  x_candidates = matrix(rnorm(2 * 2), nrow = 2),
  x_forced = matrix(rnorm(2 * 1), nrow = 2)
)

# Collection of invalid cases for saemvsData
invalid_cases <- list(
  "y and t length mismatch" = list(
    args = list(y = list(rnorm(5), rnorm(3)), t = list(1:5)),
    msg = "must have the same number of elements"
  ),
  "non-numeric y" = list(
    args = list(y = list(letters[1:5], rnorm(3)), t = list(1:5, 1:3)),
    msg = "must be numeric"
  ),
  "wrong rows in x_candidates" = list(
    args = list(
      y = list(rnorm(5), rnorm(3)),
      t = list(1:5, 1:3),
      x_candidates = matrix(rnorm(6), nrow = 3)
    ),
    msg = "must have 2 rows"
  ),
  "x_candidates non-numeric" = list(
    args = list(
      y = list(rnorm(5), rnorm(3)),
      t = list(1:5, 1:3),
      x_candidates = matrix(letters[1:4], nrow = 2)
    ),
    msg = "All entries in 'x_candidates' must be numeric"
  ),
  "identical column in x_forced" = list(
    args = list(
      y = list(rnorm(5), rnorm(3)),
      t = list(1:5, 1:3),
      x_candidates = matrix(c(1, 2, 3, 4), nrow = 2),
      x_forced = matrix(c(1, 2), nrow = 2)
    ),
    msg = "Column 1 of 'x_forced' is identical to column 1 of 'x_candidates'"
  )
)
