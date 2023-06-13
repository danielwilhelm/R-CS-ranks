#' partition vector into quantile bins
#'
#' @param x vector of values to be partitioned
#' @param n number of bins

#' @return vector of the same dimension as \code{x} containing the a bin membership indicator
#' @noRd
createbins <- function(x, n) {
  bins <- cut(x, breaks = quantile(x, probs = seq(0, 1, by = 1 / n), na.rm = TRUE), include.lowest = TRUE)
  levels(bins) <- as.character(1:n)
  return(bins)
}

#' partition vector into quartile bins
#'
#' @param x vector of values to be partitioned

#' @return vector of the same dimension as \code{x} containing the a quartile bin membership indicator
#' @noRd
createquartiles <- function(x) {
  return(createbins(x, 4))
}

#' Indices utils
#'
#' Elements of `matrix` can be accessed by double indices `M[i,j]` or single `M[k]`.
#' This function allows to switch from the latter kind of indices to the former.
#'
#' @noRd
get_double_from_single_indices <- function(indices, matrix_size) {
  row_indices <- indices %% matrix_size
  row_indices[row_indices == 0] <- matrix_size
  matrix(c(
    row_indices,
    ceiling(indices / matrix_size)
  ), ncol = 2)
}