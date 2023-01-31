#' Compute ranks
#'
#' @param x vector of values to be ranked
#' @param best logical; if \code{TRUE} (default), the rank of the j-th element is defined as the number of other elements strictly larger than the j-th. Otherwise, the rank is defined as the number of other elements strictly smaller than the j-th.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} (if any).

#' @return vector of the same dimension as \code{x} containing the ranks
#' @examples
#' xrank(c(4, 4, 4, 3, 1, 10, 7, 7))
#' @export
xrank <- function(x, best = TRUE, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  return(best * (colSums(outer(x, x, ">")) + 1) + (1 - best) * (colSums(outer(x, x, "<")) + 1))
}

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

#' Process csranks arguments
#' @noRd

process_csranks_args <- function(x, Sigma, na.rm){
  if(!is.numeric(x))
    cli::cli_abort(c("{.var x} must be a numeric vector.",
                     "x" = "Supplied {.var x} is of {.cls {typeof(x)}} type."))
  if(!is.vector(x))
    cli::cli_abort(c("{.var x} must be a numeric vector.",
                    "x" = "{.var x} is of {.cls {class(x)}} class."))
  if(!is.numeric(Sigma))
    cli::cli_abort(c("{.var Sigma} must be a numeric matrix.",
                   "x" = "{.var Sigma} is of {.cls {typeof(Sigma)}} type."))
  if(!is.matrix(Sigma))
    cli::cli_abort(c("`Sigma` must be a numeric matrix.",
                     "x" = "{.var Sigma} if of {.cls {class(Sigma)}} class.",
                     "i" = "Did you provide only a vector of variances? If so, use a diagonal covariance matrix."))
  # remove NAs
  if (na.rm) {
    ind <- !is.na(x) & apply(Sigma, 1, function(v) all(!is.na(v))) & apply(Sigma, 2, function(v) all(!is.na(v)))
    x <- x[ind]
    Sigma <- Sigma[ind, ind]
  }
  if(any(is.na(x)))
    cli::cli_abort("NA values found in `x`.")
  if(any(is.na(Sigma)))
    cli::cli_abort("NA values found in `Sigma`.")
  list(x=x,Sigma=Sigma)
}

process_indices_argument <- function(indices,p){
  if (any(is.na(indices))) indices <- 1:p
  else{
    if(!is.numeric(indices))
      cli::cli_abort(c("{.var indices} must be an integer vector.",
                       "x" = "{.var indices} is of {.cls {typeof(x)}} type."))
    deviation_from_int <- which.max(abs(indices - round(indices)))
    if(indices[deviation_from_int] < 1e-10)
      cli::cli_abort(c("{.var indices} must be an integer vector.",
                       "x" = "{.var indices[{deviation_from_int}]}={indices[deviation_from_int]} is not an integer."))
    if(!is.vector(indices))
      cli::cli_abort(c("{.var indices} must be an integer vector.",
                       "x" = "{.var indices} is of {.cls {class(x)}} class."))
    if(any(indices < 1) || any(indices > p)){
      wrong_index <- which(indices < 1 | indices > p)[1]
      cli::cli_abort(c("{.var indices} must be a vector of integer indices of x.",
                       "x" = "{.var indices[{wrong_index}]}={indices[wrong_index]} is not an integer."))
    }
  }
  indices
}

#' Indices utils
#'
#' Elements of `matrix` can be accessed by double indices `M[i,j]`
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