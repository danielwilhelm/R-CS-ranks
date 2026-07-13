#' Compute ranks
#'
#' Compute integer of fractional ranks with flexible handling of ties.
#'
#' @param x vector of values to be ranked
#' @param omega numeric value in \[0,1\], defining how ties in \code{x} (if any) are handled; default is \code{0}. See Details.
#' @param increasing logical; if \code{FALSE} (default), then large elements in \code{x} receive a small rank. Otherwise, large elements in \code{x} receive a large rank.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x}. Default: \code{FALSE}.
#' @param weights an optional vector of weights to be used in the ranking process. Should be ‘NULL’ or a numeric vector.
#'
#' @details
#' `irank` implements all possible definitions of ranks of the values in \code{x}. Different definitions of the ranks are chosen through combinations of the two arguments
#' \code{omega} and \code{increasing}. Suppose \code{x} is of length \eqn{p}. If \code{increasing=TRUE}, then the largest value in \code{x} receives the rank \eqn{p} and the smallest
#' the rank \eqn{1}. If \code{increasing=FALSE}, then the largest value in \code{x} receives the rank \eqn{1} and the smallest
#' the rank \eqn{p}.
#'
#' The value of \code{omega} indicates how ties are handled. If there are no ties in \code{x}, then the value of \code{omega} does not affect the ranks and the only choice to be made is whether
#' the ranks should be increasing or decreasing with the values in \code{x}. When there are ties in \code{x}, however, then there are infinitely
#' many possible ranks that can be assigned to a tied value.
#'
#' When \code{increasing=TRUE}, then \code{omega=0} leads to the smallest possible and \code{omega=1} to the largest possible rank of a tied value. Values of \code{omega} between
#' 0 and 1 lead to values of the rank between the largest and smallest.
#'
#' Gotcha: the definition of weighted ranks drops a constant (1-`omega`). This means that `irank` et al. called with equal weights will have different outcome from `irank` et al. called with `weights=NULL`.
#' The values will be shifted by a constant 1-`omega`.
#'
#'
#' @return Numeric vector of the same length as \code{x} containing the integer (for `irank`) or fractional (for `frank`) ranks.
#' @examples
#' # simple example without ties:
#' x <- c(3, 8, -4, 10, 2)
#' irank(x, increasing = TRUE)
#' irank(x, increasing = FALSE)
#'
#' # since there are no ties, the value of omega has no impact:
#' irank(x, increasing = TRUE, omega = 0)
#' irank(x, increasing = TRUE, omega = 0.5)
#' irank(x, increasing = TRUE, omega = 1)
#'
#' # simple example with ties:
#' x <- c(3, 4, 7, 7, 10, 11, 15, 15, 15, 15)
#' irank(x, increasing = TRUE, omega = 0) # smallest possible ranks
#' irank(x, increasing = TRUE, omega = 0.5) # mid-ranks
#' irank(x, increasing = TRUE, omega = 1) # largest possible ranks
#'
#' @export
irank <- function(x, omega = 0, increasing = FALSE, na.rm = FALSE, weights = NULL) {
  irank_against(x, x, omega = omega, increasing = increasing, na.rm = na.rm, weights = weights)
}

#' Compute integer ranks in another reference vector
#'
#' The method \code{\link{irank}} compares ranks using the same vector as reference.
#' `irank_against` returns integer ranks, that values from \code{x} would assume if (individually)
#' inserted into \code{v}. `frank_against` acts analogously, returning fractional ranks.
#'
#' @param x numeric query vector.
#' @param v numeric reference vector.
#' @inheritParams irank
#'
#' @details
#' It's useful to think about `frank_against(x,v)` as a generalization of Empirical Cumulative
#' Distribution Function, created for `v` and evaluated for points in `x`.
#' `frank_agaist(x,v,increasing=TRUE,omega=1)` is identical
#' to `ecdf(v)(x)`.
#'
#' `increasing` switches the inequality sign in ECDF definition from
#' \eqn{F_V(t) = \hat P(V <= t)} to \eqn{\hat P(V >= t)}.
#'
#' `omega=0` introduces the strict inequality (\eqn{\hat P(V < t)} instead of \eqn{\hat P(V <= t)}).
#' Any `omega` in between is a weighted average of the cases `omega=1` and `omega=0`.
#'
#' Finally, `irank_against` is equal to `frank_against` multiplied by the `length(v)`.
#'
#' This particular choice of default parameters was made for compatibility with default parameters of
#' `irank` and `frank`. `irank(x)` is always equal to `irank_against(x,x)` and `frank(x)` is always equal to `frank_against(x,x)`.
#'
#' @return Numeric vector of the same length as \code{x} containing the integer (for `irank_against`) or fractional (for `frank_against`) ranks.
#' @examples
#' irank_against(1:10, c(4, 4, 4, 3, 1, 10, 7, 7))
#' @seealso [irank()], [ecdf()]
#' @export
irank_against <- function(x, v, omega = 0, increasing = FALSE, na.rm = FALSE, weights = NULL) {
  l <- process_irank_against_args(x = x, v = v, omega = omega, increasing = increasing, na.rm = na.rm, weights = weights)
  x <- l$x
  v <- l$v
  weights <- l$weights

  v_order <- order(v)
  v_ordered <- v[v_order]

  if (!is.null(weights)) {
    weights_ordered <- weights[v_order]
    weights_cumsum <- cumsum(weights_ordered)
  } else {
    weights_cumsum <- NULL
  }

  n_lequal_lesser <- count_lequal_lesser(x, v_ordered)
  n_lequal_weighted <- get_weighted_counts(n_lequal_lesser$n_lequal, weights_cumsum)
  n_lesser_weighted <- get_weighted_counts(n_lequal_lesser$n_lesser, weights_cumsum)
  out <- omega * n_lequal_weighted + (1 - omega) * n_lesser_weighted

  # Backward compatibility with rank definition without weights
  if (is.null(weights)) {
    out <- out + 1 - omega
  }

  names(out) <- names(x)
  out
}

#' @param counts vector of integers, referring to indices of weights. Could contain 0s.
#' @param weights vector of non-negative numbers. Could be NULL, in which case `counts` is returned as-is.
#' @noRd
get_weighted_counts <- function(counts, weights = NULL) {
  if (is.null(weights)) {
    return(counts)
  }
  output <- numeric(length(counts)) # Init with 0s
  counts_in_range <- counts > 0
  output[counts_in_range] <- weights[counts[counts_in_range]]
  output
}

#' Compute minimum and maximum integer ranks in another reference vector
#'
#' For each element of query vector x:
#'     count, how many observations in the reference vector v are lesser
#'     (returned in n_lesser element)
#'     and lower or equal (returned in n_lequal element) than this element.
#'
#' @param v If NULL - set it as x. An often usecase.
#' @return A list with 2 (or 3 in case of return_inverse_ranking) elements.
#'
#' @noRd
count_lequal_lesser <- function(x, v_ordered = NULL, return_inverse_ranking = FALSE) {
  if (is.null(v_ordered)) {
    v_ordered <- x
  }
  assert_has_no_NAs(v_ordered, "v")
  n_lower_or_equal <- findInterval(x, v_ordered, left.open = FALSE)
  n_lower <- findInterval(x, v_ordered, left.open = TRUE)
  out <- list(
    n_lequal = stats::setNames(n_lower_or_equal, NULL),
    n_lesser = stats::setNames(n_lower, NULL)
  )
  return(out)
}

#' @rdname irank
#' @details
#' `frank` takes the ranking returned by \code{irank} and divides the result by \code{length(x)}. The result is a ranking with
#' ranks in the interval \[0,1\]. An important special case occurs for \code{increasing=TRUE} and \code{omega=1}: in this case, the rank
#' of the value \code{x[j]} is equal to the empirical cdf of \code{x} evaluated at \code{x[j]}.
#'
#' @examples
#' # simple example of fractional ranks without ties:
#' x <- c(3, 8, -4, 10, 2)
#' frank(x, increasing = TRUE)
#' frank(x, increasing = FALSE)
#' @export
frank <- function(x, omega = 0, increasing = FALSE, na.rm = FALSE, weights = NULL) {
  return(frank_against(x, x, omega, increasing, na.rm, weights = weights))
}

#' @rdname irank_against
#' @export
frank_against <- function(x, v, omega = 0, increasing = FALSE, na.rm = FALSE, weights = NULL) {
  l <- process_irank_against_args(x = x, v = v, omega = omega, increasing = increasing, na.rm = na.rm, weights = weights)
  divisor <- get_maximum_rank(l$v, l$weights)
  out <- irank_against(x, v, omega, increasing, na.rm)
  return(out / divisor)
}

#' @noRd
get_maximum_rank <- function(v, weights = NULL) {
  if (is.null(weights)) {
    return(length(v))
  } else {
    return(sum(weights))
  }
}
