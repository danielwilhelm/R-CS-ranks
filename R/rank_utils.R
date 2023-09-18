#' Compute ranks
#' 
#' Compute ranks with flexible handling of ties.
#'
#' @param x vector of values to be ranked
#' @param omega numeric value in [0,1], defining how ties in \code{x} (if any) are handled; default is \code{0}. See Details.
#' @param increasing logical; if \code{FALSE} (default), then large elements in \code{x} receive a small rank. Otherwise, large elements in \code{x} receive a large rank.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x}. Default: \code{FALSE}.
#' 
#' @details 
#' This function implements all possible definitions of ranks of the values in \code{x}. Different definitions of the ranks are chosen through combinations of the two arguments
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
#' 
#' @return Integer vector of the same length as \code{x} containing the ranks.
#' @examples
#' # simple example without ties:
#' x <- c(3,8,-4,10,2)
#' irank(x, increasing=TRUE)
#' irank(x, increasing=FALSE)
#' 
#' # since there are no ties, the value of omega has no impact:
#' irank(x, increasing=TRUE, omega=0)
#' irank(x, increasing=TRUE, omega=0.5)
#' irank(x, increasing=TRUE, omega=1)
#' 
#' # simple example with ties:
#' x <- c(3,4,7,7,10,11,15,15,15,15)
#' irank(x, increasing=TRUE, omega=0) # smallest possible ranks
#' irank(x, increasing=TRUE, omega=0.5) # mid-ranks
#' irank(x, increasing=TRUE, omega=1) # largest possible ranks
#' @export
irank <- function(x, omega=0, increasing=FALSE, na.rm=FALSE) {
  irank_against(x, x, omega=omega, increasing=increasing, na.rm=na.rm)
}

#' Compute integer ranks in another reference vector
#' 
#' The method \code{\link{irank}} compares ranks using the same vector as reference.
#' This method returns ranks, that values from \code{x} would assume if (individually)
#' inserted into \code{v}. 
#' 
#' @param x numeric query vector.
#' @param v numeric reference vector.
#' @inheritParams irank
#' 
#' @inherit irank details return
#' @examples 
#' irank_against(1:10, c(4,4,4,3,1,10,7,7))
#' @export
irank_against <- function(x, v, omega=0, increasing=FALSE, na.rm=FALSE){
  l <- process_irank_against_args(x=x, v=v, omega=omega, increasing=increasing, na.rm=na.rm)
  x <- l$x; v <- l$v
  n_lequal_lesser <- count_lequal_lesser(x, v)
  out <- omega * n_lequal_lesser$n_lequal + (1-omega) * n_lequal_lesser$n_lesser + 
    1 - omega
  names(out) <- names(x)
  out
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
count_lequal_lesser <- function(x, v=NULL, return_inverse_ranking=FALSE){
  if(is.null(v))
    v <- x
  assert_has_no_NAs(v, "v")
  ranking <- order(v)
  n_lower_or_equal <- findInterval(x, v[ranking], left.open = FALSE)
  n_lower <- findInterval(x, v[ranking], left.open = TRUE)
  out <- list(n_lequal=stats::setNames(n_lower_or_equal, NULL), 
       n_lesser=stats::setNames(n_lower, NULL))
  return(out)
}

#' @describeIn irank Compute fractional ranks
#' 
#' This function takes the ranking returned by \code{irank} and divides the result by \code{length(x)}. The result is a ranking with 
#' ranks in the interval [0,1]. An important special case occurs for \code{increasing=TRUE} and \code{omega=1}: in this case, the rank 
#' of the value \code{x[j]} is equal to the empirical cdf of \code{x} evaluated at \code{x[j]}.
#' 
#' @examples
#' 
#' # simple example of fractional ranks without ties:
#' x <- c(3,8,-4,10,2)
#' frank(x, increasing=TRUE)
#' frank(x, increasing=FALSE)
#' @export
frank <- function(x, omega=0, increasing=FALSE, na.rm=FALSE) 
  return(frank_against(x, x, omega, increasing, na.rm))

#' @describeIn irank_against Compute integer ranks in another reference vector
#' @export
frank_against <- function(x, v, omega=0, increasing=FALSE, na.rm=FALSE){
  if(na.rm){
    l <- sum(!is.na(v))
  } else
    l <- length(v)
  out <- irank_against(x, v, omega, increasing, na.rm)
  return(out / l)
} 
