#' Compute ranks from feature values
#' 
#' Given estimates of a certain feature for a set of populations,
#' calculate the integer ranks of populations, i.e. places in ranking done by feature
#' values. The larger (or smaller) feature value, the higher the place and the lower the integer
#' rank (lowest, 1, is the best place).
#'
#' @param x vector of values to be ranked
#' @param omega numeric; numeric value in [0,1], each corresponding to a different definition of the rank; default is \code{0}. See Details.
#' @param increasing logical; if \code{TRUE}, then large elements in \code{x} receive a large rank. 
#' In other words, larger values in \code{x} are lower in ranking. Otherwise, large elements receive small ranks. 
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x}.
#' In other case the output for NAs is NA and for other values it's for extreme
#' possibilities that NA values are actually in first or last positions of ranking.
#' 
#' @details 
#' \code{omega} (\eqn{\omega}) value determines, how equal entries in \code{x} should be ranked; 
#' in other words how to handle ex aequo cases. If there are none, then the parameter 
#' does not affect the output of this function. 
#' For example, let's say, that \eqn{n} largest entries in \code{x} are equal.
#' Those entries could receive (minimum) rank 1 or (maximum) rank \eqn{n} or some value in between.
#'
#' Suppose, that we want to assign rank to \eqn{n} equal values in an array.
#' Denote their minimum rank as \eqn{r} and maximum as \eqn{R = r + n - 1}.
#' Then the assigned rank is an average of 
#' minimum and maximum rank, weighted by \eqn{\omega}: 
#' \deqn{r(1-\omega) + R\omega} 
#' 
#' @return Integer vector of the same length as \code{x} containing the ranks.
#' @examples
#' irank(c(4,3,1,10,7))
#' irank(c(4,3,1,10,7), omega=1) # equal to previous ranks because there are no ties
#' irank(c(4,3,1,10,7), omega=0.5) # equal to previous ranks because there are no ties
#' irank(c(4,4,4,3,1,10,7,7))
#' irank(c(4,4,4,3,1,10,7,7), omega=1)
#' irank(c(4,4,4,3,1,10,7,7), omega=0.5) 
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
  l <- process_compare_args(x=x, v=v, omega=omega, increasing=increasing, na.rm=na.rm)
  x <- l$x; v <- l$v
  irank_min_max <- irank_minmax(x, v)
  out <- omega * irank_min_max[,1] + (1-omega) * irank_min_max[,2]  + 1 - omega
  names(out) <- names(x)
  out
}

#' Compute minimum and maximum integer ranks in another reference vector
#' 
#' @return A matrix of size length(x), 2. In column there are integer ranks in 
#' extreme cases of omega=0 and omega=1; increasing=TRUE.
#' 
#' @noRd
irank_minmax <- function(x, v){
  ranking <- order(v)
  n_lower_or_equal <- findInterval(x, v[ranking], left.open = FALSE)
  n_lower <- findInterval(x, v[ranking], left.open = TRUE)
  return(cbind(n_lower_or_equal, n_lower, deparse.level = 0))
}

#' @describeIn irank Compute fractional ranks
#' 
#' This method returns ranks in form of fractions from [0-1] interval.
#' Smaller values (closer to 0) indicate higher rank.
#' 
#' @examples
#' frank(c(4,3,1,10,7))
#' frank(c(4,3,1,10,7), omega=1) # equal to previous ranks because there are no ties
#' frank(c(4,3,1,10,7), omega=0.5) # mid-ranks, equal to previous ranks because there are no ties
#' frank(c(4,4,4,3,1,10,7,7))
#' frank(c(4,4,4,3,1,10,7,7), omega=1)
#' frank(c(4,4,4,3,1,10,7,7), omega=0.5) # mid-ranks
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

#' Comparator function
#' 
#' @inheritParams irank_against
#' @return Matrix M of size `length(x)`, `length(v)`. For increasing = FALSE M[i,j] = 
#' 1 if x[i] < v[j]
#' 0 if x[i] > v[j]
#' omega if x[i] == v[j]
#' For decreasing = FALSE, the `<` and `>` in above definition is swapped.
#' 
#' @noRd
compare <- function(x, v=NULL, omega=0, increasing=FALSE, na.rm=FALSE){
  l <- process_compare_args(x, v, omega, increasing, na.rm)
  x <- l$x; v <- l$v
  irank_min_max <- irank_minmax(x, v)
  ranking <- order(v)
  n_higher_or_equal <- irank_min_max[,1]
  n_higher <- irank_min_max[,2]
  n_equal <- n_higher_or_equal - n_higher
  n_lower <- length(v) - n_higher - n_equal
  
  out_for_sorted_v <- matrix(rep(
      rep(c(1, omega, 0), times = length(x)),
      times = as.vector(rbind(n_higher, n_equal, n_lower))
  ), byrow = TRUE, nrow = length(x))
  
  # return in order of original v
  original_order <- order(ranking)
  out <- out_for_sorted_v[,original_order]
  out
}