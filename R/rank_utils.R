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
#' @return vector of the same length as \code{x} containing the ranks
#' @examples
#' irank(c(4,3,1,10,7))
#' irank(c(4,3,1,10,7), omega=1) # equal to previous ranks because there are no ties
#' irank(c(4,3,1,10,7), omega=0.5) # equal to previous ranks because there are no ties
#' irank(c(4,4,4,3,1,10,7,7))
#' irank(c(4,4,4,3,1,10,7,7), omega=1)
#' irank(c(4,4,4,3,1,10,7,7), omega=0.5) 
#' @export
irank <- function(x, omega=0, increasing=FALSE, na.rm=FALSE) {
  check_irank_args(x, omega, increasing, na.rm)
  if (na.rm) x <- x[!is.na(x)]
  compares <- compare(x, omega=omega, increasing = increasing, na.rm = na.rm)
  compares[is.na(compares)] <- omega
  out <- rowSums(compares) + 1 - omega
  out[is.na(x)] <- NA
  out
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
frank <- function(x, omega=0, increasing=FALSE, na.rm=FALSE) return(irank(x, omega, increasing, na.rm) / length(x))

#' Comparator function
#' 
#' @inheritParams irank
#' @return Matrix M of size `length(x)`, `length(x)`. For increasing = FALSE M[i,j] = 
#' 1 if v[i] < v[j]
#' 0 if v[i] > v[j]
#' omega if v[i] == v[j]
#' For decreasing = FALSE, the `<` and `>` in above definition is swapped.
#' 
#' @noRd
compare <- function(x, omega, increasing, na.rm){
  if(!increasing){
    x <- -x
  }
  if (na.rm){
    was_NA <- rep(FALSE, sum(!is.na(x)))
  } else {
    was_NA <- is.na(x)
  }
  
  x <- x[!is.na(x)]
  n_NAs <- sum(was_NA)
  
  ranking <- order(x)
  sorted <- x[ranking]
  equal_to_next <- c(diff(sorted) == 0, FALSE)
  # block is a sequence of equal values
  # In the outcome matrix for sorted values, those are blocks on diagonal
  block_ends <- which(!equal_to_next)
  block_starts <- c(1, block_ends[-length(block_ends)] + 1)
  
  out_for_sorted <- matrix(0, nrow = length(x), ncol = length(x))
  out_for_sorted[lower.tri(out_for_sorted)] <- 1
  for(i in 1:length(block_ends)){
    out_for_sorted[block_starts[i]:block_ends[i],
                   block_starts[i]:block_ends[i]] <- omega
  }
  
  # return in order of original x
  original_order <- order(ranking)
  out_without_nas <- out_for_sorted[original_order,original_order]
  
  # correct for NAs
  out <- matrix(NA, nrow = length(x) + n_NAs, ncol = length(x) + n_NAs)
  out[!was_NA, !was_NA] <- out_without_nas
  out
}