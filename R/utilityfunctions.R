#' Compute ranks from feature values
#' 
#' Given estimates of a certain feature for a set of populations,
#' calculate the integer ranks of populations, i.e. places in ranking done by feature
#' values. The larger feature value, the higher the place and the lower the integer
#' rank (lowest, 1, is the best place).
#'
#' @param x vector of values to be ranked
#' @param omega numeric; numeric value in [0,1], each corresponding to a different definition of the rank; default is \code{0}. See Details.
#' @param increasing logical; if \code{TRUE}, then large elements in \code{x} receive a large rank. 
#' In other words, larger values in \code{x} are lower in ranking. Otherwise, large elements receive small ranks. 
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} (if any). 
#' 
#' @details 
#' \code{omega} value determines, how equal entries in x should be ranked; 
#' in other words how to handle ex aequo cases. If there are none, then the parameter 
#' does not affect the output of this function. 
#' For example, let's say, that two largest entries in x are equal.
#' Those entries could both receive rank 1 or rank 2 or some value in between.
#'
#' Suppose, that we want to assign rank to \eqn{n} largest values in an array and they are all equal.
#' Then the assigned rank is \eqn{1 + (n - 1)\omega}. 
#' Ranking n equal values, that are not largest, works similarly. See also examples below.
#' 
#' 
#' @return vector of the same dimension as \code{x} containing the ranks
#' @examples
#' irank(c(4,3,1,10,7))
#' irank(c(4,3,1,10,7), omega=1) # equal to previous ranks because there are no ties
#' irank(c(4,3,1,10,7), omega=0.5) # equal to previous ranks because there are no ties
#' irank(c(4,4,4,3,1,10,7,7))
#' irank(c(4,4,4,3,1,10,7,7), omega=1)
#' irank(c(4,4,4,3,1,10,7,7), omega=0.5) 
#' @section Details:
#' tba
#' @export
irank <- function(x, omega=0, increasing=FALSE, na.rm=FALSE) {
	if (na.rm) x <- x[!is.na(x)]
	if (increasing) {
	  ranking <- order(x)
	  sorted <- x[ranking]
	  equal_to_next <- c(diff(sorted) == 0, FALSE)
	  # block is a sequence of equal values
	  block_ends <- which(!equal_to_next)
	  block_sizes <- diff(c(0, block_ends))
	  
	  # how many populations are higher or equal in ranking?
	  # Connected to Nminus, Nplus
	  n_higher_or_equal <- rep(block_ends, times = block_sizes)
	  n_equal <- rep(block_sizes, times = block_sizes)
	  n_higher <- n_higher_or_equal - n_equal
	  minimum_rank <- n_higher + 1
	  maximum_rank <- n_higher_or_equal
	  corrected_ranking <- minimum_rank * (1-omega) + maximum_rank * omega
	  # return in order of original x
	  corrected_ranking[order(ranking)]
	} else {
	  return(irank(-x, omega = omega, increasing = TRUE))
	}
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