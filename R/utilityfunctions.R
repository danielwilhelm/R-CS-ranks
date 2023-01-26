#' Compute integer ranks
#'
#' @param x vector of values to be ranked
#' @param omega numeric; numeric value in [0,1], each corresponding to a different definition of the rank; default is \code{0}. See Details.
#' @param increasing logical; if \code{TRUE}, then large elements in \code{x} receive a large rank. Otherwise, large elements receive small ranks. 
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} (if any). 

#' @return vector of the same dimension as \code{x} containing the ranks
#' @examples
#' irank(c(4,4,4,3,1,10,7,7))
#' @export
irank <- function(x, omega=0, increasing=FALSE, na.rm=FALSE) {
	if (na.rm) x <- x[!is.na(x)]
	if (increasing) {
		return( omega*colSums(outer(x, x, "<=")) + (1-omega)*colSums(outer(x, x, "<")) + 1 - omega )
	} else {
		return( omega*colSums(outer(x, x, ">=")) + (1-omega)*colSums(outer(x, x, ">")) + 1 - omega )
	}
}


#' Compute fractional ranks
#'
#' @inheritParams csranks

#' @return vector of the same dimension as \code{x} containing the ranks
#' @examples
#' frank(c(4,4,4,3,1,10,7,7))
#' @export
frank <- function(x, omega=0, increasing=FALSE, na.rm=FALSE) return(irank(x, omega, increasing, na.rm) / length(x))



#' partition vector into quantile bins
#'
#' @param x vector of values to be partitioned
#' @param n number of bins

#' @return vector of the same dimension as \code{x} containing the a bin membership indicator
#' @noRd
createbins <- function(x, n) {
	bins <- cut(x, breaks=quantile(x, probs=seq(0,1, by=1/n), na.rm=TRUE), include.lowest=TRUE)	
	levels(bins) <- as.character(1:n)
	return(bins)
}

#' partition vector into quartile bins
#'
#' @param x vector of values to be partitioned

#' @return vector of the same dimension as \code{x} containing the a quartile bin membership indicator
#' @noRd
createquartiles <- function(x) return(createbins(x,4))
