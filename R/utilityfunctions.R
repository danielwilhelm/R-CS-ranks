#' Compute ranks
#'
#' @param x vector of values to be ranked
#' @param best logical; if \code{TRUE} (default), the rank of the j-th element is defined as the number of other elements strictly larger than the j-th. Otherwise, the rank is defined as the number of other elements strictly smaller than the j-th.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} (if any). 

#' @return vector of the same dimension as \code{x} containing the ranks
#' @examples
#' xrank(c(4,4,4,3,1,10,7,7))
#' @export
xrank <- function(x, best=TRUE, na.rm=FALSE) {
	if (na.rm) x <- x[!is.na(x)]
	return(best*(colSums(outer(x,x,'>'))+1) + (1-best)*(colSums(outer(x,x,'<'))+1))
}

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
