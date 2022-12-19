
#' Confidence sets for ranks based on multinomial data
#'
#' Given data with counts of successes for each category, calculate confidence set for ranks by most successes for categories.
#'
#' @param x vector of counts of successes for each category
#' @param multcorr multiplicity correction to be used: \code{Bonferroni} (default), \code{Holm} or \code{none}.
#' @inheritParams csranks
#' @inherit csranks return
#' @section Details
#' Bonferroni is TODO, Holm is TODO
#' @examples
#' x <- c(rmultinom(1, 1000, 1:10))
#' csranks_multinom(x)
#' @export
csranks_multinom <- function(x, coverage=0.95, cstype="two-sided", simul=TRUE, multcorr="Bonferroni", indices=NA, na.rm=FALSE) {

	if (simul) {
		return(csranks_multinom_simul(x, coverage=coverage, cstype=cstype, multcorr=multcorr, indices=indices, na.rm=na.rm))	
	} else {
		return(csranks_multinom_marg(x, coverage=coverage, cstype=cstype, multcorr=multcorr, indices=indices, na.rm=na.rm))
	}
}



#' Simultaneous Confidence sets for ranks based on multinomial data
#'
#' This function is called by \code{\link{csranks_multinom}} when \code{simul=TRUE}.
#'
#' @param x vector of counts of successes for each category
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param cstype type of confidence set (\code{two-sided}, \code{upper}, \code{lower}). Default is \code{two-sided}.
#' @param multcorr multiplicity correction to be used: \code{Bonferroni} (default) or \code{Holm}.
#' @param indices vector of indices of \code{x} for whose ranks the confidence sets are computed. \code{indices=NA} (default) means computation for all ranks.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{sd} (if any). 

#' @return L, U lower and upper bounds of the confidence set for ranks indicated in \code{indices}.

#' @export
#' @importFrom stats rmultinom
#' @importFrom stats pbinom
#' @importFrom stats aggregate
csranks_multinom_simul <- function(x, coverage=0.95, cstype="two-sided", multcorr="Bonferroni", indices=NA, na.rm=FALSE) {

	# check arguments
	cstype <- match.arg(cstype, c("two-sided", "upper", "lower"))
	multcorr <- match.arg(multcorr, c("Bonferroni", "Holm"))

	# remove NAs
	if (na.rm) x <- x[!is.na(x)]
	stopifnot(all(!is.na(x)))

	# initializations
	S <- outer(x,x,'+')
	p <- length(x)

	# which comparisons?
	if (any(is.na(indices))) indices <- 1:p
	ind <- matrix(TRUE,p,p)
	if (cstype=="two-sided") {
		ind[-indices,] <- FALSE
		ind[,indices] <- TRUE	
	}
	if (cstype=="upper") {
		ind[,-indices] <- FALSE	
	}
	if (cstype=="lower") {
		ind[-indices,] <- FALSE
	}
	diag(ind) <- FALSE

	# how many comparisons?
	ncomp <- sum(ind)
	
	# arrange data for comparisonsdata
	Svec <- c(S[ind])
		
	indj <- matrix(rep(1:p,p),p,p)
	indjvec <- indj[ind]
	indk <- t(indj)
	indkvec <- indk[ind]

	mat <- data.frame(xj=x[indkvec], S=Svec)
	
	# compute p-values 
	# pval <- function(x) sum(choose(x[2],x[1]:x[2])*2^(-x[2]))
	pval <- function(x) 1-pbinom(x[1]-1,x[2],1/2)
	mat$pval <- apply(mat, 1, pval)

	# add pair labels
	mat$j <- indkvec
	mat$k <- indjvec	

	# order pairs by their p-values (from smallest to largest)
	mat <- mat[order(mat$pval),]

	# compute critical value
	beta <- switch(multcorr,
		"Bonferroni" = (1-coverage)/ncomp,
		"Holm" = (1-coverage)/(ncomp+1-(1:ncomp)))	
	
	# perform test(s)
	mat$rej <- mat$pval<=beta
	if (multcorr=="Holm") {
		ind_first_false <- which(!mat$rej)[1]
		if (!is.na(ind_first_false)) mat$rej[ind_first_false:length(mat$rej)] <- FALSE
	}

	# compute confidence set for ranks
	if (cstype=="upper" | cstype=="two-sided") {
		rej_plus <- aggregate(mat["rej"], by=mat["j"], sum)
		ind_plus <- rej_plus$j %in% indices
		Nplus <- rej_plus$rej[ind_plus]
	} else {
		Nplus <- 0*indices
	}

	if (cstype=="lower" | cstype=="two-sided") {
		rej_minus <- aggregate(mat["rej"], by=mat["k"], sum)
		ind_minus <- rej_minus$k %in% indices
		Nminus <- rej_minus$rej[ind_minus]
	} else {
		Nminus <- indices*0
	}

	L <- Nminus+1
	U <- p-Nplus

	return(list(L=as.integer(L), U=as.integer(U)))
}


#' Marginal Confidence sets for ranks based on multinomial data
#'
#' This function is called by \code{\link{csranks_multinom}} when \code{simul=FALSE}.
#'
#' @param x vector of counts of successes for each category
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param cstype type of confidence set (\code{two-sided}, \code{upper}, \code{lower}). Default is \code{two-sided}.
#' @param multcorr multiplicity correction to be used: \code{Bonferroni} (default) or \code{Holm}.
#' @param indices vector of indices of \code{x} for whose ranks the confidence sets are computed. \code{indices=NA} (default) means computation for all ranks.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{sd} (if any). 

#' @return L, U lower and upper bounds of the confidence set for ranks indicated in \code{indices}.

#' @export
csranks_multinom_marg <- function(x, coverage=0.95, cstype="two-sided", multcorr="Bonferroni", indices=NA, na.rm=FALSE) {

	# remove NAs
	if (na.rm) x <- x[!is.na(x)]
	stopifnot(all(!is.na(x)))

	# initializations
	p <- length(x)
	if (any(is.na(indices))) indices <- 1:p
	L <- rep(NA,length(indices)); U <- L

	# compute marginal CS for each category indicated by indices
	for (i in 1:length(indices)) {
		CS <- csranks_multinom_simul(x, coverage=coverage, cstype=cstype, multcorr=multcorr, indices=indices[i])
		L[i] <- CS$L
		U[i] <- CS$U
	}

	return(list(L=as.integer(L),U=as.integer(U)))
}



