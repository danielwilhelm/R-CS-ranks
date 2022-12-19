#' Confidence sets for ranks
#' 
#' Given estimates and stds of a certain feature, compute confidence sets for ranks based on the feature values.
#'
#' @param x vector of estimates
#' @param sd vector of standard errors of \code{x}
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param cstype type of confidence set (\code{two-sided}, \code{upper}, \code{lower}). Default is \code{two-sided}.
#' @param stepdown logical; if \code{TRUE} (default), stepwise procedure is used, otherwise single step procedure is used. See Details section for more.
#' @param R number of bootstrap replications. Default is 1000.
#' @param simul logical; if \code{TRUE} (default), then simultaneous confidence sets are computed, which jointly cover all populations indicated by \code{indices}. 
#'		Otherwise, for each population indicated in \code{indices} a marginal confidence set is computed.
#' @param indices vector of indices of \code{x} for whose ranks the confidence sets are computed. \code{indices=NA} (default) means computation for all ranks.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{sd} (if any). 
#' @param seed seed for bootstrap random variable draws. If set to \code{NA} (default), then seed is not set.

#' @return L, U lower and upper bounds of the confidence set for ranks indicated in \code{indices}.

#' @examples
#' x <- seq(1,3,length=10)
#' sd <- rep(0.2,10)
#' csranks(x, sd)

#' @section Details:
#' Stepwise procedure is TODO, Single step procedure is TODO.
#' Parametric bootstrap used to calculate distribution for confidence sets based on the normal distribution with independent populations.
#'
#' @references {1:Mogstad, Romano, Shaikh, and Wilhelm ("Inference for Ranks with Applications to Mobility across Neighborhoods and Academic Achievements across Countries", \href{https://www.ucl.ac.uk/~uctpdwi/papers/cwp1020.pdf}{CeMMAP Working Paper CWP10/20})}
#' @export
csranks <- function(x, sd, coverage=0.95, cstype="two-sided", stepdown=TRUE, R=1000, simul=TRUE, indices=NA, na.rm=FALSE, seed=NA) {
	if (simul) {
		return(csranks_simul(x, sd, coverage=coverage, cstype=cstype, stepdown=stepdown, R=R, indices=indices, na.rm=na.rm, seed=seed))
	} else {
		return(csranks_marg(x, sd, coverage=coverage, cstype=cstype, stepdown=stepdown, R=R, indices=indices, na.rm=na.rm, seed=seed))
	}
}

#' Simultaneous confidence sets for ranks
#'
#' This function is called by \code{\link{csranks}} when \code{simul=TRUE}.
#'
#' @param x vector of estimates
#' @param sd vector of standard errors of \code{x}
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param cstype type of confidence set (\code{two-sided}, \code{upper}, \code{lower}). Default is \code{two-sided}.
#' @param stepdown logical; if \code{TRUE} (default), stepwise procedure is used, otherwise single step procedure is used.
#' @param R number of bootstrap replications. Default is 1000.
#' @param indices vector of indices of \code{x} for whose ranks the confidence sets are computed. \code{indices=NA} (default) means computation for all ranks.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{sd} (if any). 
#' @param seed seed for bootstrap random variable draws. If set to \code{NA} (default), then seed is not set.

#' @return L, U lower and upper bounds of the confidence set for ranks indicated in \code{indices}.

#' @section Details:
#' Implentation of the simultaneous confidence sets proposed in Mogstad, Romano, Shaikh, and Wilhelm ("Inference for Ranks with Applications to Mobility across Neighborhoods and Academic Achievements across Countries", \href{https://www.ucl.ac.uk/~uctpdwi/papers/cwp1020.pdf}{CeMMAP Working Paper CWP10/20}).
#' Parametric bootstrap based on the normal distribution with independent populations.
#' @export
csranks_simul <- function(x, sd, coverage=0.95, cstype="two-sided", stepdown=TRUE, R=1000, indices=NA, na.rm=FALSE, seed=NA) {

	# remove NAs
	if (na.rm) {
		ind <- !is.na(x) & !is.na(sd)
		x <- x[ind]
		sd <- sd[ind]
	}
	stopifnot(all(!is.na(x) & !is.na(sd)))

	# joint CS for difference in means
	csdifftype <- switch(cstype,
		"lower" = "upper",
		"upper" = "lower",
		"two-sided" = "symmetric")
	p <- length(x)
	if (any(is.na(indices))) indices <- 1:p
	res <- csdiffmeans(x, sd, coverage=coverage, indices=indices, cstype=csdifftype, stepdown=stepdown, R=R, seed=seed)
	L <- res$L; U <- res$U

	# compute Nminus and Nplus
	if (stepdown & cstype=="two-sided") {
		if (length(indices)==1 & all(!is.na(indices))) {
			Nplus <- sum(L[indices,]>0, na.rm=TRUE)
			Nminus <- sum(L[,indices]>0, na.rm=TRUE)
		} else {
			Nplus <- rowSums(L[indices,]>0, na.rm=TRUE)
			Nminus <- colSums(L[,indices]>0, na.rm=TRUE)
		}
	} else {
		if (length(indices)==1 & all(!is.na(indices))) {
			Nplus <- sum(L[indices,]>0, na.rm=TRUE)
			Nminus <- sum(U[indices,]<0, na.rm=TRUE)
		} else {
			Nplus <- rowSums(L[indices,]>0, na.rm=TRUE)
			Nminus <- rowSums(U[indices,]<0, na.rm=TRUE)
		}	
	}

	# return lower and upper confidence bounds for the ranks
	return(list(L=as.integer(Nminus+1),U=as.integer(p-Nplus)))
}

#' Marginal confidence sets for ranks
#'
#' This function is called by \code{\link{csranks}} when \code{simul=FALSE}.
#'
#'
#' @param x vector of estimates
#' @param sd vector of standard errors of \code{x}
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param cstype type of confidence set (\code{two-sided}, \code{upper}, \code{lower}). Default is \code{two-sided}.
#' @param stepdown logical; if \code{TRUE} (default), stepwise procedure is used, otherwise single step procedure is used.
#' @param R number of bootstrap replications. Default is 1000.
#' @param indices vector of indices of \code{x} for whose ranks the confidence sets are computed. \code{indices=NA} (default) means computation for all ranks.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{sd} (if any). 
#' @param seed seed for bootstrap random variable draws. If set to \code{NA} (default), then seed is not set.

#' @return L, U lower and upper bounds of the confidence set for ranks indicated in \code{indices}.

#' @section Details:
#' Implentation of the marginal confidence sets proposed in Mogstad, Romano, Shaikh, and Wilhelm ("Inference for Ranks with Applications to Mobility across Neighborhoods and Academic Achievements across Countries", \href{https://www.ucl.ac.uk/~uctpdwi/papers/cwp1020.pdf}{CeMMAP Working Paper CWP10/20}).
#' Parametric bootstrap based on the normal distribution with independent populations.
#' @export
csranks_marg <- function(x, sd, coverage=0.95, cstype="two-sided", stepdown=TRUE, R=1000, indices=NA, na.rm=FALSE, seed=NA) {

	# remove NAs
	if (na.rm) {
		ind <- !is.na(x) & !is.na(sd)
		x <- x[ind]
		sd <- sd[ind]
	}
	stopifnot(all(!is.na(x) & !is.na(sd)))

	# initializations
	p <- length(x)
	if (any(is.na(indices))) indices <- 1:p
	L <- rep(NA,length(indices)); U <- L

	# compute marginal CS for each population indicated by indices
	for (i in 1:length(indices)) {
		CS <- csranks_simul(x, sd, coverage=coverage, cstype=cstype, stepdown=stepdown, R=R, indices=indices[i], seed=seed)
		L[i] <- CS$L
		U[i] <- CS$U
	}
	
	return(list(L=as.integer(L),U=as.integer(U)))
}


#' Projection confidence sets for the tau-best
#'
#' @param x vector of estimates
#' @param sd vector of standard errors of \code{x}
#' @param tau the confidence set contains indicators for the elements in \code{x} whose rank is less than or equal to \code{tau}.
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param stepdown logical; if \code{TRUE} (default), stepwise procedure is used, otherwise single step procedure is used.
#' @param R number of bootstrap replications. Default is 1000.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{sd} (if any). 
#' @param seed seed for bootstrap random variable draws. If set to \code{NA} (default), then seed is not set.

#' @return logical vector indicating which of the elements of \code{x} are in the confidence set for the tau-best.

#' @section Details:
#' Implentation of the projection confidence sets for the tau-best proposed in Mogstad, Romano, Shaikh, and Wilhelm ("Inference for Ranks with Applications to Mobility across Neighborhoods and Academic Achievements across Countries", \href{https://www.ucl.ac.uk/~uctpdwi/papers/cwp1020.pdf}{CeMMAP Working Paper CWP10/20}). The confidence set contains indicators for the elements in \code{x} whose rank is less than or equal to \code{tau} with probability approximately equal to the coverage indicated in \code{coverage}. 
#' Parametric bootstrap based on the normal distribution with independent populations.

#' @examples
#' x <- seq(1,3,length=10)
#' sd <- rep(0.2,10)
#' cstaubest(x, sd, tau=3)

#' @export
cstaubest <- function(x, sd, tau=2, coverage=0.95, stepdown=TRUE, R=1000, na.rm=FALSE, seed=NA) {
	
	# return indices whose lower bound on the rank is <= tau
	L <- csranks_simul(x, sd, coverage=coverage, cstype="lower", stepdown=stepdown, R=R, indices=NA, na.rm=na.rm, seed=seed)$L
	return(L<=tau)
}


#' Projection confidence sets for the tau-worst
#'
#' @param x vector of estimates
#' @param sd vector of standard errors of \code{x}
#' @param tau the confidence set contains indicators for the elements in \code{x} whose rank is among the tau worst, i.e. larger than or equal to \code{length(x)-tau+1}.
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param stepdown logical; if \code{TRUE} (default), stepwise procedure is used, otherwise single step procedure is used.
#' @param R number of bootstrap replications. Default is 1000.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{sd} (if any). 
#' @param seed seed for bootstrap random variable draws. If set to \code{NA} (default), then seed is not set.

#' @return logical vector indicating which of the elements of \code{x} are in the confidence set for the tau-worst

#' @section Details:
#' Implentation of the projection confidence sets for the tau-worst proposed in Mogstad, Romano, Shaikh, and Wilhelm ("Inference for Ranks with Applications to Mobility across Neighborhoods and Academic Achievements across Countries", \href{https://www.ucl.ac.uk/~uctpdwi/papers/cwp1020.pdf}{CeMMAP Working Paper CWP10/20}). The confidence set contains indicators for the elements in \code{x} whose rank is larger than or equal to \code{length(x)-tau+1} with probability approximately equal to the coverage indicated in \code{coverage}. 
#' Parametric bootstrap based on the normal distribution with independent populations.

#' @examples
#' x <- seq(1,3,length=10)
#' sd <- rep(0.2,10)
#' cstauworst(x, sd, tau=3)

#' @export
cstauworst <- function(x, sd, tau=2, coverage=0.95, stepdown=TRUE, R=1000, na.rm=FALSE, seed=NA) {
	
	# return indices whose lower bound on the rank is <= tau
	U <- csranks_simul(x, sd, coverage=coverage, cstype="upper", stepdown=stepdown, R=R, indices=NA, na.rm=na.rm, seed=seed)$U
	p <- length(x)
	return(U>=p-tau+1)
}


	
