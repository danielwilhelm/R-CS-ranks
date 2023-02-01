#' Confidence sets for ranks
#'
#' Given estimates and their covariance matrix of a certain feature for a set of populations,
#' calculate confidence sets for the ranks of populations,
#' where populations are ranked by the feature values.
#'
#' @param x vector of estimates.
#' @param Sigma covariance matrix of \code{x}. Note, that it must be covariance matrix
#' of feature \bold{means}, not features themselves.
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param cstype type of confidence set (\code{two-sided}, \code{upper}, \code{lower}). Default is \code{two-sided}.
#' @param stepdown logical; if \code{TRUE} (default), stepwise procedure is used, otherwise single step procedure is used. See Details section for more.
#' @param R number of bootstrap replications. Default is 1000.
#' @param simul logical; if \code{TRUE} (default), then simultaneous confidence sets are computed, which jointly cover all populations indicated by \code{indices}.
#' 		Otherwise, for each population indicated in \code{indices} a marginal confidence set is computed.
#' @param indices vector of indices of \code{x} for whose ranks the confidence sets are computed. \code{indices=NA} (default) means computation for all ranks.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{Sigma} (if any).
#' @param seed seed for bootstrap random variable draws. If set to \code{NA} (default), then seed is not set.

#' @return A list with two items, \code{L} and \code{U} - lower and upper bounds of the confidence set for ranks indicated in \code{indices}.

#' @examples
#' # Setup example data
#' n <- 10
#' x <- seq(1, 3, length = n)
#' Sigma <- matrix(0.001, nrow = n, ncol = n)
#' diag(Sigma) <- 0.04
#' 
#' # Run csranks to get confidence sets for ranks of features
#' csranks(x, Sigma)
#' 
#' # If you assume that the feature measurements are independent 
#' # (or have access only to variances / standard errors estimates),
#' # then pass a diagonal covariance matrix.
#' Sigma <- diag(rep(0.04, 10))
#' csranks(x, Sigma)

#' @section Details:
#' IMPORTANT: make sure, that the \code{Sigma} is a (perhaps estimated) covariance matrix
#' of estimates of feature means across populations, not of features themselves.
#' For example, sample of size \eqn{n} of a feature following a standard normal distribution 
#' has variance \eqn{\sigma^2=1}, but mean from such sample has variance \eqn{1/n}.
#' We refer to the latter.
#' 
#' The command implements the procedure for construction of confidence sets for ranks described in the referenced paper below.
#' Generally, it consists of verification of a large set of hypotheses. After rejection of certain set
#' of hypotheses, one can terminate the procedure or keep verifying a smaller set of hypotheses that
#' were not rejected so far. The former corresponds to \code{stepdown=FALSE}; the latter to \code{stepdown=TRUE}.
#'
#' From a practical point of view, \code{stepdown=TRUE} takes more time, but usually results
#' in tighter (better) confidence sets.
#'
#' Parametric bootstrap used to calculate distribution for confidence sets based on the multivariate normal distribution.
#'
#' @references Mogstad, Romano, Shaikh, and Wilhelm (2023), "Inference for Ranks with Applications to Mobility across Neighborhoods and Academic Achievements across Countries", forthcoming at Review of Economic Studies
#' 
#' \href{http://dwilhelm.userweb.mwn.de/papers/cwp0323.pdf}{link to pdf}
#' @export
csranks <- function(x, Sigma, coverage = 0.95, cstype = "two-sided", stepdown = TRUE, R = 1000, simul = TRUE, indices = NA, na.rm = FALSE, seed = NA) {
  # initializations
  check_csranks_args(coverage=coverage, cstype=cstype, stepdown=stepdown, R=R,
                     simul=simul, seed=seed)
  l <- process_csranks_args(x, Sigma, indices, na.rm)
  x <- l$x; Sigma <- l$Sigma; indices <- l$indices
  
  if (simul) {
    return(csranks_simul(x, Sigma, coverage = coverage, cstype = cstype, stepdown = stepdown, R = R, indices = indices, na.rm = na.rm, seed = seed))
  } else {
    return(csranks_marg(x, Sigma, coverage = coverage, cstype = cstype, stepdown = stepdown, R = R, indices = indices, na.rm = na.rm, seed = seed))
  }
}

#' Simultaneous confidence sets for ranks 
#' 
#' This function is called by \code{csranks} when \code{simul=TRUE}.
#'
#' @noRd
csranks_simul <- function(x, Sigma, coverage = 0.95, cstype = "two-sided", stepdown = TRUE, R = 1000, indices = NA, na.rm = FALSE, seed = NA) {
  indices <- process_indices_argument(indices, length(x))
  # joint CS for difference in means
  csdifftype <- switch(cstype,
    "lower" = "upper",
    "upper" = "lower",
    "two-sided" = "symmetric"
  )
  res <- csdiffmeans(x, Sigma, coverage = coverage, indices = indices, cstype = csdifftype, stepdown = stepdown, R = R, seed = seed)
  L <- res$L
  U <- res$U

  # compute Nminus and Nplus
  # AKA the number of populations, for which the feature value
  # is larger / smaller in a statistically significant way
  if (stepdown & cstype == "two-sided") {
      Nplus <- rowSums(L[indices, , drop=FALSE] > 0, na.rm = TRUE)
      Nminus <- colSums(L[, indices, drop=FALSE] > 0, na.rm = TRUE)
  } else {
      Nplus <- rowSums(L[indices, , drop=FALSE] > 0, na.rm = TRUE)
      Nminus <- rowSums(U[indices, , drop=FALSE] < 0, na.rm = TRUE)
  }

  # return lower and upper confidence bounds for the ranks
  return(list(L = as.integer(Nminus + 1), U = as.integer(length(x) - Nplus)))
}

#' This function is called by \code{csranks} when \code{simul=FALSE}.
#' 
#' Marginal confidence sets for ranks
#' 
#' @noRd
csranks_marg <- function(x, Sigma, coverage = 0.95, cstype = "two-sided", stepdown = TRUE, R = 1000, indices = NA, na.rm = FALSE, seed = NA) {
  indices <- process_indices_argument(indices, length(x))
  L <- rep(NA, length(indices))
  U <- L

  # compute marginal CS for each population indicated by indices
  for (i in 1:length(indices)) {
    CS <- csranks_simul(x, Sigma, coverage = coverage, cstype = cstype, stepdown = stepdown, R = R, indices = indices[i], seed = seed)
    L[i] <- CS$L
    U[i] <- CS$U
  }

  return(list(L = as.integer(L), U = as.integer(U)))
}


#' Projection confidence sets for the tau-best
#'
#' Find a set of populations, which belong to tau-best populations according to
#' some feature with given confidence.
#'
#' @param tau the confidence set contains indicators for the elements in \code{x} whose rank is less than or equal to \code{tau}.
#' @inheritParams csranks
#' @return logical vector indicating which of the elements of \code{x} are in the confidence set for the tau-best.

#' @section Details:
#' The confidence set contains indicators for the elements in \code{x} whose rank is less than or equal to \code{tau} with probability approximately equal to the coverage indicated in \code{coverage}.
#' Parametric bootstrap based on the multivariate normal distribution.
#'
#' If \code{na.rm=TRUE} and NAs are present, then results are returned for tau-best (worst)
#' populations among those without NA values, i.e. after NA removal.
#' @examples
#' # Setup example data
#' n <- 10
#' x <- seq(1, 3, length = n)
#' Sigma <- matrix(0.001, nrow = n, ncol = n)
#' diag(Sigma) <- 0.04
#' 
#' # Run csranks to get confidence sets for top 3 populations
#' cstaubest(x, Sigma, tau = 3)
#' cstauworst(x, Sigma, tau = 3)
#' 
#' # If you assume that the feature measurements are independent, 
#' # (or just have access to variances / standard errors)
#' # then pass a diagonal covariance matrix.
#' Sigma <- diag(rep(0.04, 10))
#' cstaubest(x, Sigma, tau = 3)
#' cstauworst(x, Sigma, tau = 3)
#'
#' @inherit csranks references
#' @export
cstaubest <- function(x, Sigma, tau = 2, coverage = 0.95, stepdown = TRUE, R = 1000, na.rm = FALSE, seed = NA) {
  check_csranks_args(coverage=coverage, stepdown=stepdown, R=R, seed=seed, 
                     simul=TRUE, cstype="lower")
  check_tau(tau, length(x))
  l <- process_csranks_args(x, Sigma, NA, na.rm)
  x <- l$x; Sigma <- l$Sigma; indices <- l$indices
  
  # return indices whose lower bound on the rank is <= tau
  L <- csranks_simul(x, Sigma, coverage = coverage, cstype = "lower", stepdown = stepdown, R = R, indices = indices, na.rm = na.rm, seed = seed)$L
  return(L <= tau)
}


#' @describeIn cstaubest Projection confidence sets for the tau-worst
#'
#' Similar method, but for populations, which are tau-worst.
#' Equivalent to calling \code{cstaubest} with \code{-x}.
#'
#' @export
cstauworst <- function(x, Sigma, tau = 2, coverage = 0.95, stepdown = TRUE, R = 1000, na.rm = FALSE, seed = NA) {
  check_csranks_args(coverage=coverage, stepdown=stepdown, R=R, seed=seed, 
                     simul=TRUE, cstype="lower")
  check_tau(tau, length(x))
  l <- process_csranks_args(x, Sigma, NA, na.rm)
  x <- l$x; Sigma <- l$Sigma; indices <- l$indices
  # return indices whose lower bound on the rank is <= tau
  U <- csranks_simul(x, Sigma, coverage = coverage, cstype = "upper", stepdown = stepdown, R = R, indices = indices, na.rm = na.rm, seed = seed)$U
  p <- length(x)
  return(U >= p - tau + 1)
}