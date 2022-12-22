#' Confidence sets for ranks
#'
#' Given estimates and their standard errors of a certain feature for a set of populations,
#' calculate confidence sets for the ranks of populations,
#' where populations are ranked by the feature values.
#'
#' @param x vector of estimates
#' @param sd vector of standard errors of \code{x}
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param cstype type of confidence set (\code{two-sided}, \code{upper}, \code{lower}). Default is \code{two-sided}.
#' @param stepdown logical; if \code{TRUE} (default), stepwise procedure is used, otherwise single step procedure is used. See Details section for more.
#' @param R number of bootstrap replications. Default is 1000.
#' @param simul logical; if \code{TRUE} (default), then simultaneous confidence sets are computed, which jointly cover all populations indicated by \code{indices}.
#' 		Otherwise, for each population indicated in \code{indices} a marginal confidence set is computed.
#' @param indices vector of indices of \code{x} for whose ranks the confidence sets are computed. \code{indices=NA} (default) means computation for all ranks.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{sd} (if any).
#' @param seed seed for bootstrap random variable draws. If set to \code{NA} (default), then seed is not set.

#' @return A list with two items, `L` and `U` - lower and upper bounds of the confidence set for ranks indicated in \code{indices}.

#' @examples
#' x <- seq(1, 3, length = 10)
#' sd <- rep(0.2, 10)
#' csranks(x, sd)
#' csranks_simul(x, sd)
#' csranks_marg(x, sd)

#' @section Details:
#' The command implements the procedure for construction of confidence sets for ranks described in the referenced paper below.
#' Generally, it consists of verification of a large set of hypotheses. After rejection of certain set
#' of hypotheses, one can terminate the procedure or keep verifying a smaller set of hypotheses that
#' were not rejected so far. The former corresponds to \code{stepdown=FALSE}; the latter to \code{stepdown=TRUE}.
#'
#' From a practical point of view, \code{stepdown=TRUE} takes more time, but usually results
#' in tighter (better) confidence sets.
#'
#' Parametric bootstrap used to calculate distribution for confidence sets based on the normal distribution with independent populations.
#'
#' @references Mogstad, Romano, Shaikh, and Wilhelm.
#' "Inference for Ranks with Applications to Mobility across Neighborhoods and Academic Achievements across Countries"
#' \href{https://www.ucl.ac.uk/~uctpdwi/papers/cwp1020.pdf}{CeMMAP Working Paper CWP10/20}
#' @export
csranks <- function(x, sd, coverage = 0.95, cstype = "two-sided", stepdown = TRUE, R = 1000, simul = TRUE, indices = NA, na.rm = FALSE, seed = NA) {
  if (simul) {
    return(csranks_simul(x, sd, coverage = coverage, cstype = cstype, stepdown = stepdown, R = R, indices = indices, na.rm = na.rm, seed = seed))
  } else {
    return(csranks_marg(x, sd, coverage = coverage, cstype = cstype, stepdown = stepdown, R = R, indices = indices, na.rm = na.rm, seed = seed))
  }
}

#' @describeIn csranks Simultaneous confidence sets for ranks
#'
#' This function is called by \code{csranks} when \code{simul=TRUE}.
#'
#' @export
csranks_simul <- function(x, sd, coverage = 0.95, cstype = "two-sided", stepdown = TRUE, R = 1000, indices = NA, na.rm = FALSE, seed = NA) {
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
    "two-sided" = "symmetric"
  )
  p <- length(x)
  if (any(is.na(indices))) indices <- 1:p
  res <- csdiffmeans(x, sd, coverage = coverage, indices = indices, cstype = csdifftype, stepdown = stepdown, R = R, seed = seed)
  L <- res$L
  U <- res$U

  # compute Nminus and Nplus
  if (stepdown & cstype == "two-sided") {
    if (length(indices) == 1 & all(!is.na(indices))) {
      Nplus <- sum(L[indices, ] > 0, na.rm = TRUE)
      Nminus <- sum(L[, indices] > 0, na.rm = TRUE)
    } else {
      Nplus <- rowSums(L[indices, ] > 0, na.rm = TRUE)
      Nminus <- colSums(L[, indices] > 0, na.rm = TRUE)
    }
  } else {
    if (length(indices) == 1 & all(!is.na(indices))) {
      Nplus <- sum(L[indices, ] > 0, na.rm = TRUE)
      Nminus <- sum(U[indices, ] < 0, na.rm = TRUE)
    } else {
      Nplus <- rowSums(L[indices, ] > 0, na.rm = TRUE)
      Nminus <- rowSums(U[indices, ] < 0, na.rm = TRUE)
    }
  }

  # return lower and upper confidence bounds for the ranks
  return(list(L = as.integer(Nminus + 1), U = as.integer(p - Nplus)))
}

#' @describeIn csranks Marginal confidence sets for ranks
#'
#' This function is called by \code{csranks} when \code{simul=FALSE}.
#'
#' @export
csranks_marg <- function(x, sd, coverage = 0.95, cstype = "two-sided", stepdown = TRUE, R = 1000, indices = NA, na.rm = FALSE, seed = NA) {
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
  L <- rep(NA, length(indices))
  U <- L

  # compute marginal CS for each population indicated by indices
  for (i in 1:length(indices)) {
    CS <- csranks_simul(x, sd, coverage = coverage, cstype = cstype, stepdown = stepdown, R = R, indices = indices[i], seed = seed)
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
#' Parametric bootstrap based on the normal distribution with independent populations.

#' @examples
#' x <- seq(1, 3, length = 10)
#' sd <- rep(0.2, 10)
#' cstaubest(x, sd, tau = 3)
#' cstauworst(x, sd, tau = 3)
#'
#' @inherit csranks references
#' @export
cstaubest <- function(x, sd, tau = 2, coverage = 0.95, stepdown = TRUE, R = 1000, na.rm = FALSE, seed = NA) {
  # return indices whose lower bound on the rank is <= tau
  L <- csranks_simul(x, sd, coverage = coverage, cstype = "lower", stepdown = stepdown, R = R, indices = NA, na.rm = na.rm, seed = seed)$L
  return(L <= tau)
}


#' @describeIn cstaubest Projection confidence sets for the tau-worst
#'
#' Similar method, but for populations, which are tau-worst.
#' Equivalent to calling \code{cstaubest} with \code{-x}.
#'
#' @export
cstauworst <- function(x, sd, tau = 2, coverage = 0.95, stepdown = TRUE, R = 1000, na.rm = FALSE, seed = NA) {
  # return indices whose lower bound on the rank is <= tau
  U <- csranks_simul(x, sd, coverage = coverage, cstype = "upper", stepdown = stepdown, R = R, indices = NA, na.rm = na.rm, seed = seed)$U
  p <- length(x)
  return(U >= p - tau + 1)
}
