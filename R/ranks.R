#' Confidence sets for ranks
#'
#' Marginal and simultaneous confidence sets for ranks.
#'
#' @param x vector of estimates containing estimated features by which the \code{length} populations are to be ranked.
#' @param Sigma estimated covariance matrix of \code{x}.
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param cstype type of confidence set (\code{two-sided}, \code{upper}, \code{lower}). Default is \code{two-sided}.
#' @param stepdown logical; if \code{TRUE} (default), stepwise procedure is used, otherwise single step procedure is used. See Details section for more.
#' @param R number of bootstrap replications. Default is 1000.
#' @param simul logical; if \code{TRUE} (default), then simultaneous confidence sets are computed, which jointly cover all populations indicated by \code{indices}.
#' 		Otherwise, for each population indicated in \code{indices} a marginal confidence set is computed.
#' @param indices vector of indices of \code{x} for whose ranks the confidence sets are computed. \code{indices=NA} (default) means computation of confidence sets for all populations.
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{Sigma} (if any).
#' @param seed seed for bootstrap random variable draws. If set to \code{NA} (default), then seed is not set.

#' @return A \code{csranks} object, which is a list with three items:
#' \describe{
#'  \item{\code{L}}{Lower bounds of the confidence sets for ranks indicated in \code{indices}}
#'  \item{\code{rank}}{Estimated ranks from \code{\link{irank}} with default parameters}
#'  \item{\code{U}}{Upper bounds of the confidence sets.}
#' }

#' @examples
#' # simple simulated example:
#' n <- 100
#' p <- 10
#' X <- matrix(rep(1:p,n)/p, ncol=p, byrow=TRUE) + matrix(rnorm(n*p), 100, 10)
#' thetahat <- colMeans(X)
#' Sigmahat <- cov(X) / n
#' csranks(thetahat, Sigmahat)
#' 
#' # PISA example:
#' attach(pisa)
#' math_cov_mat <- diag(math_se^2)
#' 
#' # marginal confidence set for each country:
#' csranks(math_score, math_cov_mat, simul=FALSE)
#' 
#' # simultaneous confidence set for all countries:
#' csranks(math_score, math_cov_mat, simul=TRUE)


#' @section Details:
#' Suppose \eqn{j=1,\ldots,p} populations (e.g., schools, hospitals, political parties, countries) are to be ranked according to 
#' some measure \eqn{\theta=(\theta_1,\ldots,\theta_p)}. We do not observe the true values \eqn{\theta_1,\ldots,\theta_p}. Instead, for each population, 
#' we have data from which we have estimated these measures, \eqn{\hat{\theta}=(\hat{\theta}_1,\ldots,\hat{\theta}_p)}. The values \eqn{\hat{\theta}_1,\ldots,\hat{\theta}_p}
#' are estimates of the true values \eqn{\theta_1,\ldots,\theta_p} and thus contain statistical uncertainty. In consequence, a ranking of the populations by
#' the values \eqn{\hat{\theta}_1,\ldots,\hat{\theta}_p} contains statistical uncertainty and is not necessarily equal to the true ranking of \eqn{\theta_1,\ldots,\theta_p}.
#' 
#' The function computes confidence sets for the rank of one, several or all of the populations (\code{indices} indicates which of the \eqn{1,\ldots,p} populations are of interest). \code{x} is a vector containing the estimates 
#' \eqn{\hat{\theta}_1,\ldots,\hat{\theta}_p} and \code{Sigma} is an estimate of the covariance matrix of \code{x}. The method assumes that the estimates are asymptotically normal and the sample sizes of the datasets 
#' are large enough so that \eqn{\hat{\theta}-\theta} is approximately distributed as \eqn{N(0,\Sigma)}. The argument \code{Sigma} should contain an estimate of the covariance matrix \eqn{\Sigma}. For instance, if for each population \eqn{j}
#' \deqn{\sqrt{n_j} (\hat{\theta}_j-\theta_j) \to_d N(0, \sigma_j^2)}
#' and the datasets for each population are drawn independently of each other, then \code{Sigma} is a diagonal matrix \deqn{diag(\hat{\sigma}_1^2/n_1,\ldots,\hat{\sigma}_p^2/n_p)}
#' containing estimates of the asymptotic variances divided by the sample size. More generally, the estimates in \code{x} may be dependent, but then \code{Sigma}
#' must be an estimate of its covariance matrix including off-diagonal terms. 
#' 
#' Marginal confidence sets (\code{simul=FALSE}) are such that the confidence set for a population \eqn{j} contains the true rank of that population \eqn{j} with probability approximately
#' equal to the nominal coverage level. Simultaneous confidence sets (\code{simul=TRUE}) on the other hand are such that the confidence sets for populations indicated in \code{indices} cover the true ranks
#' of all of these populations simultaneously with probability approximately equal to the nominal coverage level. For instance, in the PISA example below, a marginal confidence set of a country \eqn{j} covers the true
#' rank of country \eqn{j} with probability approximately equal to 0.95. A simultaneous confidence set for all countries covers the true ranks of all countries simultaneously with probability approximately equal to 0.95.
#' 
#' The function implements the procedures developed and described in more detail in Mogstad, Romano, Shaikh, and Wilhelm (2023). The procedure is based on
#' on testing a large family of hypotheses for pairwise comparisons. Stepwise methods can be used to improve the power of the procedure by, potentially,
#' rejecting more hypotheses without violating the desired coverage property of the resulting confidence set. These are employed when
#' \code{stepdown=TRUE}. From a practical point of view, \code{stepdown=TRUE} is computationally more demanding, but often results
#' in tighter confidence sets.
#'
#' The procedure uses a parametric bootstrap procedure based on the above approximate multivariate normal distribution.
#'
#' @references Mogstad, Romano, Shaikh, and Wilhelm (2023), "Inference for Ranks with Applications to Mobility across Neighborhoods and Academic Achievements across Countries", forthcoming at Review of Economic Studies
#' \href{http://dwilhelm.userweb.mwn.de/papers/cwp0323.pdf}{cemmap working paper}
#' \doi{10.1093/restud/rdad006}
#' @export
csranks <- function(x, Sigma, coverage = 0.95, cstype = "two-sided", stepdown = TRUE, R = 1000, simul = TRUE, indices = NA, na.rm = FALSE, seed = NA) {
  # initializations
  check_csranks_args(coverage=coverage, cstype=cstype, stepdown=stepdown, R=R,
                     simul=simul, seed=seed)
  l <- process_csranks_args(x, Sigma, indices, na.rm)
  x <- l$x; Sigma <- l$Sigma; indices <- l$indices
  x_ranks <- irank(x)[indices]
  if (simul) {
    confidence_set <- csranks_simul(x, Sigma, coverage = coverage, cstype = cstype, stepdown = stepdown, R = R, indices = indices, na.rm = na.rm, seed = seed)
  } else {
    confidence_set <- csranks_marg(x, Sigma, coverage = coverage, cstype = cstype, stepdown = stepdown, R = R, indices = indices, na.rm = na.rm, seed = seed)
  }
  structure(list(L = confidence_set$L,
            rank = x_ranks,
            U = confidence_set$U),
            class = "csranks")
}

#' Simultaneous confidence sets for ranks 
#' 
#' This function is called by \code{csranks} when \code{simul=TRUE}.
#'
#' @noRd
csranks_simul <- function(x, Sigma, coverage = 0.95, cstype = "two-sided", stepdown = TRUE, R = 1000, indices = NA, na.rm = FALSE, seed = NA) {
  indices <- process_indices_argument(indices, length(x))
  # joint CS for difference in means
  # switching cause lower rank is associated with higher values in x
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
  
  # compute marginal CS for each population indicated by indices
  LU <- sapply(indices, function(i){
    CS <- csranks_simul(x, Sigma, coverage = coverage, cstype = cstype, stepdown = stepdown, R = R, indices = i, seed = seed)
    c(L = CS$L, U = CS$U)
  })

  return(list(L = as.integer(LU["L",]), U = as.integer(LU["U",])))
}


#' Confidence sets for the tau-best
#'
#' Computation of confidence sets for the identities of populations among the tau best.
#'
#' @param tau the confidence set contains indicators for the elements in \code{x} whose rank is less than or equal to \code{tau}.
#' @inheritParams csranks
#' @return logical vector indicating which of the elements of \code{x} are in the confidence set for the tau-best.

#' @section Details:
#' The function computes a confidence set containing indicators for the elements in \code{x} whose rank is less than or equal to \code{tau} with probability approximately equal to the nominal coverage (\code{coverage}).
#' 
#' The function implements the projection confidence set for the tau-best developed and described in more detail in Mogstad, Romano, Shaikh, and Wilhelm (2023).

#' @references Mogstad, Romano, Shaikh, and Wilhelm (2023), "Inference for Ranks with Applications to Mobility across Neighborhoods and Academic Achievements across Countries", forthcoming at Review of Economic Studies
#' \href{http://dwilhelm.userweb.mwn.de/papers/cwp0323.pdf}{cemmap working paper}, \doi{10.1093/restud/rdad006}
#' @examples
#' # simple simulated example:
#' n <- 100
#' p <- 10
#' X <- matrix(rep(1:p,n)/p, ncol=p, byrow=TRUE) + matrix(rnorm(n*p), 100, 10)
#' thetahat <- colMeans(X)
#' Sigmahat <- cov(X) / n
#' 
#' # confidence set for the populations that may be among the top-3 
#' # (with probability approximately 0.95):
#' cstaubest(thetahat, Sigmahat, tau=3)
#' 
#' # confidence set for the populations that may be among the bottom-3 
#' # (with probability approximately 0.95):
#' cstauworst(thetahat, Sigmahat, tau=3)
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


#' @describeIn cstaubest Confidence sets for the tau-worst
#'
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
