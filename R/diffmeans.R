#' Confidence sets for vector of differences
#'
#' @param x vector of estimates
#' @param sd vector of standard errors of \code{x}
#' @param coverage nominal coverage of the confidence set. Default is 0.95.
#' @param indices vector of indices of \code{x} whose differences with all others the confidence set jointly covers. \code{indices=NA} (default) means confidence set jointly covers all possible differences.
#' @param cstype type of confidence set (\code{symmetric}, \code{upper}, \code{lower}). Default is \code{symmetric}.
#' @param stepdown logical; if \code{TRUE} (default), stepwise procedure is used, otherwise single step procedure is used.
#' @param R number of bootstrap replications. Default is 1000.
#' @param seed seed for bootstrap random variable draws. If set to \code{NA} (default), then seed is not set.

#' @return L,U lower and upper bounds of the confidence set for ranks indicated in \code{indices}.

#' @importFrom stats quantile
#' @importFrom stats rnorm
csdiffmeans <- function(x, sd, coverage=0.95, indices=NA, cstype="symmetric", stepdown=TRUE, R=1000, seed=NA) {

  # check arguments
  cstype <- match.arg(cstype, c("symmetric", "upper", "lower"))

  # initializations
  p <- length(x)
  if (any(is.na(indices))) indices <- 1:p
  thetadiff <- outer(x, x, '-')
  sigmadiff <- sqrt(outer(sd^2, sd^2, '+'))
  anyrejections <- TRUE
  ChatL <- matrix(0,p,p); ChatU <- ChatL
  maxabs <- function(x) max(abs(x))
  mmax <- function(x) max(-x)

  # which differences?
  I0 <- matrix(TRUE, p, p); I0[-indices,] <- FALSE;
  if (stepdown & cstype=="symmetric") {
    cstype <- "lower"
    I0[,indices] <- TRUE
  }
  diag(I0) <- FALSE; I1 <- I0; Istar <- I0

  # parametric bootstrap
  diffFn <- function(Z, I, fn) {
    Zdiff <- outer(Z, Z, '-')
    fn(Zdiff[I]/sigmadiff[I])
  }

  # beta-quantiles from Ln and Un
  LLowerInv <- function(I, beta) quantile(replicate(R, diffFn(sd*rnorm(p), I, max)), probs=beta)
  LUpperInv <- function(I, beta) quantile(replicate(R, diffFn(sd*rnorm(p), I, mmax)), probs=beta)
  LSymmInv <- function(I, beta) quantile(replicate(R, diffFn(sd*rnorm(p), I, maxabs)), probs=beta)

  # compute CS
  if (!is.na(seed)) set.seed(seed)
  while (anyrejections) {

    # compute upper and lower confidence bounds
    if (cstype=="lower") {
      ChatL[I1] <- thetadiff[I1] - sigmadiff[I1]*LLowerInv(I1, coverage)
      ChatU <- matrix(Inf, p, p)
    }
    if (cstype=="upper") {
      ChatL <- matrix(-Inf, p, p)
      ChatU[I1] <- thetadiff[I1] + sigmadiff[I1]*LUpperInv(I1, coverage)
    }
    if (cstype=="symmetric") {
      z <- LSymmInv(I1, coverage)
      ChatL[I1] <- thetadiff[I1] - sigmadiff[I1]*z
      ChatU[I1] <- thetadiff[I1] + sigmadiff[I1]*z
    }

    # remove indices for which interval does not contain 0
    I1[ChatU < 0 | ChatL > 0] <- FALSE

    # stepdown improvements?
    if (stepdown) {
      # any new rejections?
      if (sum(I1)==sum(I0) | sum(I1)==0) {
        anyrejections <- FALSE
      } else {
        I0 <- I1
      }
    } else anyrejections <- FALSE
  }

  ChatL[!Istar] <- NA
  ChatU[!Istar] <- NA

  return(list(L=ChatL,U=ChatU))
}
