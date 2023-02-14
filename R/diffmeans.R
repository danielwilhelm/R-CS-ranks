#' Confidence sets for vector of differences
#' @inheritParams csranks
#' 
#' @return Confidence intervals for diferences between feature values in x.
#' In a form of a list with two items, \code{L} and \code{U}.
#' Each is a matrix of size length(x) times length(x).
#' The confidence interval for a difference between feature from ith and jth 
#' population is [out$L[i,j], out$U[i,j]].
#' NAs are returned whenever the value is not of interest (one-sided CI or 
#' i,j not in indices)
#' 
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @noRd
csdiffmeans <- function(x, cov_mat, coverage = 0.95, indices = NA, cstype = "symmetric", stepdown = TRUE, R = 1000, seed = NA) {
  # check arguments
  cstype <- match.arg(cstype, c("symmetric", "upper", "lower"))
  
  # initializations
  p <- length(x)
  if (any(is.na(indices))) indices <- 1:p
  thetadiff <- outer(x, x, "-")
  sigmadiff <- calculate_difference_sds(cov_mat)
  anyrejections <- TRUE
  # Bounds for confidence intervals of differences
  ChatL <- matrix(0, p, p)
  ChatU <- ChatL
  maxabs <- function(x) max(abs(x))
  mmax <- function(x) max(-x)

  # which differences?
  l <- initialize_I0(p=p, indices=indices, stepdown = stepdown, cstype=cstype)
  I0 <- l$I0; cstype <- l$cstype
  I1 <- I0
  Istar <- I0

  # beta-quantiles from Ln and Un
  # using parametric bootstrap
  LInv <- function(I, beta, fn){
    educed_I <- reduce_I(I)
    needed_variables <- reduced_I$needed_variables
    requested_diffrences <- reduced_I$requested_diffrences
    # Marginal distribution of MVTNormal is MVTNormal 
    # with correct block of covariance matrix
    needed_cov_mat <- cov_mat[needed_variables, needed_variables]
    Z <- MASS::mvrnorm(R, mu = rep(0, nrow(needed_cov_mat)),
                       Sigma = needed_cov_mat)
    Zdiff_scaled <- calculate_scaled_differences_in_samples(Z, requested_diffrences,
                                                            sigmadiff[I])
    bootstrap_estimates <- apply(Zdiff_scaled, 1, fn)
    quantile(bootstrap_estimates, probs=beta)
  }
  LLowerInv <- function(I, beta) LInv(I, beta, max)
  LUpperInv <- function(I, beta) LInv(I, beta, mmax)
  LSymmInv <- function(I, beta) LInv(I, beta, maxabs)

  # compute CS
  if (!is.na(seed)) set.seed(seed)
  while (anyrejections) {
    # compute upper and lower confidence bounds
    if (cstype == "lower") {
      ChatL[I1] <- thetadiff[I1] - sigmadiff[I1] * LLowerInv(I1, coverage)
      ChatU <- matrix(Inf, p, p)
    }
    if (cstype == "upper") {
      ChatL <- matrix(-Inf, p, p)
      ChatU[I1] <- thetadiff[I1] + sigmadiff[I1] * LUpperInv(I1, coverage)
    }
    if (cstype == "symmetric") {
      z <- LSymmInv(I1, coverage)
      ChatL[I1] <- thetadiff[I1] - sigmadiff[I1] * z
      ChatU[I1] <- thetadiff[I1] + sigmadiff[I1] * z
    }

    # remove indices for which interval does not contain 0
    I1[ChatU < 0 | ChatL > 0] <- FALSE

    # stepdown improvements?
    if (stepdown) {
      # any new rejections?
      if (sum(I1) == sum(I0) | sum(I1) == 0) {
        anyrejections <- FALSE
      } else {
        I0 <- I1
      }
    } else {
      anyrejections <- FALSE
    }
  }

  ChatL[!Istar] <- NA
  ChatU[!Istar] <- NA

  return(list(L = ChatL, U = ChatU))
}


#' Calculate standard deviations of diffrences in an MVTNormal distribution
#' @param cov_mat a covariance matrix of multivariate normal distribution
#' @return matrix pxp with standard deviations of differences of variables
#' @noRd
calculate_difference_sds <- function(cov_mat){
  # Var[X_1 - X_2] = Var[X_1] + Var[X_2] - 2Cov[X_1,X_2]
  # And a difference of correlated Gaussians is still Gaussian
  sqrt(outer(diag(cov_mat), diag(cov_mat), "+") - 2 * cov_mat)
}

initialize_I0 <- function(p, indices, stepdown, cstype){
  I0 <- matrix(TRUE, p, p)
  I0[-indices,] <- FALSE
  if (stepdown & cstype == "symmetric") {
    # for other cstypes, correction is unnecessary
    I0[, indices] <- TRUE
    cstype <- "lower"
  }
  diag(I0) <- FALSE
  list(I0 = I0, cstype = cstype)
}

#' Convert information about needed differences
#' from a boolean matrix to more useful form
#' @return a list with `needed_variables` - boolean
#' and `requested differences` - matrix with 2 columns. Each row corresponds
#' to single TRUE entraince in I. Its contents are the indices of the entry.
#' 
#' @examples
#' I <- matrix(c(FALSE, FALSE, FALSE,
#'               FALSE, FALSE, TRUE,
#'               FALSE, FALSE, FALSE), byrow = TRUE, ncol = 3)
#' r <- reduce_I(I)
#' r$needed_variables # c(FALSE, TRUE, TRUE)
#' r$requested_diffrences # matrix(c(2,3), ncol = 2)
#' @noRd
reduce_I <- function(I){
  needed_variables <- sapply(1:nrow(I), function(i){
    any(I[i,]) || any(I[,i])
  })
  needed_I <- I[needed_variables, needed_variables]
  requested_diffrences <- get_double_from_single_indices(which(needed_I), 
                                                         nrow(needed_I))
  list(needed_variables = needed_variables,
       requested_diffrences = requested_diffrences)
}

#' @param Z sample from mvt normal with observations in rows
#' @param requested_differences matrix with 2 columns
#' we want to calculate differences between z_i and z_j variables
#' for i,j in rows of requested_differences
#' for each observation z in Z 
#' And then scale them with
#' @param scales: numeric of length nrow(requested_differences)
#' @noRd
calculate_scaled_differences_in_samples <- function(Z, requested_diffrences, scales){
  Zdiff <- Z[, requested_diffrences[, 1]] - Z[, requested_diffrences[, 2]]
  # Vectorized division goes over rows
  Zdiff_scaled <- t(t(Zdiff) / scales) 
  Zdiff_scaled
}


