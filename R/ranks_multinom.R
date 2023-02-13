#' Confidence sets for ranks based on multinomial data
#'
#' Given data on counts of successes for each category, calculate confidence sets for the ranks of categories, where categories are ranked by their success probabilities.
#'
#' @param x vector of counts of successes for each category
#' @param multcorr multiplicity correction to be used: \code{Holm} (default), \code{Bonferroni} or \code{none}. See Details section for more.
#' @inheritParams csranks
#' @inherit csranks return
#' @section Details:
#' The command implements the procedure for construction of confidence sets for ranks described in the referenced paper below.
#' 
#' It involves testing multiple hypotheses. The `\code{multcorr}` states, how the 
#' p-values should be corrected to control the Family Wise Error Rate (FWER).
#' Not applying correction (\code{multcorr=none}) is not advised.
#'
#' From a practical point of view, \code{multcorr=Holm} takes more time, but usually results
#' in tighter (better) confidence sets than \code{multcorr=Bonferroni}.
#'
#' @references
#' Bazylik, Mogstad, Romano, Shaikh, and Wilhelm.
#' "Finite-and large-sample inference for ranks using multinomial data with an application to ranking political parties".
#'
#' @examples
#' x <- c(rmultinom(1, 1000, 1:10))
#' csranks_multinom(x)
#' @export
csranks_multinom <- function(x, coverage = 0.95, cstype = "two-sided", simul = TRUE, multcorr = "Holm", indices = NA, na.rm = FALSE) {
  # initializations
  check_csranks_multinom_args(coverage=coverage, cstype=cstype, simul=simul, multcorr=multcorr)
  l <- process_csranks_multinom_args(x, indices, na.rm)
  x <- l$x; indices <- l$indices
  
  if (simul) {
    return(csranks_multinom_simul(x, coverage = coverage, cstype = cstype, multcorr = multcorr, indices = indices, na.rm = na.rm))
  } else {
    return(csranks_multinom_marg(x, coverage = coverage, cstype = cstype, multcorr = multcorr, indices = indices, na.rm = na.rm))
  }
}

#' Simultaneous Confidence sets for ranks based on multinomial data
#'
#' This function is called by \code{csranks_multinom} when \code{simul=TRUE}.
#'
#' @noRd
#' @importFrom stats rmultinom
#' @importFrom stats pbinom
#' @importFrom stats aggregate
csranks_multinom_simul <- function(x, coverage = 0.95, cstype = "two-sided", multcorr = "Holm", indices = NA, na.rm = FALSE) {
  
  indices <- process_indices_argument(indices, length(x))
  p <- length(x)
  
  which_pairs <- which_pairs_are_compared(p, indices, cstype)
  pairwise_p_values <- calculate_pairwise_p_values(x, which_pairs)
  rejection_results <- reject_or_accept(pairwise_p_values, multcorr, coverage)
  Nlist <- calculate_N_plus_minus(rejection_results, cstype, indices)
  convert_N_plus_minus_to_csrank(Nlist$Nminus, Nlist$Nplus, p)
}

#' @return boolean matrix M. M[j,k] == TRUE means that we want to test hypothesis
#' that x_j <= x_k
#' @noRd
which_pairs_are_compared <- function(p, indices, cstype){
  ind <- matrix(TRUE, p, p)
  if (cstype == "two-sided") {
    ind[-indices, ] <- FALSE
    ind[, indices] <- TRUE
  }
  if (cstype == "lower") {
    ind[, -indices] <- FALSE
  }
  if (cstype == "upper") {
    ind[-indices, ] <- FALSE
  }
  diag(ind) <- FALSE
  ind
}

calculate_pairwise_p_values <- function(x, which_pairs){
  l <- reduce_I(which_pairs)
  # H0 - under null hypothesis
  H0_smaller <- l$requested_diffrences[,1]
  H0_larger <- l$requested_diffrences[,2]
  
  x_parwise_sums <- outer(x, x, "+")
  relevant_parwise_sums <- as.vector(x_parwise_sums[which_pairs])
  
  pval <- calculate_p_value(x[H0_smaller], relevant_parwise_sums)
  
  data.frame(j = H0_smaller, 
             k = H0_larger,
             pval = pval)
}

#' x | S ~ binom(S, p) 
#' Test with H0: p <= 0.5
#' @noRd
calculate_p_value <- function(x, S){
  # X ~ binom(S, 0.5)
  # Probability of observing X >= x under H0, written as P(X >= x), is
  # 1 - P(X < x) = 1 - P(X <= x-1)
  1 - pbinom(x - 1, S, 0.5)
}

reject_or_accept <- function(df, multcorr, coverage){
  # order pairs by their p-values (from smallest to largest)
  df <- df[order(df$pval), ]
  
  # compute critical value
  ncomp <- nrow(df)
  beta <- switch(multcorr,
                 "Bonferroni" = (1 - coverage) / ncomp,
                 "Holm" = (1 - coverage) / (ncomp + 1 - (1:ncomp)),
                 "none" = 1 - coverage
  )
  
  # perform test(s)
  df$rej <- df$pval <= beta
  if (multcorr == "Holm") {
    ind_first_false <- which(!df$rej)[1]
    if (!is.na(ind_first_false)) df$rej[ind_first_false:length(df$rej)] <- FALSE
  }
  df
}

#' @return list with Nminus and Nplus, numeric vectors
#' Nminus[j] is the amount of populations, which have significantly larger 
#' ranking feature that population j
#' Thus, it's the number of populations that will be confidently higher in ranking
#' than j
#' Nplus[j] - smaller
#' @noRd
calculate_N_plus_minus <- function(df, cstype, indices){
  if (cstype == "upper" | cstype == "two-sided") {
    rej_plus <- aggregate(df["rej"], by = df["j"], sum)
    ind_plus <- rej_plus$j %in% indices
    Nplus <- rej_plus$rej[ind_plus]
  } else {
    Nplus <- 0 * indices
  }
  
  if (cstype == "lower" | cstype == "two-sided") {
    rej_minus <- aggregate(df["rej"], by = df["k"], sum)
    ind_minus <- rej_minus$k %in% indices
    Nminus <- rej_minus$rej[ind_minus]
  } else {
    Nminus <- indices * 0
  }
  list(Nminus = Nminus, Nplus = Nplus)
}

convert_N_plus_minus_to_csrank <- function(Nminus, Nplus, p){
  L <- Nminus + 1
  U <- p - Nplus
  list(L = as.integer(L), U = as.integer(U))
}


#' Marginal Confidence sets for ranks based on multinomial data
#'
#' This function is called by \code{csranks_multinom} when \code{simul=FALSE}.
#'
#' @noRd
csranks_multinom_marg <- function(x, coverage = 0.95, cstype = "two-sided", multcorr = "Holm", indices = NA, na.rm = FALSE) {
  
  L <- rep(NA, length(indices))
  U <- L

  # compute marginal CS for each category indicated by indices
  for (i in 1:length(indices)) {
    CS <- csranks_multinom_simul(x, coverage = coverage, cstype = cstype, multcorr = multcorr, indices = indices[i])
    L[i] <- CS$L
    U[i] <- CS$U
  }

  return(list(L = as.integer(L), U = as.integer(U)))
}
