#' Confidence sets for ranks based on multinomial data
#'
#' Marginal and simultaneous confidence sets for ranks of categories, where categories are ranked by the probabilities of being chosen.
#'
#' @param x vector of counts indicating how often each category was chosen.
#' @param multcorr multiplicity correction to be used: \code{Holm} (default) or \code{Bonferroni}. See Details section for more.
#' @inheritParams csranks
#' @inherit csranks return
#' @section Details:
#' This function computes confidence sets for ranks similarly as \code{\link{csranks}}, but it is tailored to the special case of 
#' multinomial data. Suppose there are \eqn{p} populations (for the case of multinomial data, we will refer to them as "categories") such 
#' as political parties, for example, that one wants to rank by the probabilities of them being chosen. For political parties, this would
#' correspond to the share of votes each party obtains. Here, the underlying data are multinomial: each observation corresponds to a choice
#' among the \eqn{p} categories. The vector \code{x} contains the counts of how often each category was chosen in the data.
#' 
#' In this setting, \code{link{csranks}} could be applied to compute confidence sets for the ranks of each category, but instead this function
#' implements a different method proposed by Bazylik, Mogstad, Romano, Shaikh, and Wilhelm (2023), which exploits the 
#' multinomial structure of the problem and yields confidence sets for the ranks that are valid in finite samples (whereas \code{\link{csranks}} produces 
#' confidence sets that are valid only asymptotically).
#' 
#' The procedure involves testing multiple hypotheses. The `\code{multcorr}` indicates a method for multiplicity correction. See the paper for
#' details.
#'
#' @references
#' Bazylik, Mogstad, Romano, Shaikh, and Wilhelm.
#' "Finite-and large-sample inference for ranks using multinomial data with an application to ranking political parties".
#' \href{http://dwilhelm.userweb.mwn.de/papers/cwp4021.pdf}{cemmap working paper}
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
  x_ranks <- irank(x)[indices]
  if (simul) {
    confidence_set <- csranks_multinom_simul(x, coverage = coverage, cstype = cstype, multcorr = multcorr, indices = indices, na.rm = na.rm)
  } else {
    confidence_set <- csranks_multinom_marg(x, coverage = coverage, cstype = cstype, multcorr = multcorr, indices = indices, na.rm = na.rm)
  }
  structure(list(L = confidence_set$L,
                 rank = x_ranks,
                 U = confidence_set$U),
            class = "csranks")
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

#' @return logical matrix M. M[j,k] == TRUE means that we want to test hypothesis
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
  H0_smaller <- as.integer(l$requested_differences[,1])
  H0_larger <- as.integer(l$requested_differences[,2])
  
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
                 "Holm" = (1 - coverage) / seq(ncomp, 1, by = -1)
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
#' than j, so will have smaller rank
#' Nplus[j] - smaller, lower, larger
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
  indices <- process_indices_argument(indices, length(x))
  # compute marginal CS for each category indicated by indices
  LU <- sapply(indices, function(i){
    CS <- csranks_multinom_simul(x, coverage = coverage, cstype = cstype, indices = i)
    c(L = CS$L, U = CS$U)
  })
  
  return(list(L = as.integer(LU["L",]), U = as.integer(LU["U",])))
}
