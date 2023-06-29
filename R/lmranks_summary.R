#' Summarizing fits of linear models for ranks
#' 
#' Summary method for class "\code{lmranks}". It returns theoretically valid standard
#' errors, in comparison to naively running \code{summary(lm(...))}.
#' 
#' @param object A \code{lmranks} object.
#' @inheritParams stats::summary.lm
#' @export
summary.lmranks <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...){
  if(symbolic.cor){
    cli::cli_abort("{.var symbolic.cor} are not yet implemented for {.class lmranks}.")
  }
  # call summary.lm
  object$df.residual <- nrow(object$model) - ncol(object$model)
  outcome <- NextMethod()
  object$df.residual <- NA
  
  # Mark what is unknown (for now)
  outcome$coefficients[,2:4] <- NA
  outcome$sigma <- NA
  # This one causes errors in print.summary.lm
  # If needed, we could Ctrl-C Ctrl-V and adapt the method
  #outcome$df <- c(NA, nrow(object$model) - ncol(object$model), NA)
  outcome$fstatistic <- NULL
  # Remember to handle this, once fstatistic is known
  outcome$adj.r.squared <- NA
  cov_matrix <- vcov(object, complete=FALSE)
  outcome$cov.unscaled <- matrix(NA, nrow = nrow(outcome$cov.unscaled),
                                 ncol = ncol(outcome$cov.unscaled))
  
  outcome$coefficients[, 2] <- sqrt(diag(cov_matrix))
  outcome$coefficients[, 3] <- outcome$coefficients[, 1] / outcome$coefficients[, 2]
  outcome$coefficients[, 4] <- 2*stats::pnorm(-abs(outcome$coefficients[, 3]))
  
  if(correlation)
    outcome$correlation <- stats::cov2cor(cov_matrix)
  cli::cli_warn(c("The number of residual degrees of freedom is not correct.", 
                "Also, z-value, not t-value, since the distribution used for p-value calculation is standard normal."))
  class(outcome) <- c("summary.lmranks", class(outcome))
  outcome
}

#' @export
confint.lmranks <- function(object, parm, level = 0.95, ...){
  # As is the case with confint.lm, this method returns *marginal* CIs for coefficients
  # not simultaneous
  if(missing(parm))
    stats::confint.default(object=object, level=level, ...)
  else
    stats::confint.default(object=object, parm=parm, level=level, ...)
}

#' @describeIn lmranks Calculate Variance-Covariance Matrix for a Fitted \code{lmranks} object
#' 
#' Returns the variance-covariance matrix of the regression coefficients 
#' (main parameters) of a fitted \code{lmranks} object.
#' 
#' @param complete logical indicating if the full variance-covariance matrix 
#' should be returned also in case of an over-determined system where 
#' some coefficients are undefined and \code{coef(.)} contains NAs correspondingly. 
#' When \code{complete = TRUE}, \code{vcov()} is compatible with \code{coef()} also in this singular case.
#' @importFrom stats vcov
#' @export
vcov.lmranks <- function(object, complete = TRUE, ...){
  projection_residual_matrix <- get_projection_residual_matrix(object)
  X <- model.matrix(object)
  projection_residuals <- X %*% projection_residual_matrix
  original_resids <- resid(object)
  H1 <- calculate_H1(object, projection_residuals)
  
  H2 <- calculate_H2(object, projection_residuals)
  
  H1_mean <- colMeans(H1)
  H3 <- calculate_H3(object, projection_residual_matrix, H1_mean)
  
  projection_variances <- apply(projection_residuals, 2, var)
  psi <- t(t(H1 + H2 + H3) / projection_variances)
  
  sigmahat <- (t(psi) %*% psi) / (nrow(psi) ^ 2)
  colnames(sigmahat) <- names(coef(object))
  rownames(sigmahat) <- colnames(sigmahat)
  if(!complete){
    sigmahat <- sigmahat[!is.na(coef(object)),
                         !is.na(coef(object))]
  }
  return(sigmahat)
}

get_projection_residual_matrix <- function(object){
  regressor_dropped <- is.na(coef(object))
  if(any(regressor_dropped)){
    X <- model.matrix(object)[,!regressor_dropped]
    R <- qr.R(qr(X))
  } else if(is.null(object$qr)){
    R <- qr.R(qr(model.matrix(object)))
  } else {
    R <- qr.R(object$qr)
  }
  
  XTX_inv <- chol2inv(R)
  prec <- 1/diag(XTX_inv)
  out <- t(t(XTX_inv) * prec)
  
  if(!any(regressor_dropped)){
    return(out)
  }
  full_out <- matrix(NA, nrow=length(coef(object)),
                     ncol=length(coef(object)))
  full_out[!regressor_dropped, !regressor_dropped] <- out
  full_out[regressor_dropped, !regressor_dropped] <- 0
  return(full_out)
}

calculate_H1 <- function(object, projection_residuals){
  original_resids <- resid(object)
  projection_residuals * original_resids # n x p matrix  
}

#' Calculate H2 component for covariance estimation
#' 
#' Originally defined as h_2(x,y) = E[I(y,Y)-rhoI(x,X)-Wbeta)(R_X(X) - Wgamma)]
#' Estymator in matrix notation:
#' (R_X(X)-Wgamma) %*% (I_Y-rhoI_X-Wbeta) / n
#' Equal to
#' (R_X(X)-Wgamma) %*% I_Y / n - 
#' (R_X(X)-Wgamma) %*% I_X \* rho / n - 
#' (R_X(X)-Wgamma) %*% Wbeta) / n
#' @noRd
calculate_H2 <- function(object, projection_residuals){
  get_x_ineq <- get_ineq_indicator_function(object, "X")
  get_y_ineq <- get_ineq_indicator_function(object, "Y")
  H2 <- sapply(1:nobs(object), function(i){
    I_x_i <- get_x_ineq(i)
    I_Y_i <- get_y_ineq(i)
    ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
      object, I_x_i, I_Y_i)
    t(ineq_resids) %*% projection_residuals / nobs(object) 
    # ineq_resids could be delivered in batches nxp
  }) 
  t(H2)
}

calculate_H3 <- function(object, projection_residual_matrix, H1_mean){
  l <- get_and_separate_regressors(object)
  rank_column_index <- l$rank_column_index; RX <- l$RX
  if(length(rank_column_index)==0)
    return(0)
  get_x_ineq <- get_ineq_indicator_function(object, "X")
  X_projection_coef <- projection_residual_matrix[rank_column_index,]
  original_resids <- resid(object)
  weighted_X_diff <- sapply(1:nobs(object), function(i){
    ineq_X <- get_x_ineq(i)
    t(original_resids) %*% (ineq_X-RX)
  })
  H3_minus_H1_mean <- weighted_X_diff %o% 
    X_projection_coef / 
    nobs(object)
  
  t(t(H3_minus_H1_mean) + H1_mean)
}

#' Extract regressors from a model object and separate rank- from usual ones
#' 
#' @return a list with entries:
#' - RX: vector of ranks of ranked regressor. May be empty.
#' - W: matrix of non-ranked regressors.
#' - rank_column_index: which column in model.matrix corresponds to ranked regressor?
#' @noRd
get_and_separate_regressors <- function(model){
  if(length(model$rank_terms_indices) > 1) cli::cli_abort("Not implemented yet")
  rank_column_index <- get_ranked_column_index(model)
  if(length(rank_column_index) > 1) cli::cli_abort("Not implemented yet")
  if(length(rank_column_index) > 0){
    RX <- stats::model.matrix(model)[,rank_column_index]
  } else {
    RX <- integer(0)
  }
  return(list(RX=RX,
              rank_column_index=rank_column_index))
}

#' Which column in model.matrix corresponds to ranked regressor?
#' @return integer vector. Could be empty.
#' @noRd
get_ranked_column_index <- function(model){
  rank_column_index <- which(model$assign %in% model$rank_terms_indices)
}

get_ineq_indicator_function <- function(object, var_name){
  if(var_name=="X"){
    l <- get_and_separate_regressors(object)
    rank_column_index <- l$rank_column_index
    v <- l$RX
    is_ranked <- length(rank_column_index) == 1
  } else if(var_name == "Y"){
    v <- stats::model.response(stats::model.frame(object))
    is_ranked <- object$ranked_response
  } else cli::abort("Wrong var_name argument")
  
  if(!is_ranked){
    f <- function(i){
      return(v)
    }
    return(f)
  } 
  
  n_lequal_lesser <- count_lequal_lesser(v, return_inverse_ranking=TRUE)
  omega <- object$omega
  f <- function(i){
      return(get_ineq_indicator(n_lequal_lesser = n_lequal_lesser,
                                i=i,
                                omega=omega))
  }
  return(f)
}

#' Compare a certain (ith) observation against all others in (*the same*) vector 
#' of interest. The vector of interest is not passed explicitly;
#' All the required info about it is in `count_lequal_lesser`. Here the vector of interest
#' might actually be a ranked regressor or response.
#' 
#' @param count_lequal_lesser output of count_lequal_lesser function for the vector of interest.
#' To reiterate: its third element, called `inverse_ranking`, is a vector J of indices that rearrange 
#' a sorted vector of interest into its original order.
#' In other words, it is an *inverse* of sorting permutation of vector of interest
#' @param i single integer; index of the vector of interest.
#' @param omega for the definition of rank, as in irank.
#' 
#' @return A vector I of same length as vector of interest x. 
#' I[j] = 1 if x[i] < x[j];
#' I[j] = omega if x[i] == x[j]; and
#' I[j] = 0 if x[i] > x[j].
#' In other words, I[j] = I(x[i], x[j]), where the I() is the indicator function
#' as described in references
#' @noRd
get_ineq_indicator <- function(n_lequal_lesser, i, omega){
  n_lesser_or_equal <- n_lequal_lesser$n_lequal[i] # Number of observations lower or equal than ith observation
  n_lesser <- n_lequal_lesser$n_lesser[i]
  n_equal <- n_lesser_or_equal - n_lesser
  n_higher <- length(n_lequal_lesser$inverse_ranking) - n_lesser_or_equal
  # If the vector of interest were sorted (in increasing order), 
  # then the output would look like this:
  out <- rep(c(0, omega, 1), times = c(n_lesser, n_equal, n_higher))[
    # However, for generality, it has to be rearranged
               n_lequal_lesser$inverse_ranking]
  out
}

#' Here we consider inequality indicators I for a single observation (x_i, y_i). 
#' 
#' @return a numeric vector V s.t.
#' V[j] = I(y_i, Y_j) - rho*I(x_i, X_j) - beta %*% W[,j]
#' where rho and beta - respective coefficients from `model`. 
#' Note, that this function has a lot of usecases; the `model` might be the original
#' or one of the projection models (in which case, a different notation is used in references).
#' 
#' @param nonrank_predictor Numeric vector; the linear predictor part calculated from non-rank regressors.
#' @param I_X Numeric vector; Indicator whether jth X (ranked regressor) value is larger/equal/smaller 
#' than a certain (implicit & fixed for a single funtion call) ith value. Can be NULL.
#' @param I_Y Same as X, but for response. Can be simply the response, if it's not ranked. TODO: make a test for this case 
#' @noRd
replace_ranks_with_ineq_indicator_and_calculate_residuals <- function(
    model, I_X, I_Y){
  l <- get_and_separate_regressors(model)
  RX <- l$RX
  RY <- stats::model.response(stats::model.frame(model))
  rank_column_index <- l$rank_column_index
  has_ranked_regressors <- length(rank_column_index) > 0
  if(has_ranked_regressors){
    rhohat <- coef(model)[rank_column_index]
    X_diff <- (RX - I_X) * rhohat
  }
  else
    X_diff <- 0
  out <- resid(model) - RY + I_Y - X_diff
  return(out)
}

#' @importFrom stats sigma
#' @export
sigma.lmranks <- function(object, ...){
  cli::cli_abort("Not theoretically developped yet.")
}