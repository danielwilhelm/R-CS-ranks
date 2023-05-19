# correct already:
# coerce, dummy.coef, family, formula, kappa, model.frame, model.matrix,
# nobs, print, qr, residuals, show, update, effects
# alias, hatvalues, proj, case.names, variable.names, labels

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
  l <- get_and_separate_regressors(object)
  RX <- l$RX
  RY <- stats::model.response(stats::model.frame(object))
  if(object$ranked_response){
    irank_minmax_Y <- irank_minmax(RY, return_inverse_ranking=TRUE)
  } else {
    irank_minmax_Y <- NULL
  }
  
  if(length(object$rank_terms_indices) == 1){
    irank_minmax_X <- irank_minmax(RX, return_inverse_ranking=TRUE)
  } else {
    irank_minmax_X <- NULL
  }
  
  psi_sample <- sapply(1:length(coef(object)), function(i){
    if(is.na(coef(object))[i]){
      return(rep(NA, length(RY)))
    }
    proj_model <- get_projection_model(object, i)
    g_l_1 <- calculate_g_l_1(object, proj_model)
    g_l_2 <- calculate_g_l_2(object, proj_model, irank_minmax_X=irank_minmax_X, irank_minmax_Y=irank_minmax_Y)
    g_l_3 <- calculate_g_l_3(object, proj_model, irank_minmax_X=irank_minmax_X)
    (g_l_1 + g_l_2 + g_l_3) / var(resid(proj_model))
  })
  
  sigmahat <- (t(psi_sample) %*% psi_sample) / (nrow(psi_sample) ^ 2)
  colnames(sigmahat) <- names(coef(object))
  rownames(sigmahat) <- colnames(sigmahat)
  if(!complete){
    sigmahat <- sigmahat[!is.na(coef(object)),
                         !is.na(coef(object))]
  }
  return(sigmahat)
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
    W <- stats::model.matrix(model)[,-rank_column_index,drop=FALSE]
    RX <- stats::model.matrix(model)[,rank_column_index]
  } else {
    W <- stats::model.matrix(model)
    attr(W, "assign") <- NULL
    RX <- integer(0)
  }
  return(list(RX=RX,
              W=W,
              rank_column_index=rank_column_index))
}

get_ranked_column_index <- function(model){
  rank_column_index <- which(model$assign %in% model$rank_terms_indices)
}

#' Construct a *projection* onto a selected regressor
#' i.e. a linear model with one of regressors as response variable
#' in terms of the rest of the regressors (but not the original response)
#' Used heavily in vcov
#' @noRd
get_projection_model <- function(original_model, projected_regressor_index){
  l <- get_and_separate_regressors(original_model)
  RX <- l$RX; W <- l$W; rank_column_index <- l$rank_column_index
  if(length(rank_column_index) == 0){
    l <- projected_regressor_index
    W_minus_l <- W[,-l,drop=FALSE]
    W_l <- W[,l]
    if(ncol(W_minus_l) > 0)
      proj_model <- lm(W_l ~ W_minus_l - 1)
    else
      cli::cli_abort("Not theoretically developped yet.")
    proj_model$rank_terms_indices <- integer(0)
    proj_model$ranked_response <- FALSE
    proj_model$omega <- original_model$omega
    return(proj_model)
  }
  if(projected_regressor_index %in% rank_column_index){
    proj_model <- lm(RX ~ W-1)
    proj_model$rank_terms_indices <- numeric(0)
    proj_model$ranked_response <- TRUE  
  } else {
    l <- projected_regressor_index - sum(projected_regressor_index > rank_column_index)
    W_minus_l <- W[,-l,drop=FALSE]
    W_l <- W[,l]
    if(ncol(W_minus_l) > 0)
      proj_model <- lm(W_l ~ RX + W_minus_l - 1)
    else
      proj_model <- lm(W_l ~ RX - 1)
    proj_model$rank_terms_indices <- 1
    proj_model$ranked_response <- FALSE
  }
  proj_model$omega <- original_model$omega
  return(proj_model)
}

calculate_g_l_1 <- function(original_model, proj_model){
  epsilonhat <- resid(original_model)
  nuhat <- resid(proj_model)
  return(epsilonhat * nuhat)
}

calculate_g_l_2 <- function(original_model, proj_model, irank_minmax_X, irank_minmax_Y){
 nuhat <- resid(proj_model)
  return(calculate_weighted_ineq_resid_means(
    original_model, weights=nuhat, irank_minmax_X=irank_minmax_X, irank_minmax_Y=irank_minmax_Y
  ))
}

calculate_g_l_3 <- function(original_model, proj_model, irank_minmax_X){
  if(is.null(irank_minmax_X)){
    return(rep(0, stats::nobs(original_model)))
  }
  epsilonhat <- resid(original_model)
  if(proj_model$ranked_response){
    return(calculate_weighted_ineq_resid_means(
      proj_model, weights=epsilonhat, irank_minmax_Y=irank_minmax_X))
  } else {
    return(calculate_weighted_ineq_resid_means(
      proj_model, weights=epsilonhat, irank_minmax_X=irank_minmax_X))
  }
}

#' Sometimes, instead of regular residuals, we are interested in residuals in a model,
#' where we replace ranks with indicator of relation (>=).
#' So instead of sth like RY = a*RX + b*W + c
#' we will have I_Y = a*I_X + b*W + c
#' 
#' Easier done in paper than here, hence we have a separate method that does just that.
#' It has several use cases (ranked response and regressor / only response / only regressor). 
#'
#' @param weights Numeric vector. Usually the ordinary residuals of "the other model"
#' (if `model` is the main model, the "other" is one of the projected models and vice versa).
#' @param irank_minmax_X output of irank_minmax function for the ranked regressor. Could be NULL.
#' @param irank_minmax_Y Same as above, but for the response. Also could be NULL.
#' @return A vector of length equal to number of observations in `model`.
#' Its ith element is an estimate of expected weighted difference between prediction and response
#' In a situation when we use the `model` to predict inequality indicator I(y_i,Y) instead of response's rank
#' Using I(x_i, X) instead of a rank of regressor X.
#' @noRd
calculate_weighted_ineq_resid_means <- function(
    model, weights, irank_minmax_X=NULL, irank_minmax_Y=NULL){
  predictor <- extract_nonrank_predictor(model)
  has_ranked_regressors <- !is.null(irank_minmax_X)
  has_ranked_response <- !is.null(irank_minmax_Y)
  if(has_ranked_regressors && has_ranked_response){
    return(sapply(1:(stats::nobs(model)), function(i){
      i_X <- get_ineq_indicator(irank_minmax_X, i, model$omega)
      Y <- get_ineq_indicator(irank_minmax_Y, i, model$omega)
      ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
        model, predictor, i_X, Y)
      mean(ineq_resids * weights)
    }))
  } else if(has_ranked_regressors && !has_ranked_response){
    Y <- stats::model.response(stats::model.frame(model))
    return(sapply(1:(stats::nobs(model)), function(i){
      i_X <- get_ineq_indicator(irank_minmax_X, i, model$omega)
      ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
        model, predictor, i_X, Y)
      mean(ineq_resids * weights)
    }))
  } else if(!has_ranked_regressors && has_ranked_response){
    i_X <- rep(0, stats::nobs(model))
    return(sapply(1:(stats::nobs(model)), function(i){
      Y <- get_ineq_indicator(irank_minmax_Y, i, model$omega)
      ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
        model, predictor, i_X, Y)
      mean(ineq_resids * weights)
    }))
  } else {
    i_X <- rep(0, stats::nobs(model))
    Y <- stats::model.response(stats::model.frame(model))
    ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
      model, predictor, i_X, Y)
    return(rep(mean(ineq_resids * weights),
               stats::nobs(model))) 
  }
}

#' @return Numeric vector; the linear predictor part calculated from non-rank regressors.
#' @noRd
extract_nonrank_predictor <- function(model){
  l <- get_and_separate_regressors(model)
  W <- l$W; rank_column_index <- l$rank_column_index
  has_ranked_regressors <- length(rank_column_index) > 0
  if(has_ranked_regressors)
    betahat <- coef(model)[-rank_column_index]
  else
    betahat <- coef(model)
  
  # in singular fit case, some coefficients are NA
  predictor <- as.vector(W[,!is.na(betahat), drop=FALSE] %*% betahat[!is.na(betahat)])
  return(predictor)
}

#' Compare a certain (ith) observation against all others in (*the same*) vector 
#' of interest. The vector of interest is not passed explicitly;
#' All the required info about it is in `irank_min_max`. Here the vector of interest
#' might actually be a ranked regressor or response.
#' 
#' @param irank_min_max output of irank_minmax function for the vector of interest.
#' @param i single integer; index of the vector of interest.
#' @param omega for the definition of rank, as in irank.
#' 
#' @return A vector v of same length as vector of interest x. 
#' v[j] = 1 if x[i] < x[j];
#' v[j] = omega if x[i] == x[j]; and
#' v[j] = 0 if x[i] > x[j].
#' @noRd
get_ineq_indicator <- function(irank_min_max, i, omega){
  # Note: we actually do not need to know the exact ordering in the vector of interest
  n_lower_or_equal <- irank_min_max[i,1]
  n_lower <- irank_min_max[i,2]
  n_equal <- n_lower_or_equal - n_lower
  n_higher <- nrow(irank_min_max) - n_lower_or_equal
  out <- rep(c(0, omega, 1), times = c(n_lower, n_equal, n_higher))[irank_min_max[,3]]
  out
}

#' Here we consider inequality indicators I for a single observation (x_i, y_i). 
#' 
#' @return a numeric vector V s.t.
#' V[j] = I(y_i, Y_j) - rho*I(x_i, X_j) - beta %*% W[,j]
#' 
#' @param nonrank_predictor Numeric vector; the linear predictor part calculated from non-rank regressors.
#' @param I_X Numeric vector; Indicator whether jth X (ranked regressor) value is larger/equal/smaller 
#' than a certain (implicit & fixed for a single funtion call) ith value. Can be NULL.
#' @param I_Y Same as X, but for response. Can be simply the response, if it's not ranked. TODO: make a test for this case 
#' @noRd
replace_ranks_with_ineq_indicator_and_calculate_residuals <- function(
    model, nonrank_predictor, I_X, I_Y){
  rank_column_index <- get_ranked_column_index(model)
  has_ranked_regressors <- length(rank_column_index) > 0
  if(has_ranked_regressors){
    rhohat <- coef(model)[rank_column_index]
    rank_predictor <- I_X * rhohat
  }
  else
    rank_predictor <- 0
  return(I_Y - rank_predictor - nonrank_predictor)
}

#' @importFrom stats sigma
#' @export
sigma.lmranks <- function(object, ...){
  cli::cli_abort("Not theoretically developped yet.")
}