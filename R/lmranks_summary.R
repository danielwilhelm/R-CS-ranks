# correct already:
# coerce, dummy.coef, family, formula, kappa, model.frame, model.matrix,
# nobs, print, qr, residuals, show, update, effects
# alias, hatvalues, proj, case.names, variable.names, labels

#' Summarizing fits of linear models for ranks
#' 
#' Summary method for class "\code{lmranks}". It returns theoretically valid standard
#' errors.
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
  outcome$coefficients[, 4] <- 2*pnorm(-abs(outcome$coefficients[, 3]))
  
  if(correlation)
    outcome$correlation <- cov2cor(cov_matrix)
  cli::cli_warn(c("The number of residual degrees of freedom is not correct.", 
                "Also, z-value, not t-value, since the distribution used for p-value calculation is standard normal."))
  class(outcome) <- c("summary.lmranks", class(outcome))
  outcome
}

calculate_rank_coef_std <- function(object){
  if(length(object$rank_terms_indices) > 1) cli::cli_abort("Not implemented yet")
  
  rank_column_index <- which(object$assign %in% object$rank_terms_indices)
  if(length(rank_column_index) > 1) cli::cli_abort("Not implemented yet")
  RX <- model.matrix(object)[,rank_column_index]
  W <- model.matrix(object)[,-rank_column_index]
  RY <- model.response(model.frame(object))
  I_Y <- compare(RY, omega=object$omega, increasing=TRUE, na.rm=FALSE)
  I_X <- compare(RX, omega=object$omega, increasing=TRUE, na.rm=FALSE)
  
  RX_by_Ws <- lm(RX ~ W-1)
  Wgammahat <- fitted.values(RX_by_Ws)
  nuhat <- resid(RX_by_Ws)
  gammahat <- coef(RX_by_Ws)
  
  rhohat <- coef(object)[rank_column_index]
  betahat <- coef(object)[-rank_column_index]
  epsilonhat <- resid(object)
  
  h1 <- epsilonhat * nuhat
  
  # vectorized `-` and `*` goes over rows
  h2 <- colMeans((t(I_Y - rhohat * I_X) - c(W %*% betahat)) * nuhat)
  
  h3 <- colMeans(epsilonhat * (t(I_X) - Wgammahat))
  
  sigmahat <- mean((h1+h2+h3)^2) / var(nuhat)^2
  se <- sqrt(sigmahat / length(RY))
  return(se)
}

#' @export
confint.lmranks <- function(object, parm, level = 0.95, ...){
  # As is the case with confint.lm, this method returns *marginal* CIs for coefficients
  # not simultaneous
  if(missing(parm))
    confint.default(object=object, level=level, ...)
  else
    confint.default(object=object, parm=parm, level=level, ...)
}

#' @export
vcov.lmranks <- function(object, complete = TRUE, ...){
  l <- get_and_separate_regressors(object)
  RX <- l$RX
  RY <- model.response(model.frame(object))
  if(ranked_response){
    I_Y <- compare(RY, omega=object$omega, increasing=TRUE, na.rm=FALSE)
  } else {
    I_Y <- NULL
  }
  
  if(length(object$rank_terms_indices) == 1){
    I_X <- compare(RX, omega=object$omega, increasing=TRUE, na.rm=FALSE)
  } else {
    I_X <- NULL
  }
  
  psi_sample <- sapply(1:length(coef(object)), function(i){
    if(is.na(coef(object))[i]){
      return(rep(NA, length(RY)))
    }
    proj_model <- get_projection_model(object, i)
    g_l_1 <- calculate_g_l_1(object, proj_model)
    g_l_2 <- calculate_g_l_2(object, proj_model, I_X=I_X, I_Y=I_Y)
    g_l_3 <- calculate_g_l_3(object, proj_model, I_X=I_X)
    (g_l_1 + g_l_2 + g_l_3) / var(resid(proj_model))
  })
  
  sigmahat <- (t(psi_sample) %*% psi_sample) / (nrow(psi_sample) ^ 2)
  if(!complete){
    sigmahat <- sigmahat[!is.na(coef(object)),
                         !is.na(coef(object))]
  }
  return(sigmahat)
}

get_and_separate_regressors <- function(model){
  if(length(model$rank_terms_indices) > 1) cli::cli_abort("Not implemented yet")
  rank_column_index <- which(model$assign %in% model$rank_terms_indices)
  if(length(rank_column_index) > 1) cli::cli_abort("Not implemented yet")
  if(length(rank_column_index) > 0){
    W <- model.matrix(model)[,-rank_column_index,drop=FALSE]
    RX <- model.matrix(model)[,rank_column_index]
  } else {
    W <- model.matrix(model)
    attr(W, "assign") <- NULL
    RX <- integer(0)
  }
  return(list(RX=RX,
              W=W,
              rank_column_index=rank_column_index))
}

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
  return(proj_model)
}

calculate_g_l_1 <- function(original_model, proj_model){
  epsilonhat <- resid(original_model)
  nuhat <- resid(proj_model)
  return(epsilonhat * nuhat)
}

calculate_g_l_2 <- function(original_model, proj_model, I_X, I_Y){
 ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
    original_model, I_X=I_X, I_Y=I_Y)
 nuhat <- resid(proj_model)
  return(colMeans(ineq_resids * nuhat))
}

calculate_g_l_3 <- function(original_model, proj_model, I_X){
  if(is.null(I_X)){
    return(rep(0, nobs(original_model)))
  }
  epsilonhat <- resid(original_model)
  if(proj_model$ranked_response){
  ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
    proj_model, I_Y=I_X)
  } else {
    ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
      proj_model, I_X=I_X)
  }
  return(colMeans(epsilonhat * ineq_resids))
}

replace_ranks_with_ineq_indicator_and_calculate_residuals <- function(
    model, I_X=NULL, I_Y=NULL){
  l <- get_and_separate_regressors(model)
  RX <- l$RX; W <- l$W; rank_column_index <- l$rank_column_index
  RY <- model.response(model.frame(model))
  has_ranked_regressors <- length(rank_column_index) > 0
  has_ranked_response <- model$ranked_response
  rhohat <- coef(model)[rank_column_index]
  if(length(rank_column_index) > 0)
    betahat <- coef(model)[-rank_column_index]
  else
    betahat <- coef(model)
  
  # in singular fit case, some coefficients are NA
  predictor <- as.vector(W[,!is.na(betahat), drop=FALSE] %*% betahat[!is.na(betahat)])
  if(has_ranked_regressors && has_ranked_response){
    return(t(I_Y) - t(I_X * rhohat) - predictor)
  } else if(has_ranked_response && !has_ranked_regressors){
    return(t(I_Y) - predictor)
  } else if(has_ranked_regressors && !has_ranked_response){
    return(RY - t(I_X * rhohat) - predictor)
  } else {
    return(matrix(rep(RY - predictor, times = length(RY)),
                  byrow = TRUE, nrow = length(RY)))
  }
}

#' @export
sigma.lmranks <- function(object, ...){
  cli::cli_abort("Not theoretically developped yet.")
}