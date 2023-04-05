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
  if(correlation || symbolic.cor){
    cli::cli_abort("{.var correlation} and {.var symbolic.cor} are not yet implemented for {.class lmranks}.")
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
  cov_matrix <- vcov(object)
  outcome$cov.unscaled <- matrix(NA, nrow = nrow(outcome$cov.unscaled),
                                 ncol = ncol(outcome$cov.unscaled))
  
  outcome$coefficients[rank_predictor_index, 2] <- sqrt(diag(cov_matrix))
  outcome$coefficients[rank_predictor_index, 3] <- 
    outcome$coefficients[rank_predictor_index, 1] / outcome$coefficients[rank_predictor_index, 2]
  outcome$coefficients[rank_predictor_index, 4] <- 2*pnorm(-abs(outcome$coefficients[rank_predictor_index, 3]))
  
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
vcov.lmranks <- function(object, ...){
  l <- get_and_separate_regressors(object)
  RX <- l$RX; W <- l$W; rank_column_index <- l$rank_column_index
  RY <- model.response(model.frame(object))
  I_Y <- compare(RY, omega=object$omega, increasing=TRUE, na.rm=FALSE)
  I_X <- compare(RX, omega=object$omega, increasing=TRUE, na.rm=FALSE)
  
  psi_W_sample <- sapply(1:ncol(W), function(l){
    W_minus_l <- W[,-l,drop=FALSE]
    W_l <- W[,l,drop=FALSE]
    proj_model <- lm(W_l ~ RX + W_minus_l - 1)
    proj_model$rank_terms_indices <- 1
    proj_model$ranked_response <- FALSE
    
    g_l_1 <- calculate_g_l_1(object, proj_model)
    g_l_2 <- calculate_g_l_2(object, proj_model, I_X=I_X, I_Y=I_Y)
    g_l_3 <- calculate_g_l_3(object, proj_model, I_X=I_X)
    (g_l_1 + g_l_2 + g_l_3) / var(resid(proj_model))
  })
  
  proj_model <- lm(RX ~ W-1)
  proj_model$rank_terms_indices <- numeric(0)
  proj_model$ranked_response <- TRUE
  h1 <- calculate_g_l_1(object, proj_model)
  h2 <- calculate_g_l_2(object, proj_model, I_X=I_X, I_Y=I_Y)
  h3 <- calculate_g_l_3(object, proj_model, I_X=I_X)
  psi_ranked_sample <- (h_1 + h_2 + h_3) / var(resid(proj_model))
  
  psi_sample <- matrix(nrow = length(RY), ncol=length(coef(model)))
  psi_sample[,rank_column_index] <- psi_ranked_sample
  psi_sample[,-rank_column_index] <- t(psi_W_sample)
  
  sigmahat <- (t(psi_sample) %*% psi_sample) / (nrow(psi_sample) ^ 2)
  return(sigmahat)
}

get_and_separate_regressors <- function(model){
  if(length(object$rank_terms_indices) > 1) cli::cli_abort("Not implemented yet")
  rank_column_index <- which(object$assign %in% object$rank_terms_indices)
  if(length(rank_column_index) > 1) cli::cli_abort("Not implemented yet")
  RX <- model.matrix(object)[,rank_column_index]
  W <- model.matrix(object)[,-rank_column_index,drop=FALSE]
  return(list(RX=RX,
              W=W,
              rank_column_index = rank_column_index))
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
  epsilonhat <- resid(original_model)
  ineq_resids <- replace_ranks_with_ineq_indicator_and_calculate_residuals(
    proj_model, I_Y=I_X)
  return(colMeans(epsilonhat * ineq_resids))
}

replace_ranks_with_ineq_indicator_and_calculate_residuals <- function(
    model, I_X=NULL, I_Y=NULL){
  l <- get_and_separate_regressors(object)
  RX <- l$RX; W <- l$W; rank_column_index <- l$rank_column_index
  RY <- model.response(model.frame(object))
  has_ranked_regressors <- length(rank_column_index) > 0
  has_ranked_response <- model$ranked_response
  rhohat <- coef(model)[rank_column_index]
  betahat <- coef(model)[-rank_column_index]
  
  if(has_ranked_regressors && has_ranked_response){
    return(t(I_Y) - t(I_X * rhohat) - W %*% betahat)
  } else if(has_ranked_response && !has_ranked_regressors){
    return(t(I_Y) - W %*% betahat)
  } else if(has_ranked_regressors && !has_ranked_response){
    return(RY - t(I_X * rhohat) - W %*% betahat)
  } else {
    return(resid(model))
  }
}

#' @export
sigma.lmranks <- function(object, ...){
  cli::cli_abort("Not theoretically developped yet.")
}