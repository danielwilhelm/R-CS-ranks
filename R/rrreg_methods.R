# correct already:
# coerce, dummy.coef, family, formula, kappa, model.frame, model.matrix,
# nobs, print, qr, residuals, show, update, effects
# alias, hatvalues, proj, case.names, variable.names, labels

#' @export
logLik.lmranks <- function(object, ...){
  cli::cli_abort("This method does not return correct results.")
}

#' @export
rstandard.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' @export
simulate.lmranks <- function(object, nsim = 1, seed = NULL, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' @export
slotsFromS3.lmranks <- function(object){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' @export
anova.lmranks <- function(object, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

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
  if(length(object$rank_terms_indices) > 1) cli::cli_abort("Not implemented yet")
  
  rank_column_index <- which(object$assign %in% object$rank_terms_indices)
  if(length(rank_column_index) > 1) cli::cli_abort("Not implemented yet")
  RX <- model.matrix(object)[,rank_column_index]
  W <- model.matrix(object)[,-rank_column_index]
  RY <- model.response(model.frame(object))
  I_Y <- compare(RY, omega=object$omega, increasing=TRUE, na.rm=FALSE)
  I_X <- compare(RX, omega=object$omega, increasing=TRUE, na.rm=FALSE)
  
  rhohat <- coef(object)[rank_column_index]
  betahat <- coef(object)[-rank_column_index]
  epsilonhat <- resid(object)
  
  psi_W_sample <- sapply(1:ncol(W), function(i){
    W_minus_l <- W[,-i]
    W_l <- W[,i]
    proj_model <- lm(W_l ~ RX + W_minus_l - 1)
    
    v_l_hat <- resid(proj_model)
    delta_l_hat <- coef(proj_model)[-1]
    tau_l_hat <- coef(proj_model)[1]
    
    g_l_1 <- epsilonhat * v_l_hat
    g_l_2 <- colMeans((t(I_Y - rhohat * I_X) - c(W %*% betahat)) * v_l_hat)
    g_l_3 <- colMeans(epsilonhat * (W_l - tau_l_hat * t(I_X) - 
                                      W_minus_l %*% delta_l_hat))
    (g_l_1 + g_l_2 + g_l_3) / var(v_l_hat)
  })
  
  RX_by_Ws <- lm(RX ~ W-1)
  Wgammahat <- fitted.values(RX_by_Ws)
  nuhat <- resid(RX_by_Ws)
  gammahat <- coef(RX_by_Ws)
  
  h1 <- epsilonhat * nuhat
  
  # vectorized `-` and `*` goes over rows
  h2 <- colMeans((t(I_Y - rhohat * I_X) - c(W %*% betahat)) * nuhat)
  
  h3 <- colMeans(epsilonhat * (t(I_X) - Wgammahat))
  
  psi_sample <- matrix(nrow = length(RY), ncol=ncol(W) + 1)
  psi_sample[,rank_column_index] <- (h1 + h2 + h3)/var(nuhat)
  psi_sample[,-rank_column_index] <- t(psi_W_sample)
  
  sigmahat <- (t(psi_sample) %*% psi_sample) / (nrow(psi_sample) ^ 2)
  return(sigmahat)
}

#' @export
predict.lmranks <- function(object, newdata, se.fit = FALSE, scale = NULL, 
                            df = Inf, interval = "none", level = 0.95, type = "response",
                            terms = NULL, na.action = na.pass, pred.var = res.var/weights,
                            weights = 1, ...){
  # Disable not (yet) supported arguments 
  if(se.fit || !is.null(scale) || !is.infinite(df) || interval != "none" || 
     type != "response"){
    cli::cli_error(c("Only {.var object}, {.var newdata} and {.var na.action} arguments are currently suppoerted.",
                   "i" = "Currently, only basic prediction is supported, wihtout calculation of standard error, confidence intervals, or per-term breakdown."))
  }
  
  rank_env <- environment(object$terms)
  assign(".r_predict", TRUE, rank_env)
  out <- NextMethod()
  assign(".r_predict", FALSE, rank_env)
  out
}

#' @export
cooks.distance.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then."))
}

#' @export
dfbeta.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then."))
}

#' @export
dfbetas.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then."))
}

#' @export
influence.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then."))
}

#' @export
rstudent.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then."))
}

#' @export
extractAIC.lmranks <- function(model, ...){
  cli::cli_warn(c("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear."))
  c(NA, NA)
}

#' @export
add1.lmranks <- function(object, scope, ...){
  cli::cli_warn(c("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear."))
  NextMethod()
}

#' @export
drop1.lmranks <- function(object, scope, ...){
  cli::cli_warn(c("NAs returned.",
                 "i" = "The notion of effective degrees of freedom for rank models is not clear."))
  NextMethod()
}

#' @export
sigma.lmranks <- function(object, ...){
  cli::cli_warn(c("NA returned.",
                "i" = "The notion of effective degrees of freedom for rank models is not clear."))
  return(NA)
}

#' @export
plot.lmranks <- function(x,...){
  cli::cli_warn("Some plots are not correct.")
  NextMethod()
}

# AIC, BIC