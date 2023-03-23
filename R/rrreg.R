#' Linear model for ranks
#' 
#' Fit a linear model with a single rank response and
#' a single rank covariate (and possibly other usual covariates).
#'
#' @param Y Numeric vector. Based on it, ranks for response will be constructed
#' and used as final response in the linear model.
#' @param X Numeric vector. Based on it, ranks for this covariate will be constructed
#' and used as final covariate in the linear model.
#' @param W Matrix. Other covariates to include in the model. Should not include an intercept column.
#' @param omega numeric; numeric value in [0,1], each corresponding to a different definition of the rank; default is \code{0}. See \code{\link{frank}} for details.
#' @param increasing logical; if \code{TRUE}, then large elements in \code{X} and \code{Y} receive a large rank. Otherwise, large elements receive small ranks. 
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{X} and \code{Y} (if any). 

#' @return A list with two items: 
#' \itemize{
#' \item{\code{rhohat} - the linear coefficient for rank of X.}
#' \item{\code{se} - Estimated standard error of the coefficient rhohat.}
#' }
#' 
#' @examples
#' Y <- c(3,1,2,4,5)
#' y_frank <- c(0.6, 1.0, 0.8, 0.4, 0.2)
#' X <- 1:5
#' omega <- 0.5
#' x_frank <- c(1.0, 0.8, 0.6, 0.4, 0.2)
#' W <- matrix(y_frank * 0.1 + 5 + rnorm(5, sd = 0.1), ncol = 1)
#'
#' simple_lmranks(Y, X, W)
#' 
#' @importFrom stats lm predict resid coef var
#' @export
simple_lmranks <- function(Y, X, W=NULL, omega=0, increasing=FALSE, na.rm=FALSE) {
  l <- process_simple_lmranks_args(Y, X, W, omega, increasing, na.rm)
  X <- l$X; Y <- l$Y; W <- l$W
  I_Y <- compare(Y, omega=omega, increasing=increasing, na.rm=na.rm)
  I_X <- compare(X, omega=omega, increasing=increasing, na.rm=na.rm)
  RY <- (rowSums(I_Y) + 1 - omega) / length(Y)
  RX <- (rowSums(I_X) + 1 - omega) / length(X)
  
  RX_by_Ws <- lm(RX ~ W-1)
  Wgammahat <- predict(RX_by_Ws)
  nuhat <- resid(RX_by_Ws)
  gammahat <- coef(RX_by_Ws)
  
  res <- lm(RY~RX+W-1)
  rhohat <- coef(res)[1]
  betahat <- coef(res)[-1]
  epsilonhat <- resid(res)
  
  h1 <- epsilonhat * nuhat
  
  # vectorized `-` and `*` goes over rows
  h2 <- colMeans((t(I_Y - rhohat * I_X) - c(W %*% betahat)) * nuhat)
  
  h3 <- colMeans(epsilonhat * (t(I_X) - Wgammahat))
  
  sigmahat <- mean((h1+h2+h3)^2) / var(nuhat)^2
  se <- sqrt(sigmahat / length(Y))
  names(rhohat) <- NULL
	return(list(rhohat=rhohat, se=se))
}

#' @describeIn simple_lmranks Calculate the standard error of rank regression coefficient

simple_lmranks_rho_se <- function(Y, X, W=NULL, omega=0, increasing=FALSE, na.rm=FALSE){
  lmranks_outcome <- simple_lmranks(Y, X, W=W, omega=omega, increasing=increasing, na.rm=na.rm)
  return(lmranks_outcome$se)
}

#' Linear model for ranks
#' 
#' Fit a linear model with a single rank response and
#' a single rank covariate (and possibly other usual covariates).
#' 
#' @param formula An object of class "\code{\link{formula}}": a symbolic description
#' of the model to be fitted. Exactly like the formula for linear model except that
#' rank terms, \code{r()}, can be added to specify that the linear predictor depends on ranks of predictors
#' and to specify a rank response. See Details.
#' @inheritParams stats::lm
#' @param omega as in \code{\link{frank}}.
#' 
#' @details 
#' This function is useful in case when relationship not between variables themselves, but their rank
#' (or, put differently, their ECDF value) is of interest. The variables to be ranked
#' (both dependent and independent) can be marked with \code{r()}. 
#' 
#' The \code{r()} is a private alias for \code{\link{frank}} and can accept 
#' \code{increasing} argument. However, \code{omega} argument must be specified globally 
#' (as a specification of rank definition) in the call to \code{lmranks}.
#' 
#' Currently, only models with single rank response, single rank covariate and
#' (possibly) other usual covariates is available.
#' 
#' 
#' @return 
#' An object of class \code{lmranks}, inheriting (as well as possible) from class \code{lm}.
#' See the \code{\link{lm}} documentation for more.
#' 
#' A number of methods defined for \code{lm} does not yield theoretically correct 
#' results when applied to \code{lmranks} objects; errors or warnings are raised consciously.
#' Also, the \code{df.residual} component is set to NA, since the notion of effects of freedom
#' for the rank models is not theoretically developed.
#' 
#' @seealso 
#' \code{\link{lm}} For details about other arguments; \code{\link{frank}}.
#' 
#' @examples 
#' Y <- c(3,1,2,4,5)
#' y_frank <- c(0.6, 1.0, 0.8, 0.4, 0.2)
#' X <- 1:5
#' omega <- 0.5
#' x_frank <- c(1.0, 0.8, 0.6, 0.4, 0.2)
#' W <- matrix(y_frank * 0.1 + 5 + rnorm(5, sd = 0.1), ncol = 1)
#'
#' lmranks(r(Y) ~ r(X) + W)
#' # equivalent:
#' lm(y_frank ~ x_frank + W)
#' 
#' data(mtcars)
#' lmranks(r(mpg) ~ r(hp) + ., data = mtcars)
#' 
#' @export
lmranks <- function(formula, data, subset, 
                    weights, # TODO?
                    na.action, #TODO?
                    method = "qr", model = TRUE, x = FALSE, qr = TRUE,
                    singular.ok = TRUE, constrasts = NULL, offset = offset,
                    omega=0, na.rm=FALSE, ...){
  process_lmranks_formula(formula)
  # I would gladly move that to separate function, if not for non standard evaluation
  original_call <- match.call() # for the final output
  lm_call <- match.call()       # to be sent to lm(); it will handle missing arguments etc
  lm_call[[1]] <- quote(stats::lm)
  lm_call$omega <- NULL              # remove local parameters
  lm_call$na.rm <- NULL
  if(!is.null(lm_call$weights))
    cli::cli_abort("{.var weights} argument is not yet supported. ")
  if(!is.null(lm_call$na.action))
    cli::cli_abort("{.var na.action} argument is not yet supported.")
  rank_env <- new.env(parent = parent.frame()) # From this environment lm will take the definition of r()
  r <- function(x, increasing=FALSE) x
  body(r) <- bquote({
    csranks::frank(x, increasing=increasing, omega=.(omega), na.rm=.(na.rm))
  })
  assign("r", r, envir = rank_env)
  # It will mask "r" objects from higher frames inside lm, but not modify them
  
  main_model <- eval(lm_call, rank_env)
  # Correct the output
  main_model$call <- original_call
  main_model$df.residual <- NA
  class(main_model) <- c("lmranks", class(main_model))
  # Phew.
  main_model
}

process_lmranks_formula <- function(formula){
  # for now: one rank regressor, one rank outcome, no interactions
  formula_terms <- terms(formula, specials="r", allowDotAsName = TRUE)
  rank_variables_indices <- attr(formula_terms, "specials")[["r"]]
  response_variable_index <- attr(formula_terms, "response")
  if(length(response_variable_index) != 1 && !response_variable_index %in% rank_variables_indices || length(rank_variables_indices) != 2){
    cli::cli_abort("In formula there must be exactly one ranked response and exactly one ranked regressor.")
  }
  regressor_variable_index <- setdiff(rank_variables_indices, response_variable_index)
  variable_table <- attr(formula_terms, "factors")
  rank_regressor_occurances <- variable_table[regressor_variable_index,]
  occured_exactly_once <- sum(rank_regressor_occurances == 1) == 1 && sum(rank_regressor_occurances == 0) == (length(rank_regressor_occurances) - 1)
  if(!occured_exactly_once){
    cli::cli_abort("In formula, the ranked regressor must occur exactly once. No interactions are supported.")
  }
}



