#' Linear model for ranks
#' 
#' Fit a linear model with a single rank response and
#' a single rank regressor (and possibly other usual regressors).
#'
#' @param Y Numeric vector. Based on it, ranks for response will be constructed
#' and used as final response in the linear model.
#' @param X Numeric vector. Based on it, ranks for this regressor will be constructed
#' and used as final regressor in the linear model.
#' @param W Matrix. Other regressors to include in the model. Should not include an intercept column.
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
#' a single rank regressor (and possibly other usual regressors).
#' 
#' @param formula An object of class "\code{\link{formula}}": a symbolic description
#' of the model to be fitted. Exactly like the formula for linear model except that
#' rank terms, \code{r()}, can be added to specify that the linear regressor depends on ranks of regressors
#' and to specify a rank response. See Details.
#' @param subset an optional vector specifying a subset of observations to be used 
#' in the fitting process. The ranks will be calculated using full data. 
#' (See additional details about how this argument interacts with data-dependent 
#' bases in the ‘Details’ section of the \code{\link{model.frame}} documentation.)
#' @param weights currently not supported.
#' @inheritParams stats::lm
#' @param omega as in \code{\link{frank}}.
#' @param na.rm If \code{FALSE}, raises an error is any \code{NA} values in ranked regressors or response
#' are encountered. If \code{TRUE}, ranks for non-\code{NA} entries are calculated by ignoring \code{NA} values.
#' For \code{NA} values, \code{NA} ranks are returned and handled later by \code{na.action}.
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
#' As a consequence of the order, in which model.frame applies operations, \code{subset} 
#' and \code{na.action} are applied after evaluation of \code{r()}. This means, that
#' 1) the ranks will be calculated using full data. In order to calculate them on subsetted data,
#' one may subset the data outside of lm and pass it simply as new \code{data} argument.
#' 2) \code{na.action} will not handle NA values in ranked regressors. This means,
#' that they have to be handled separately by the user.
#' 
#' Currently, only models with single rank response, single rank regressor and
#' (possibly) other usual regressors is available.
#' 
#' @section Warning:
#' Wrapping \code{r()} with other functions (like \code{log(r(x))}) will not 
#' recognize correctly the mark (because it will not be caught in \code{terms(formula, specials = "r")}).
#' The ranks will be calculated correctly, but their transformation will be treated later in \code{lm} as a regular
#' regressor. This means, that the corresponding regression coefficient will be calculated correctly,
#' but the standard errors, statistics etc. will not. 
#' 
#' \code{r}, \code{.r_predict} and \code{.r_cache} are special expressions, used
#' internally to interpret \code{r} mark correctly. Do not use them in \code{formula}.
#' 
#' @return 
#' An object of class \code{lmranks}, inheriting (as well as possible) from class \code{lm}.
#' See the \code{\link{lm}} documentation for more.
#' 
#' Additionally, it has an \code{omega} entry, corresponding to \code{omega} argument,
#' and a \code{rank_terms_indices} - an integer vector with indices of entries of \code{terms.labels} attribute
#' of \code{terms(formula)}, which correspond to ranked regressors.
#' 
#' A number of methods defined for \code{lm} does not yield theoretically correct 
#' results when applied to \code{lmranks} objects; errors or warnings are raised consciously.
#' Also, the \code{df.residual} component is set to NA, since the notion of effects of freedom
#' for the rank models is not theoretically established.
#' 
#' @seealso 
#' \code{\link{lm}} for details about other arguments; \code{\link{frank}}.
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
#' # naive version with same regression coefficients, but incorrect 
#' # standard errors and statistics:
#' lm(y_frank ~ x_frank + W)
#' 
#' data(mtcars)
#' lmranks(r(mpg) ~ r(hp) + ., data = mtcars)
#' 
#' @export
lmranks <- function(formula, data, subset, 
                    weights, 
                    na.action, 
                    method = "qr", model = TRUE, x = FALSE, qr = TRUE,
                    singular.ok = TRUE, constrasts = NULL, offset = offset,
                    omega=0, na.rm=FALSE, ...){
  rank_terms_indices <- process_lmranks_formula(formula)
  original_call <- match.call() # for the final output
  lm_call <- prepare_lm_call(original_call)# to be sent to lm(); it will handle missing arguments etc
  # From this environment lm will take the definition of r()
  rank_env <- create_env_to_interpret_r_mark(omega, na.rm)
  # It will mask "r" objects from higher frames inside lm, but not modify them
  # It is also inheriting from parent.frame, so evaluations of all other expressions
  # will be taken from there
  
  # Call (evaluate) lm
  main_model <- eval(lm_call, rank_env)
  
  # Correct the output
  main_model$call <- original_call
  main_model$df.residual <- NA
  main_model$rank_terms_indices <- rank_terms_indices
  main_model$omega <- omega
  main_model$ranked_response <- TRUE
  class(main_model) <- c("lmranks", class(main_model))
  
  # Phew.
  return(main_model)
}

#' Check validity of passed formula and identify ranked terms
#' 
#' For now only formulas with one rank regressor and one rank outcome are allowed.
#' Additionally, the rank regressor cannot be part of interactions.
#'
#' @return An integer vector with indices of entries of \code{terms.labels} attribute
#' of \code{terms(formula)}, which correspond to ranked regressors.
#' 
#' @note 
#' * It allows to pass r(W), where W is a matrix. This is caught later in frank.
#' In order to catch this here, we would have to know what W is.
#' * It allows to pass r(.). This is again caught later with error ". not defined".
#' Same error occurs in lm(y ~ x + log(.), data=data). Acceptable.
#' * It will not detect func(r(expr)). 
#' 
#' @noRd
process_lmranks_formula <- function(formula){
  if(!inherits(formula, "formula")){
    cli::cli_abort(c("{.var formula} must be a {.class formula} object.",
                   "x" = "The passed {.var formula} is of {.cls {class(formula)}} class."))
  }
  formula_terms <- terms(formula, specials="r",
                         keep.order = TRUE,
                         allowDotAsName = TRUE)
  rank_variables_indices <- attr(formula_terms, "specials")[["r"]]
  response_variable_index <- attr(formula_terms, "response")
  if(length(response_variable_index) != 1){
    cli::cli_abort(c("In formula there must be exactly one ranked response and exactly one ranked regressor.",
                     "x" = "There are multiple or no responses."))
  }
  if(!response_variable_index %in% rank_variables_indices){
    cli::cli_abort(c("In formula there must be exactly one ranked response and exactly one ranked regressor.",
                   "x" = "The response is not ranked."))
  }
  if(length(rank_variables_indices) != 2){
    cli::cli_abort(c("In formula there must be exactly one ranked response and exactly one ranked regressor.",
                   "x" = "There are multiple or no ranked regressors."))
  }
  regressor_variable_index <- setdiff(rank_variables_indices, response_variable_index)
  variable_table <- attr(formula_terms, "factors")
  rank_regressor_occurances <- variable_table[regressor_variable_index,]
  occured_exactly_once <- sum(rank_regressor_occurances == 1) == 1 && sum(rank_regressor_occurances == 0) == (length(rank_regressor_occurances) - 1)
  if(!occured_exactly_once){
    cli::cli_abort("In formula, the ranked regressor must occur exactly once. No interactions are supported.")
  }
  
  rank_terms_names <- colnames(variable_table)[rank_regressor_occurances == 1]
  rearranged_formula_terms <- terms(formula, allowDotAsName = TRUE, 
                                    keep.order = FALSE) # default used later inside lm
  which(attr(rearranged_formula_terms, "term.labels") %in% rank_terms_names)
}

prepare_lm_call <- function(lm_call){
  lm_call[[1]] <- quote(stats::lm)
  lm_call$omega <- NULL              # remove local parameters
  lm_call$na.rm <- NULL
  if(!is.null(lm_call$weights))
    cli::cli_abort("{.var weights} argument is not yet supported. ")
  lm_call
}

#' Create environment to interpret lmranks formula
#' 
#' A formula in lmranks has the ranked variables (regressors and response) marked with
#' `r()`. A way to interpret this mark is needed. In R for this purpose we have
#' *environments* and *non-standard evaluation*.
#' 
#' For all intents and purposes it is enough to know, which variables have been marked
#' (that's known from `process_lmranks_formula` output) 
#' and values of ranks - that's done by treating `r()` as regular function and 
#' evaluating it in an environment containing its correct definition. 
#' 
#' One difficulty is that the `r()` must return different values depending on whether
#' the linear model is being fitted or used for prediction. In *both* cases we want
#' to use fitting (training) data. For prediction, it is stored in cache.
#' 
#' Every time the `lmranks` is called, a new environment of this kind is created 
#' and "carried along" with an `lmranks` object (accessible with `environment(model$terms)`).
#' 
#' The advantage of this solution is better reuse of `lm.fit` and `predict.lm`.
#' 
#' @return an R environment, used later to call the `lm` function in it.
#' In total, it has 3 elements: 
#' - `r`, the correct definition of `r` function, a call to frank_against 
#' with preprocessed, correct arguments.
#' - `.r_cache`, a list with values saved for possible prediction later.
#' - `.r_predict`, a logical indicating whether we are in fitting (FALSE) or prediction mode.
#' In fitting mode, input data will be used for ranking and will be saved in cache;
#' in prediction mode, the data in cache will be used.
#' 
#' Its parent environment is the parent frame of the caller of this function.
#' In the primary use case it is the environment of caller of `lmranks`.
#' In this way we ensure correct evaluation of other formula terms and variables. 
#' 
#' @seealso 
#' [H. Wickham, Advanced R, Environments chapter](https://adv-r.hadley.nz/environments.html)
#' [environment()]
#' [csranks::frank_against()], [csranks:::compare]
#' @noRd
create_env_to_interpret_r_mark <- function(omega, na.rm){
  rank_env <- new.env(parent = parent.frame(2))
  r <- function(x, increasing=FALSE) x
  body(r) <- bquote({
    predict <- get(".r_predict", envir = environment(r), inherits=FALSE)
    cache <- get(".r_cache", envir = environment(r), inherits=FALSE)
    was_na <- is.na(x)
    var_name <- paste0(as.character(substitute(x)), collapse = "")
    if(!predict){
      cache[[var_name]] <- x
      assign(".r_cache", cache, envir = environment(r))
    }
    else if(is.null(cache[[var_name]]))
      cli::cli_warn("New variable at predict time. Ranks will be calculated from scratch.")
    v <- cache[[var_name]]
    out <- rep(NA, length(was_na))
    out[!was_na] <- csranks::frank_against(x, v, increasing=increasing, omega=.(omega), na.rm=.(na.rm))
    out
  })
  environment(r) <- rank_env
  assign("r", r, envir = rank_env)
  assign(".r_cache", list(), envir = rank_env)
  assign(".r_predict", FALSE, envir = rank_env)
  return(rank_env)
}

#' @export
slotsFromS3.lmranks <- function(object){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}