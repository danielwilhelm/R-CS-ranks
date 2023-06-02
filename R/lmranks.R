#' Linear model for ranks
#' 
#' Fit a linear model with a single rank response and
#' a single rank regressor (and possibly other usual regressors).
#' 
#' @param formula An object of class "\code{\link{formula}}": a symbolic description
#' of the model to be fitted. Exactly like the formula for linear model except that
#' rank terms, \code{r()}, can be added to specify that the linear regressor depends on ranks of regressors
#' and to specify a rank response. See Details.
#' @param subset currently not supported.
#' @param weights currently not supported.
#' @param na.action currently not supported. User is expected to handle NA values on his own.
#' @inheritParams stats::lm
#' @param model,y,qr logicals. If TRUE the corresponding components of the fit (the model frame, the response, the QR decomposition) are returned.
#' @param x \itemize{
#' \item{For \code{lmranks}: }{Logical. Should model matrix be returned?}
#' \item{For \code{plot} method: }{An \code{lmranks} object.}
#' }
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
#' The \code{r()} is a private alias for \code{\link{frank}} with fixed
#' \code{increasing} argument. \code{omega} argument may be specified globally 
#' (as a specification of rank definition) in the call to \code{lmranks}.
#' 
#' As a consequence of the order, in which model.frame applies operations, \code{subset} 
#' and \code{na.action} would be applied after evaluation of \code{r()}. 
#' In such a case, returned coefficients and standard errors might no longer be correct.
#' The user must handle the NA values and may filter the data on his own.
#' 
#' Currently, only models at most one rank regressor are available. The single 
#' response might be either ranked or continuous.
#' 
#' Many functions defined for \code{lm} also work correctly with \code{lmranks}.
#' This includes \code{\link[stats]{coef}}, \code{\link[stats]{model.frame}},
#' \code{\link[stats]{model.matrix}}, \code{\link[stats]{resid}} and \code{\link[stats]{update}}. 
#' On the other hand, some would not return correct results. 
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
#' \code{\link{summary.lmranks}}
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
#' @importFrom stats coef lm resid predict var
#' 
lmranks <- function(formula, data, subset, 
                    weights, 
                    na.action = stats::na.fail, 
                    method = "qr", model = TRUE, x = FALSE, qr = TRUE, y = FALSE,
                    singular.ok = TRUE, contrasts = NULL, offset = offset,
                    omega=1, na.rm=FALSE, ...){
  l <- process_lmranks_formula(formula)
  rank_terms_indices <- l$rank_terms_indices; ranked_response <- l$ranked_response
  original_call <- match.call() # for the final output
  if(length(rank_terms_indices) == 0 && !ranked_response){
    cli::cli_warn("{.var lmranks} called with no ranked terms. Using regular lm...")
    lm_call <- prepare_lm_call(original_call, check_weights = FALSE)
    out <- eval(lm_call, parent.frame())
    return(out)
  }
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
  main_model$ranked_response <- ranked_response
  class(main_model) <- c("lmranks", class(main_model))
  
  # Phew.
  return(main_model)
}

#' Check validity of passed formula and identify ranked terms
#' 
#' For now only formulas with (at most) one rank regressor are allowed.
#' The outcome can be either ranked or usual, continuous.
#' Additionally, the rank regressor cannot be part of interactions.
#'
#' @return Alist with two entries:
#' - `rank_terms_indices`, integer vector with indices of entries of \code{terms.labels} attribute
#' of \code{terms(formula)}, which correspond to ranked regressors.
#' This vector might be empty, which indicates no ranked regressors.
#' - `ranked_response`, logical.
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
  formula_terms <- stats::terms(formula, specials="r",
                         keep.order = TRUE,
                         allowDotAsName = TRUE)
  rank_variables_indices <- attr(formula_terms, "specials")[["r"]]
  response_variable_index <- attr(formula_terms, "response")
  regressor_variable_index <- setdiff(rank_variables_indices, response_variable_index)
  if(length(regressor_variable_index) > 1){
    cli::cli_abort(c("In formula there may be at most one ranked regressor.",
                   "x" = "There are multiple ranked regressors."))
  }
  ranked_response <- response_variable_index %in% rank_variables_indices
  if(length(regressor_variable_index) == 0){
    return(list(rank_terms_indices = integer(0), 
                ranked_response = ranked_response))
  } else {
    variable_table <- attr(formula_terms, "factors")
    rank_regressor_occurances <- variable_table[regressor_variable_index,]
    occured_exactly_once <- sum(rank_regressor_occurances == 1) == 1 && sum(rank_regressor_occurances == 0) == (length(rank_regressor_occurances) - 1)
    if(!occured_exactly_once){
      cli::cli_abort("In formula, the ranked regressor may occur only once. No interactions are supported.")
    }
    rank_terms_names <- colnames(variable_table)[rank_regressor_occurances == 1]
    rearranged_formula_terms <- stats::terms(formula, allowDotAsName = TRUE, 
                                      keep.order = FALSE) # default used later inside lm
    rank_terms_indices <- which(attr(rearranged_formula_terms, "term.labels") %in% rank_terms_names)
    return(list(rank_terms_indices = rank_terms_indices, 
                ranked_response = ranked_response))
  }
}

prepare_lm_call <- function(lm_call, check_weights = TRUE){
  lm_call[[1]] <- quote(stats::lm)
  lm_call$omega <- NULL              # remove local parameters
  lm_call$na.rm <- NULL
  if(check_weights && !is.null(lm_call$weights))
    cli::cli_abort("{.var weights} argument is not yet supported. ")
  if(!is.null(lm_call$na.action))
    cli::cli_abort("{.var na.action} argument is not yet supported. ")
  if(!is.null(lm_call$subset))
    cli::cli_abort("{.var subset} argument is not yet supported. ")
  lm_call$na.action <- str2lang("stats::na.fail")
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
  r <- function(x) x
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
    out[!was_na] <- csranks::frank_against(x, v, increasing=TRUE, omega=.(omega), na.rm=.(na.rm))
    out
  })
  environment(r) <- rank_env
  assign("r", r, envir = rank_env)
  assign(".r_cache", list(), envir = rank_env)
  assign(".r_predict", FALSE, envir = rank_env)
  return(rank_env)
}

# Inherited `lm` methods:
# coerce, dummy.coef, family, formula, kappa, model.frame, model.matrix,
# nobs, print, qr, residuals, show, update, effects
# alias, hatvalues, proj, case.names, variable.names, labels

slotsFromS3.lmranks <- function(object){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' @describeIn lmranks Plot diagnostics for an \code{lmranks} object
#' 
#' Displays plots useful for assessing quality of model fit. Currently, only one
#' plot is available, which plots fitted values against residuals (for homoscedacity check).
#' 
#' @param which As in \code{\link{plot.lm}}. Currently only no.1 is available.
#' @export
plot.lmranks <- function(x,which = 1,...){
  if(length(which) != 1 || which != 1)
    cli::cli_abort('For now, only basic "residuals against fitted" plot is supported.')
  NextMethod(which = which)
}