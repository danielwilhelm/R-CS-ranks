#' Regressions Involving Ranks
#' 
#' Estimation and inference for regressions involving ranks, i.e. regressions in which the dependent and/or the independent
#' variable has been transformed into ranks before running the regression.
#' 
#' @param formula An object of class "\code{\link{formula}}": a symbolic description
#' of the model to be fitted. Exactly like the formula for linear model except that
#' variables to be ranked can be indicated by \code{r()}. See Details and Examples below.
#' @param subset currently not supported.
#' @param weights currently not supported.
#' @param na.action currently not supported. User is expected to handle NA values prior to the use of this function.
#' @inheritParams stats::lm
#' @param model,y,qr logicals. If TRUE the corresponding components of the fit (the model frame, the response, the QR decomposition) are returned.
#' @param x 
#' 
#' * For `lmranks`: Logical. Should model matrix be returned?
#' * For `plot` method: An `lmranks` object.
#' @param omega real number in the interval \[0,1\] defining how ties are handled (if there are any). 
#' The value of \code{omega} is passed to \code{\link{frank}} for computation of ranks. 
#' The default is 1 so that the rank of a realized value is defined as the 
#' empirical cdf evaluated at that realized value. See Details below.
#'
#' @details 
#' This function performs estimation and inference for regressions involving ranks. Suppose there is a dependent variable \eqn{Y_i} and independent
#' variables \eqn{X_i} and \eqn{W_i}, where \eqn{X_i} is a scalar and \eqn{W_i} a vector (possibly including a constant). Instead of running a linear regression of \eqn{Y_i} on \eqn{X_i} and \eqn{W_i}, we want to first transform
#' \eqn{Y_i} and/or \eqn{X_i} into ranks. Denote by \eqn{R_i^Y} the rank of \eqn{Y_i} and \eqn{R_i^X} the rank of \eqn{X_i}. Then, a 
#' \strong{rank-rank regression}, \deqn{R_i^Y = \rho R_i^X + W_i'\beta + \varepsilon_i,} is run using the formula \code{r(Y)~r(X)+W}. Similarly, a 
#' **regression of the raw dependent variable on the ranked regressor**, 
#' \deqn{Y_i = \rho R_i^X + W_i'\beta + \varepsilon_i,} can be implemented by the formula \code{Y~r(X)+W}, and a 
#' \strong{regression of the ranked dependent variable on the raw regressors}, \deqn{R^Y_i = W_i'\beta + \varepsilon_i,} can be implemented by the formula \code{r(Y)~W}.
#' 
#' The function works, in many ways, just like \code{lm} for linear regressions. Apart from some smaller details, there are two important differences: 
#' first, in \code{lmranks}, the mark \code{r()} can be used in formulas to indicate variables to be ranked before running the regression and, second, 
#' subsequent use of \code{summary} produces a summary table with the correct standard errors, t-values and p-values (while those of the \code{lm} are not correct for
#' regressions involving ranks). See Chetverikov and Wilhelm (2023) for more details.
#' 
#' Many other aspects of the function are similar to \code{lm}. For instance, 
#' \code{.} in a formula means 'all columns not otherwise in the formula' just as in \code{lm}. An
#' intercept is included by default.
#' In a model specified as \code{r(Y)~r(X)+.}, both \code{r(X)} and \code{X} will be
#' included in the model - as it would have been in \code{lm} and, say, 
#' \code{log()} instead of \code{r()}. 
#' One can exclude \code{X} with a \code{-}, i.e. \code{r(Y)~r(X)+.-X}. See
#' \code{\link{formula}} for more about model specification.
#' 
#' The \code{r()} is a private alias for \code{\link{frank}}.
#' The \code{increasing} argument, provided at individual regressor level,
#' specifies whether the ranks should increase or decrease as regressor values increase.
#' The \code{omega} argument of \code{\link{frank}}, provided at \code{lmranks} function level,
#' specifies how ties in variables are to be handled and
#' can be supplied as argument in \code{lmranks}. For more details, see \code{\link{frank}}. 
#' By default \code{increasing} is set to \code{TRUE} and \code{omega} is set equal to \code{1},
#' which means \code{r()} computes ranks by transforming a variable through its empirical cdf.
#' 
#' 
#' Many functions defined for \code{lm} also work correctly with \code{lmranks}.
#' These include \code{\link[stats]{coef}}, \code{\link[stats]{model.frame}},
#' \code{\link[stats]{model.matrix}}, \code{\link[stats]{resid}}, 
#' \code{\link[stats]{update}} and others. 
#' On the other hand, some would return incorrect results if they treated
#' \code{lmranks} output in the same way as \code{lm}'s. The central contribution of this package
#' are \code{vcov}, \code{summary} and \code{confint} implementations using the correct asymptotic theory for regressions involving ranks.
#' 
#' See the \code{\link{lm}} documentation for more.
#'
#' @section Rank-rank regressions with clusters:
#' 
#' Sometimes, the data is divided into clusters (groups) and one is
#' interested in running rank-rank regressions separately within each cluster, where the ranks are not computed
#' within each cluster, but using all observations pooled across all clusters. Specifically, let \eqn{G_i=1,\ldots,n_G} denote 
#' a variable that indicates the cluster to which the i-th observation belongs. Then, the regression model of interest is
#' \deqn{R_i^Y = \sum_{g=1}^{n_G} 1\{G_i=g\}(\rho_g R_i^X + W_i'\beta_g) + \varepsilon_i,}
#' where \eqn{\rho_g} and \eqn{\beta_g} are now cluster-specific coefficients, but the ranks \eqn{R_i^Y} and \eqn{R_i^X} are computed as 
#' ranks among all observations \eqn{Y_i} and \eqn{X_i}, respectively. That means the rank of an observation is not computed among the other observations
#' in the same cluster, but rather among all available observations across all clusters.
#' 
#' This type of regression is implemented in the \code{lmranks} function using interaction notation: \code{r(Y)~(r(X)+W):G}. Here, the variable
#' G \strong{must} be a \code{\link{factor}}.
#' 
#' Since the theory for clustered regression mixing grouped and ungrouped (in)dependent variables is not yet developed, such a model will raise an error. 
#' Also, by default the function includes a cluster-specific intercept, i.e. \code{r(Y)~(r(X)+W):G} is internally interpreted as \code{r(Y)~(r(X)+W):G+G-1}.
#' 
#' \code{\link[stats]{contrasts}} of \code{G} must be of \code{contr.treatment} kind, 
#' which is the default.
#' 
#' @section Warning:
#' As a consequence of the order, in which \code{\link[stats]{model.frame}} applies operations, 
#' \code{subset} and \code{na.action} would be applied after evaluation of \code{r()}. 
#' That would drop some rank values from the final model frame and returned coefficients 
#' and standard errors could no longer be correct.
#' The user must handle NA values and filter the data on their own prior to usage in \code{lmranks}.
#' 
#' Wrapping \code{r()} with other functions (like \code{log(r(x))}) will not 
#' recognize correctly the mark (because it will not be caught in \code{terms(formula, specials = "r")}).
#' The ranks will be calculated correctly, but their transformation will be treated later in \code{lm} as a regular
#' regressor. This means that the corresponding regression coefficient will be calculated correctly,
#' but the standard errors, statistics etc. will not. 
#' 
#' \code{r}, \code{.r_predict} and \code{.r_cache} are special expressions, used
#' internally to interpret \code{r} mark correctly. Do not use them in \code{formula}.
#' 
#' A number of methods defined for \code{lm} do not yield theoretically correct 
#' results when applied to \code{lmranks} objects; errors or warnings are raised in those instances.
#' Also, the \code{df.residual} component is set to NA, since the notion of effects of freedom
#' for the rank models is not theoretically established (at time of 1.2 release).
#' 
#' @return 
#' An object of class \code{lmranks}, inheriting (as much as possible) from class \code{lm}.
#' 
#' Additionally, it has an \code{omega} entry, corresponding to the \code{omega} argument,
#' a \code{ranked_response} logical entry, and 
#' a \code{rank_terms_indices} - an integer vector with indices of entries of \code{terms.labels} attribute
#' of \code{terms(formula)}, which correspond to ranked regressors.
#' 
#' @references Chetverikov and Wilhelm (2023), "Inference for Rank-Rank Regressions". 
#' \href{http://arxiv.org/pdf/2310.15512}{arXiv preprint arXiv:2310.15512}
#' 
#' @seealso 
#' \code{\link{lm}} for details about other arguments; \code{\link{frank}}.
#' 
#' Generic functions \code{\link[stats]{coef}}, \code{\link[stats]{effects}}, 
#' \code{\link[stats]{residuals}},
#' \code{\link[stats]{fitted}}, \code{\link[stats]{model.frame}},
#' \code{\link[stats]{model.matrix}}, \code{\link[stats]{update}} .
#' 
#' 
#' @examples 
#' # rank-rank regression:
#' X <- rnorm(500)
#' Y <- X + rnorm(500)
#' rrfit <- lmranks(r(Y) ~ r(X))
#' summary(rrfit)
#' 
#' # naive version of the rank-rank regression:
#' RY <- frank(Y, increasing=TRUE, omega=1)
#' RX <- frank(X, increasing=TRUE, omega=1)
#' fit <- lm(RY ~ RX)
#' summary(fit)
#' # the coefficient estimates are the same as in the lmranks function, but
#' # the standard errors, t-values, p-values are incorrect
#' 
#' # support of `data` argument:
#' data(mtcars)
#' lmranks(r(mpg) ~ r(hp) + ., data = mtcars)
#' # Same as above, but use the `hp` variable only through its rank
#' lmranks(r(mpg) ~ r(hp) + . - hp, data = mtcars)
#' 
#' # rank-rank regression with clusters:
#' G <- factor(rep(LETTERS[1:4], each=nrow(mtcars) / 4))
#' lmr <- lmranks(r(mpg) ~ r(hp):G, data = mtcars)
#' summary(lmr)
#' model.matrix(lmr)
#' # Include all columns of mtcars as usual covariates:
#' lmranks(r(mpg) ~ (r(hp) + .):G, data = mtcars)
#'
#' @export
#' @importFrom stats coef lm resid predict var
#' 
lmranks <- function(formula, data, subset, 
                    weights, 
                    na.action = stats::na.fail, 
                    method = "qr", model = TRUE, x = FALSE, qr = TRUE, y = FALSE,
                    singular.ok = TRUE, contrasts = NULL, offset = offset,
                    omega=1, ...){
  # From this environment lm will take the definition of r()
  rank_env <- create_env_to_interpret_r_mark(omega)
  # It will mask "r" objects from higher frames inside lm, but not modify them
  # It is also inheriting from parent.frame, so evaluations of all other expressions
  # will be taken from there
  
  l <- process_lmranks_formula(formula, rank_env)
  rank_terms_indices <- l$rank_terms_indices; ranked_response <- l$ranked_response
  corrected_formula <- l$formula
  original_call <- match.call() # for the final output
  if(length(rank_terms_indices) == 0 && !ranked_response){
    cli::cli_warn("{.var lmranks} called with no ranked terms. Using regular lm...")
    lm_call <- prepare_lm_call(original_call, check_lm_args = FALSE)
    out <- eval(lm_call, parent.frame())
    return(out)
  }
  lm_call <- prepare_lm_call(original_call)# to be sent to lm(); it will handle missing arguments etc
  lm_call$formula <- substitute(corrected_formula)
  
  # Call (evaluate) lm
  main_model <- eval(lm_call, rank_env)
  if(method == "model.frame"){
    return(main_model)
  }
  
  # Correct the output
  main_model$rank_terms_indices <- rank_terms_indices
  check_grouping_variable(main_model)
  
  main_model$call <- original_call
  main_model$df.residual <- NA
  main_model$omega <- omega
  main_model$ranked_response <- ranked_response
  class(main_model) <- c("lmranks", class(main_model))
  
  return(main_model)
}

#' Check validity of passed formula and identify ranked terms
#' 
#' For now only formulas with (at most) one rank regressor are allowed.
#' The outcome can be either ranked or usual, continuous.
#' Additionally, the rank regressor cannot be part of interactions.
#'
#' @return Alist with three entries:
#' - `rank_terms_indices`, integer vector with indices of entries of \code{terms.labels} attribute
#' of \code{terms(formula)}, which correspond to ranked regressors.
#' This vector might be empty, which indicates no ranked regressors.
#' - `ranked_response`, logical.
#' - `formula`, corrected formula.
#' 
#' @note 
#' * It allows to pass r(W), where W is a matrix. This is caught later in frank.
#' In order to catch this here, we would have to know what W is.
#' * It allows to pass r(.). This is again caught later with error ". not defined".
#' Same error occurs in lm(y ~ x + log(.), data=data). Acceptable.
#' * It will not detect func(r(expr)). 
#' 
#' @noRd
process_lmranks_formula <- function(formula, rank_env=NULL){
  if(!inherits(formula, "formula")){
    cli::cli_abort(c("{.var formula} must be a {.class formula} object.",
                   "x" = "The passed {.var formula} is of {.cls {class(formula)}} class."))
  }
  if(is.null(rank_env))
    rank_env <- environment(formula)
  formula_terms <- stats::terms(formula, specials="r",
                         keep.order = TRUE,
                         allowDotAsName = TRUE)
  rank_variables_indices <- attr(formula_terms, "specials")[["r"]]
  response_variable_index <- attr(formula_terms, "response")
  ranked_regressor_variable_index <- setdiff(rank_variables_indices, response_variable_index)
  if(length(ranked_regressor_variable_index) > 1){
    cli::cli_abort(c("In formula there may be at most one ranked regressor.",
                   "x" = "There are multiple ranked regressors."))
  }
  is_response_ranked <- response_variable_index %in% rank_variables_indices
  if(length(ranked_regressor_variable_index) == 0){
    environment(formula) <- rank_env
    return(list(rank_terms_indices = integer(0), 
                ranked_response = is_response_ranked,
                formula=formula))
  } 
  variables_terms_table <- attr(formula_terms, "factors")
  rank_regressor_present_in_term <- variables_terms_table[ranked_regressor_variable_index,] != 0
  rank_regressor_present_in_only_1_term <- sum(rank_regressor_present_in_term) == 1
  if(!rank_regressor_present_in_only_1_term){
    cli::cli_abort(c("In formula, the ranked regressor may occur only once, as a standalone term or interacting with a global, grouping variable.",
                     "x" = "There are multiple terms involving the ranked regressor."))
  }
  is_variable_present_in_ranked_term <- variables_terms_table[,rank_regressor_present_in_term] != 0
  is_ranked_regressor_alone_in_term <- sum(is_variable_present_in_ranked_term) == 1
  if(!is_ranked_regressor_alone_in_term){
    one_interacting_var_present <- sum(is_variable_present_in_ranked_term) == 2
    if(!one_interacting_var_present){
      cli::cli_abort(c("In formula, the ranked regressor may occur only once, as a standalone term or interacting with a global, grouping variable.",
                       "x" =  "The ranked regressor interacts with multiple variables in a single term."))
    }
    interacting_var <- is_variable_present_in_ranked_term 
    interacting_var[ranked_regressor_variable_index] <- FALSE
    interacring_var_present_in_every_term <- all(variables_terms_table[interacting_var,] != 0)
    if(!interacring_var_present_in_every_term){
      cli::cli_abort("In formula, the ranked regressor may occur only once, as a standalone term or interacting with a global, grouping variable.",
                     "x" = "The grouping variable does not interact with every other term in the formula.")
    }

    # We need to replace intercept with the grouping factor
    if(attr(formula_terms, "intercept")){
      new_terms <- attr(formula_terms, "term.labels")
      if(!rownames(variables_terms_table)[interacting_var] %in% attr(formula_terms, "term.labels"))
        new_terms <- c(new_terms, rownames(variables_terms_table)[interacting_var])
      formula <- stats::reformulate(new_terms, 
                            response = formula_terms[[2]], intercept = FALSE,
                            env = rank_env)
    }
  }
  environment(formula) <- rank_env
  rank_terms_names <- colnames(variables_terms_table)[rank_regressor_present_in_term]
  rearranged_formula_terms <- stats::terms(formula, allowDotAsName = TRUE, 
                                           keep.order = FALSE) # default used later inside lm
  rank_terms_indices <- which(attr(rearranged_formula_terms, "term.labels") %in% rank_terms_names)
  return(list(rank_terms_indices = rank_terms_indices, 
              ranked_response = is_response_ranked,
              formula=formula)) 
}

prepare_lm_call <- function(lm_call, check_lm_args = TRUE){
  lm_call[[1]] <- quote(stats::lm)
  lm_call$omega <- NULL              # remove local parameters
  lm_call$na.rm <- NULL
  if(!check_lm_args){
    return(lm_call)
  }
  
  if(!is.null(lm_call$weights))
    cli::cli_abort("{.var weights} argument is not yet supported. ")
  if(!is.null(lm_call$na.action))
    cli::cli_abort("{.var na.action} argument is not yet supported. ")
  if(!is.null(lm_call$subset))
    cli::cli_abort("{.var subset} argument is not yet supported. ")
  lm_call$na.action <- str2lang("stats::na.fail")
  return(lm_call)
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
#' [csranks::frank_against()]
#' @noRd
create_env_to_interpret_r_mark <- function(omega){
  rank_env <- new.env(parent = parent.frame(2))
  r <- function(x, increasing=TRUE) x
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
    out <- csranks::frank_against(x, v, increasing=increasing, omega=.(omega), na.rm=FALSE)
    out
  })
  environment(r) <- rank_env
  assign("r", r, envir = rank_env)
  assign(".r_cache", list(), envir = rank_env)
  assign(".r_predict", FALSE, envir = rank_env)
  return(rank_env)
}

#' Return grouping variable index
#' @return integer i s.t. model.frame(object)[,i] gives the grouping variable.
#' @noRd

get_grouping_var_index <- function(object){
  rank_terms_indices <- object$rank_terms_indices
  formula_terms <- stats::terms(stats::formula(object), specials = "r")
  rank_variables_indices <- attr(formula_terms, "specials")[["r"]]
  response_variable_index <- attr(formula_terms, "response")
  regressor_variable_index <- setdiff(rank_variables_indices, response_variable_index)
  variable_table <- attr(formula_terms, "factors")
  grouping_variable_index <- setdiff(which(variable_table[,rank_terms_indices] != 0),
                                     regressor_variable_index)
  return(grouping_variable_index)
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
