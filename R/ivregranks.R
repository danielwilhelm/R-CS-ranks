#' Instrumental-Variable Regression by 2SLS Estimation Involving
#' Ranks
#'
#' Fit instrumental-variable regression involving ranks by two-stage least
#' squares (2SLS). This is equivalent to direct instrumental-variables
#' estimation when the number of instruments is equal to the number of
#' regressors.
#'
#' Regressors and instruments for \code{ivregranks} are most easily specified
#' in a formula with two parts on the right-hand side, e.g.,
#' \code{r(y) ~ x1 + r(x2) | r(z1) + z2 + z3}, where \code{x1} and \code{r(x2)}
#' are the explanatory variables and \code{r(z1)}, \code{z2}, and \code{z3} are
#' the instrumental variables. Note that exogenous regressors have to be
#' included as instruments for themselves.
#'
#' For example, if there is
#' one exogenous regressor \code{ex} and one endogenous regressor \code{r(en)}
#' with instrument \code{r(in)}, the appropriate formula would be \code{r(y) ~
#' r(en) + ex | r(in) + ex}. Alternatively, a formula with three parts on the
#' right-hand side can also be used: \code{r(y) ~ ex | r(en) | r(in)}.
#' The latter is typically more convenient, if there is a large number of
#' exogenous regressors.
#'
#' Moreover, two further equivalent specification strategies are possible that
#' are typically less convenient compared to the strategies above. One option
#' is to use an update formula with a \code{.} in the second part of the formula
#' is used: \code{r(y) ~ r(en) + ex | . - r(en) + r(in)}. Another option is to
#' use a separate formula for the instruments (only for backward compatibility
#' with earlier versions):
#' \code{formula = r(y) ~ r(en) + ex, instruments = ~ r(in) + ex}.
#'
#' Internally, all specifications are converted to the version with two parts
#' on the right-hand side.
#'
#' @param formula formula specification(s) of the regression
#' relationship and the instruments. \code{formula} has three parts as in \code{r(y) ~ x1 + r(x2) | r(z1) + z2 +
#' z3}.
#' @param data an optional data frame containing the variables in the model.
#' By default the variables are taken from the environment of the
#' \code{formula}.
#' @param subset currently not supported.
#' @param na.action currently not supported.
#' @param weights currently not supported.
#' @param offset an optional offset that can be used to specify an a priori
#' known component to be included during fitting.
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#' \code{\link[stats:model.matrix]{model.matrix.default}}.
#' @param model,x,y logicals.  If \code{TRUE} the corresponding components of
#' the fit (the model frame, the model matrices, the response) are returned.
#' These components are necessary for computing regression diagnostics.
#' @param method the method used to fit the stage 1 and 2 regression:
#' \code{"OLS"} for traditional 2SLS regression (the default and only option).
#' @param omega real number in the interval \[0,1\] defining how ties are
#' handled (if there are any).
#' @param \dots further arguments passed to \code{\link[ivreg]{ivreg.fit}}.
#'
#' @return \code{ivregranks} returns an object of class \code{"ivregranks"} that
#' inherits as much as possible from class \code{\link[ivreg]{ivreg}},
#' with the following additional components:
#' \item{rank_terms_indices}{an integer vector with indices of entries of
#' \code{terms.labels} attribute of \code{terms(formula)} for the outcome
#' equation which correspond to ranked regressors.}
#' \item{rank_instruments_indices}{an integer vector with indices of entries
#' of the ranked instrumental variables.}
#' \item{ranked_response}{a logical entry.}
#' \item{omega}{an entry corresponding to the \code{omega} argument.}
#'
#' @references Chetverikov and Wilhelm (2023), "Inference for Rank-Rank Regressions".
#' \href{http://arxiv.org/pdf/2310.15512}{arXiv preprint arXiv:2310.15512}
#'
#' @seealso \code{\link[ivreg]{ivreg.fit}}, \code{\link[csranks]{lmranks}}
#'
#' Generic functions \code{\link[stats]{coef}}, \code{\link[stats]{residuals}},
#' \code{\link[stats]{fitted}}, \code{\link[stats]{model.frame}},
#' \code{\link[stats]{model.matrix}}, \code{\link[stats]{update}} .
#'
#' @examples
#' # rank-rank regression:
#' Z <- rnorm(500)
#' X <- Z + rnorm(500)
#' Y <- X + rnorm(500)
#' rrfit <- ivregranks(r(Y) ~ r(X) | r(Z))
#' summary(rrfit)
#'
#' # naive version of the rank-rank regression:
#' RZ <- frank(Z, increasing = TRUE, omega = 1)
#' RX <- frank(X, increasing = TRUE, omega = 1)
#' RY <- frank(Y, increasing = TRUE, omega = 1)
#' fit <- ivreg::ivreg(RY ~ RX | RZ)
#' summary(fit)
#' # the coefficient estimates are the same as in the ivregranks function, but
#' # the standard errors, t-values, p-values are incorrect.
#'
#' # support for `data` argument:
#' ivr <- ivregranks(r(mpg) ~ r(hp) + cyl | r(disp) + cyl, data = mtcars)
#' summary(ivr)
#'
#' @export
ivregranks <- function(formula, data, subset, na.action, weights,
                       offset, contrasts = NULL, model = TRUE, y = TRUE,
                       x = FALSE, method = "OLS", omega = 1,
                       ...) {
  method <- rlang::try_fetch(match.arg(method), error = function(e) {
    cli::cli_abort(c("Only the OLS method is supported.",
      "x" = "Estimation method {method} is not supported."
    ))
  })
  rank_env <- create_env_to_interpret_r_mark(omega)
  l <- process_ivregranks_formula(formula,
    data = if (missing(data)) NULL else data,
    rank_env = rank_env
  )
  rank_instruments_indices <- l$rank_instruments_indices
  rank_terms_indices <- l$rank_terms_indices
  ranked_response <- l$ranked_response
  formula <- l$formula
  original_call <- match.call()

  if (length(rank_terms_indices) == 0 &&
    length(rank_instruments_indices) == 0 &&
    !ranked_response) {
    cli::cli_warn("{.var ivregranks} called with no ranked terms.
      Using regular ivreg...")
    ivreg_call <- prepare_ivreg_call(original_call, check_ivreg_args = FALSE)
    out <- eval(ivreg_call, parent.frame())
    return(out)
  }
  ivreg_call <- prepare_ivreg_call(original_call)
  ivreg_call$formula <- substitute(formula)

  main_model <- eval(ivreg_call, rank_env)

  corrected_formula <- Formula::as.Formula(main_model$formula)
  endogenous_target <- names(main_model$endogenous)
  formula_update <- stats::as.formula(paste0(endogenous_target, "~."))
  formula_fs <- stats::update(main_model$terms$instruments, formula_update) # nolint: object_usage_linter

  # needed to correctly compute the vcov for the first-stage
  # TODO: call lm() if nothing is ranked
  lmranks_call_fs <- prepare_lmranks_call_for_stage_1(original_call)
  lmranks_call_fs$formula <- substitute(formula_fs)

  object_fs <- eval(lmranks_call_fs, rank_env)

  main_model$formula <- corrected_formula
  main_model$rank_terms_indices <- rank_terms_indices
  main_model$rank_instruments_indices <- rank_instruments_indices
  main_model$call <- original_call
  main_model$df.residual <- NA
  main_model$omega <- omega
  main_model$ranked_response <- ranked_response
  main_model$object_fs <- object_fs
  class(main_model) <- c("ivregranks", class(main_model))

  return(main_model)
}

#' @param formula a Formula::Formula in canonical form.
#' @noRd
get_exogenous_variable_names <- function(formula, data) {
  stage_2_terms <-stats::terms(formula, lhs = 0, rhs = 1, data = data)
  stage_2_variables <- rownames(attr(stage_2_terms, "factors"))
  stage_1_terms <- stats::terms(formula, lhs = 0, rhs = 2, data = data)
  stage_1_variables <- rownames(attr(stage_1_terms, "factors"))

  return(setdiff(stage_2_variables, stage_1_variables))
}

get_instrument_variable_names <- function(formula, data) {
  stage_2_terms <- stats::terms(formula, lhs = 0, rhs = 1, data = data)
  stage_2_variables <- rownames(attr(stage_2_terms, "factors"))
  stage_1_terms <- stats::terms(formula, lhs = 0, rhs = 2, data = data)
  stage_1_variables <- rownames(attr(stage_1_terms, "factors"))

  return(setdiff(stage_1_variables, stage_2_variables))
}

#' Check validity of passed formula and identify ranked terms
#'
#' For now only formulas with (at most) one rank regressor and one
#' ranked instrument are allowed.
#' The outcome can be either ranked or not.
#' Additionally, the rank regressor/instrument cannot be part of interactions.
#'
#' @return A list with four entries:
#' - `rank_terms_indices`, integer vector with indices of entries of
#' \code{terms.labels} attribute of \code{terms(formula)}, which correspond to
#' ranked regressors for the outcome equation.
#' This vector might be empty, which indicates no ranked regressors.
#' - `rank_instruments_indices`, integer vector with indices of entries of the
#' ranked instrumental variables
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
process_ivregranks_formula <- function(formula, data, rank_env = NULL) {
  assert_is_formula(formula)
  if (is.null(rank_env)) {
    rank_env <- environment(formula)
  }

  canonical_formula <- convert_formula_to_canonical_form(formula)

  # TODO: call lmranks ourselves?
  if (length(canonical_formula)[2] == 1) {
    cli::cli_abort(
      c("{.var formula} must at have least two/at most three regressor parts",
        "x" = "The passed {.var formula} has a single part regressor",
        "i" = "Use lmranks."
      )
    )
  }
  formula <- canonical_formula

  # Process structural equation (stage 2)
  # One instrument for 1 endogenous variable is allowed.
  # The instrument, as well as endogenous variable and final target, must be ranked.
  structural_formula <- formula(formula, lhs = 1, rhs = 1)
  # Here we check whether the final response is ranked

  l1 <- ivregranks_processor(structural_formula, "Error while processing structural formula (for stage 2): ")

  instruments_formula <- get_instruments_formula(formula, data)
  instrument_variable <- get_instrument_variable_names(formula, data)

  # Here we check whether the endogeneous variable is ranked
  l2 <- ivregranks_processor(instruments_formula, "Error while processing instruments formula (for stage 1): ", instrument_variable)

  rank_terms_indices <- l1$rank_terms_indices
  rank_instruments_indices <- l2$rank_terms_indices

  environment(formula) <- rank_env

  return(list(
    rank_terms_indices = rank_terms_indices,
    rank_instruments_indices = rank_instruments_indices,
    ranked_response = l1$ranked_response,
    formula = formula
  ))
}

#' @noRd
has_dot <- function(formula) {
  inherits(
    try(stats::terms(formula), silent = TRUE),
    "try-error"
  )
}

#' @noRd
convert_formula_to_canonical_form <- function(formula) {
  formula <- Formula::as.Formula(formula)
  # Following logic is a copy-paste from iverg.R
  if (length(formula)[2L] == 3L) {
    canonical_formula <- Formula::as.Formula(
      formula(formula, rhs = c(2L, 1L), collapse = TRUE),
      formula(formula, lhs = 0L, rhs = c(3L, 1L), collapse = TRUE)
    )
  } else {
    canonical_formula <- formula
  }

  # More graceful than stopifnot
  if (length(canonical_formula)[1] != 1 || length(canonical_formula)[2] >= 3) {
    cli::cli_abort(c("{.var formula} must contain a single outcome and at least
      an instrument part.",
      "x" = "The passed {.var formula} has either a multi-part response
      or more than three-part regressors."
    ))
  }
  formula <- canonical_formula

  ## try to handle dots in formula
  if (has_dot(formula)) {
    f1 <- formula(formula, rhs = 1L)
    f2 <- formula(formula, lhs = 0L, rhs = 2L)
    if (!has_dot(f1) && has_dot(f2)) {
      formula <- Formula::as.Formula(
        f1,
        stats::update(formula(formula, lhs = 0L, rhs = 1L), f2)
      )
    }
  }

  return(formula)
}

get_instruments_formula <- function(formula, data) {
  instruments_formula <- formula(formula, rhs = 2)
  endogeneous_variables <- get_exogenous_variable_names(formula, data)
  formula_update <- stats::as.formula(paste0(paste(endogeneous_variables, collapse = "+"), "~."))
  stats::update(instruments_formula, formula_update)
}

#' @noRd
internal_processor <- function(formula, expected_ranked_variable_str = NULL) {
  assert_is_formula(formula)

  parsed_formula <- parse_lmranks_formula(formula)

  if (has_no_ranked_terms_nor_response(parsed_formula)) {
    return(list(
      rank_terms_indices = integer(0),
      ranked_response = FALSE,
      formula = formula
    ))
  }

  assert_has_ranked_response(parsed_formula)
  if (is.null(expected_ranked_variable_str)) {
    assert_has_exactly_one_ranked_regressor(parsed_formula)
  } else {
    assert_has_exactly_one_ranked_regressor_equal_to(parsed_formula, expected_ranked_variable_str)
  }

  assert_has_at_most_one_term_with_ranked_regressor(parsed_formula)
  assert_ranked_regressor_does_not_interact_with_other_variables(parsed_formula)

  # Find rank terms for the final output
  rank_terms_names <- parsed_formula[["ranked_terms_labels"]]
  rank_terms_indices <- get_rank_terms_indices_after_reordering(formula, rank_terms_names)

  return(list(
    rank_terms_indices = rank_terms_indices,
    ranked_response = parsed_formula[["is_response_ranked"]],
    formula = formula
  ))
}

#' @noRd
ivregranks_processor <- function(formula, error_message, expected_ranked_variable_str = NULL) {
  error_message <- glue::glue(error_message, deparse(formula))
  rlang::try_fetch(
    internal_processor(formula, expected_ranked_variable_str),
    error = function(e) {
      cli::cli_abort(error_message, parent = e)
    }
  )
}

#' @noRd
prepare_ivreg_call <- function(ivreg_call, check_ivreg_args = TRUE) {
  ivreg_call[[1]] <- quote(ivreg::ivreg)
  ivreg_call$omega <- NULL
  ivreg_call$na.rm <- NULL

  if (!check_ivreg_args) {
    return(ivreg_call)
  }

  if (!is.null(ivreg_call$weights)) {
    cli::cli_abort("{.var weights} argument is not yet supported.")
  }
  if (!is.null(ivreg_call$na.action)) {
    cli::cli_abort("{.var na.action} argument is not yet supported.")
  }
  if (!is.null(ivreg_call$subset)) {
    cli::cli_abort("{.var subset} argument is not yet supported.")
  }

  ivreg_call$na.action <- str2lang("stats::na.fail")

  return(ivreg_call)
}

#' @noRd
prepare_lmranks_call_for_stage_1 <- function(ivreg_call) {
  # no checking, since it will be checked downstream by lmranks
  ivreg_call[[1]] <- quote(csranks::lmranks)
  ivreg_call$method <- NULL
  ivreg_call
}

#' @describeIn ivregranks Plot diagnostics for an \code{ivregranks} object
#'
#' Displays plots useful for assessing quality of model fit. Currently, only one
#' plot is available, which plots fitted values against residuals
#' (for homoscedacity check).
#'
#' @param which As in \code{\link[ivreg]{plot.ivreg}}. Currently only no. 1 is
#' available.
#' @export
plot.ivregranks <- function(x, which = 1, ...) {
  if (length(which) != 1 || which != 1) {
    cli::cli_abort("For now, only basic 'residuals against fitted'
      plot is supported.")
  }
  NextMethod(which = which)
}

#' @noRd
suppress_no_rank_lmranks <- function(expr,
                                     pattern = "no ranked terms") {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(pattern, conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning") # note the capital W
      }
    }
  )
}
