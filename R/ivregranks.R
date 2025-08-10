#' Title
#'
#' @param formula
#' @param instruments
#' @param data
#' @param subset
#' @param na.action
#' @param weights
#' @param offset
#' @param contrasts
#' @param model
#' @param y
#' @param x
#' @param method
#' @param omega
#' @param ...
#'
#' @return
#' @export
ivregranks <- function(formula, instruments, data, subset, na.action, weights,
                       offset, contrasts = NULL, model = TRUE, y = TRUE,
                       x = FALSE, method = c("OLS", "M", "MM"), omega = 1,
                       ...) {
  rank_env <- create_env_to_interpret_r_mark(omega)
  l <- process_ivregranks_formula(formula, rank_env)
  rank_terms_indices <- l$rank_terms_indices
  ranked_response <- l$ranked_response
  corrected_formula <- l$formula
  original_call <- match.call()
  if (length(rank_terms_indices) == 0 && !ranked_response) {
    cli::cli_warn("{.var ivregranks} called with no ranked terms.
      Using regular ivreg...")
    ivreg_call <- prepare_call(original_call, check_ivreg_args = FALSE)
    out <- eval(ivreg_call, parent.frame())
    return(out)
  }
  ivreg_call <- prepare_call(original_call)
  ivreg_call$formula <- substitute(corrected_formula)

  main_model <- eval(ivreg_call, rank_env)
  if (method == "model.frame") {
    return(main_model)
  }

  main_model$rank_terms_indices <- rank_terms_indices

  main_model$call <- original_call
  main_model$df.residual <- NA
  main_model$omega <- omega
  main_model$ranked_response <- ranked_response
  class(main_model) <- c("ivreg", class(main_model))

  return(main_model)
}

#' Title
#'
#' @param formula
#' @param rank_env
#'
#' @return
process_ivregranks_formula <- function(formula, rank_env = NULL) {
  if (!inherits(formula, "formula")) {
    cli::cli_abort(c("{.var formula} must be a {.class formula} object.",
      "x" = "The passed {.var formula} is of {.cls {class(formula)}} class."
    ))
  }

  if (is.null(rank_env)) {
    rank_env <- environment(formula)
  }

  formula <- Formula::as.Formula(formula)
  if (length(formula)[2] == 1) {
    cli::cli_abort(c("{.var formula} must contain two regressor parts"),
      "x" = "The passed {.var formula} has a single part regressor",
      "i" = "Use lmranks."
    )
  }
  if (length(formula)[1] != 1 || length(formula)[2] > 3) {
    cli::cli_abort(c("{.var formula} must contain a single outcome and at least
      an instrument part.",
      "x" = "The passed {.var formula} has either a multi-part response
      or more than three-part regressors."
    ))
  }
  formula_terms <- stats::terms(formula,
    specials = "r", keep.order = TRUE,
    allowDotAsName = TRUE
  )
  outcome_eq_terms <- stats::terms(formula,
    rhs = 1,
    specials = "r",
    allowDotAsName = TRUE
  )

  l <- process_lmranks_formula(
    Formula::as.Formula(outcome_eq_terms),
    rank_env
  )

  rank_variables_indices <- attr(formula_terms, "specials")[["r"]]
  ranked_instrument_index <- setdiff(
    rank_variables_indices,
    l$rank_terms_indices
  )
  if (length(ranked_instrument_index) > 1) {
    cli::cli_abort(c("In formula there may be at most one ranked instrument."),
      "x" = "There are mulple ranked instruments."
    )
  }

  environment(formula) <- rank_env
  return(list(
    rank_terms_indices = rank_variables_indices,
    ranked_response = l$ranked_response, formula = formula
  ))
}

#' Title
#'
#' @param ivreg_call
#' @param check_ivreg_args
#'
#' @return
prepare_call <- function(ivreg_call, check_ivreg_args = TRUE) {
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

slotsFromS3.ivregranks <- function(object) {
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
}

#' Title
#'
#' @param x
#' @param which
#' @param ...
#'
#' @return
#' @export
plot.ivregranks <- function(x, which = 1, ...) {
  if (length(which) != 1 || which != 1) {
    cli::cli_abort("For now, only basic 'residuals against fitted'
      plot is supported.")
  }
  NextMethod(which = which)
}
