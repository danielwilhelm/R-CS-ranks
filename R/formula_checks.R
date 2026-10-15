assert_is_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    cli::cli_abort(c("{.var formula} must be a {.class formula} object.",
      "x" = "The passed {.var formula} is of {.cls {class(formula)}} class."
    ))
  }
}

assert_has_ranked_response <- function(parsed_formula) {
  if (!parsed_formula[["is_response_ranked"]]) {
    cli::cli_abort("In formula the response must be ranked.")
  }
}

assert_has_at_most_one_ranked_regressor <- function(parsed_formula) {
  if (length(parsed_formula[["ranked_regressor_variable_indices"]]) > 1) {
    cli::cli_abort(c("In formula there may be at most one term with ranked regressor.",
      "x" = "There are multiple terms with ranked regressors."
    ))
  }
}

assert_has_exactly_one_ranked_regressor <- function(parsed_formula) {
  if (length(parsed_formula[["ranked_regressor_variable_indices"]]) != 1) {
    cli::cli_abort(c("In formula there must be exactly one term with ranked regressor.",
      "x" = "There are multiple ranked regressors."
    ))
  }
}

assert_has_exactly_one_ranked_regressor_equal_to <- function(parsed_formula, expected) {
  variables <- as.character(attr(parsed_formula[["formula_terms"]], "variables"))[-1]
  ranked_variables <- variables[parsed_formula[["ranked_regressor_variable_indices"]]]

  if (length(ranked_variables) != length(expected) || any(ranked_variables != expected)) {
    cli::cli_abort(c("The following variable is expected to be ranked: {.var {expected}}",
      "x" = "The following variables are acutally ranked: {.var {ranked_variables}}"
    ))
  }
}

assert_has_at_most_one_term_with_ranked_regressor <- function(parsed_formula) {
  rank_regressor_present_in_only_1_term <- sum(parsed_formula[["ranked_regressor_present_in_term"]]) == 1
  if (!rank_regressor_present_in_only_1_term) {
    cli::cli_abort(c("In formula, the ranked regressor may occur only once, as a standalone term or interacting with a global, grouping variable.",
      "x" = "There are multiple terms involving the ranked regressor."
    ))
  }
}

assert_ranked_regressor_does_not_interact_with_other_variables <- function(parsed_formula) {
  # Assume single ranked regressor
  order_of_term_with_ranked_regressor <- sum(parsed_formula[["does_variable_interact_with_ranked_regressor"]])
  if (order_of_term_with_ranked_regressor >= 1) {
    cli::cli_abort(c("Ranked regressors cannot be part of any interactions."))
  }
}

assert_ranked_regressor_interacts_only_with_global_grouping_variable <- function(parsed_formula) {
  one_interacting_var_present <- sum(parsed_formula[["does_variable_interact_with_ranked_regressor"]]) == 1
  if (!one_interacting_var_present) {
    cli::cli_abort(c("In formula, the ranked regressor may occur only once, as a standalone term or interacting with a global, grouping variable.",
      "x" = "The ranked regressor interacts with multiple variables in a single term."
    ))
  }

  interacting_var <- parsed_formula[["does_variable_interact_with_ranked_regressor"]]
  interacting_var_present_in_every_term <- all(is_variable_present_in_term(parsed_formula$formula_terms, which(interacting_var)))

  if (!interacting_var_present_in_every_term) {
    cli::cli_abort(c("In formula, the ranked regressor may occur only once, as a standalone term or interacting with a global, grouping variable.",
      "x" = "The grouping variable does not interact with every other term in the formula."
    ))
  }
}
