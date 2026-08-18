update_formula_intercept <- function(parsed_formula) {
  # We need to replace intercept with the grouping factor
  formula_terms <- parsed_formula$formula_terms
  if (attr(formula_terms, "intercept")) {
    new_terms <- attr(formula_terms, "term.labels")
    interacting_var <- parsed_formula[["does_variable_interact_with_ranked_regressor"]]
    grouping_var_name <- rownames(attr(formula_terms, "factors"))[interacting_var]

    if (!grouping_var_name %in% attr(formula_terms, "term.labels")) {
      new_terms <- c(new_terms, grouping_var_name)
    }

    # Reformulate without intercept, adding grouping variable explicitly
    modified_formula <- stats::reformulate(new_terms,
      response = formula_terms[[2]], intercept = FALSE
    )
    return(modified_formula)
  } else {
    return(parsed_formula$raw_formula)
  }
}

get_rank_terms_indices_after_reordering <- function(formula, rank_terms_names) {
  rearranged_formula_terms <- stats::terms(formula,
    allowDotAsName = TRUE,
    keep.order = FALSE
  )
  which(attr(rearranged_formula_terms, "term.labels") %in% rank_terms_names)
}
