parse_lmranks_formula <- function(formula) {
    formula_terms <- stats::terms(formula,
        specials = "r",
        keep.order = TRUE,
        allowDotAsName = TRUE
    )
    variables_terms_table <- attr(formula_terms, "factors")

    rank_variables_indices <- attr(formula_terms, "specials")[["r"]]
    response_variable_index <- attr(formula_terms, "response")
    ranked_regressor_variable_indices <- setdiff(rank_variables_indices, response_variable_index)
    is_response_ranked <- response_variable_index %in% rank_variables_indices
    ranked_regressor_present_in_term <- variables_terms_table[ranked_regressor_variable_indices, ] != 0
    ranked_terms_labels <- colnames(variables_terms_table)[ranked_regressor_present_in_term]
    is_variable_present_in_ranked_term <- variables_terms_table[, ranked_regressor_present_in_term] != 0
    is_variable_present_in_ranked_term[rank_variables_indices] <- FALSE

    list(
        raw_formula = formula,
        formula_terms = formula_terms,
        is_response_ranked = is_response_ranked,
        ranked_regressor_variable_indices = ranked_regressor_variable_indices,
        ranked_regressor_present_in_term = ranked_regressor_present_in_term,
        ranked_terms_labels = ranked_terms_labels,
        does_variable_interact_with_ranked_regressor = is_variable_present_in_ranked_term
    )
}

is_variable_present_in_term <- function(formula_terms, variable_index) {
    variables_table <- attr(formula_terms, "factors")
    variables_table[variable_index, ] != 0
}

has_no_ranked_terms_nor_response <- function(parsed_formula) {
    return(!parsed_formula[["is_response_ranked"]] && length(parsed_formula[["ranked_regressor_variable_indices"]]) == 0)
}
