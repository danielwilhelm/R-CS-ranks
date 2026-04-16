#' @title Interaction Validation Functions for Rank Regressions
#' @description Functions to validate and process interaction patterns in rank regression formulas
#' @noRd

#' @title Prohibit All Interactions
#' @description Validator that prohibits any interactions involving ranked regressors
#' @param formula_terms Terms object from stats::terms()
#' @param ranked_regressor_index Index of ranked variable (could be empty)
#' @param formula Original formula
#' @return Original formula if validation passes, throws error otherwise
#' @noRd
prohibit_interactions <- function(formula, ranked_regressor_index) {
    # If no ranked regressors, validation passes
    if (length(ranked_regressor_index) == 0) {
        return(formula)
    }

    formula_terms <- stats::terms(formula,
        specials = "r",
        keep.order = TRUE,
        allowDotAsName = TRUE
    )

    variables_table <- attr(formula_terms, "factors")
    ranked_regressor_present_in_term <- variables_table[ranked_regressor_index, ] != 0
    order_of_terms_with_ranked_regressor <- attr(formula_terms, "order")[ranked_regressor_present_in_term]

    # Assume single ranked regressor
    if (length(order_of_terms_with_ranked_regressor) > 1 || order_of_terms_with_ranked_regressor > 1) {
        violating_terms <- attr(formula_terms, "term.labels")[ranked_regressor_present_in_term & (attr(formula_terms, "order") > 1)]
        cli::cli_abort(c("Ranked regressors cannot be part of any interactions.",
            "x" = "The following interactions contain the ranked regresor: {.var violatin_terms}"
        ))
    }

    return(formula)
}



#' @title Allow Only Grouping Interactions
#' @description Validator that allows ranked regressors only in interactions with grouping variables
#' @param formula Original formula
#' @param ranked_regressor_index Index of ranked variable (could be empty)
#' @return Modified formula if grouping interaction detected, original formula otherwise
#' @noRd
allow_only_grouping_interaction <- function(formula, ranked_regressor_index) {
    if (length(ranked_regressor_index) == 0) {
        return(formula)
    }

    formula_terms <- stats::terms(formula,
        specials = "r",
        keep.order = TRUE,
        allowDotAsName = TRUE
    )
    variables_table <- attr(formula_terms, "factors")

    # Check interaction patterns
    rank_regressor_present_in_term <- variables_table[ranked_regressor_index, ] != 0
    rank_regressor_present_in_only_1_term <- sum(rank_regressor_present_in_term) == 1

    if (!rank_regressor_present_in_only_1_term) {
        cli::cli_abort(c("In formula, the ranked regressor may occur only once, as a standalone term or interacting with a global, grouping variable.",
            "x" = "There are multiple terms involving the ranked regressor."
        ))
    }

    is_variable_present_in_ranked_term <- variables_table[, rank_regressor_present_in_term] != 0
    is_ranked_regressor_alone_in_term <- sum(is_variable_present_in_ranked_term) == 1

    if (!is_ranked_regressor_alone_in_term) {
        one_interacting_var_present <- sum(is_variable_present_in_ranked_term) == 2
        if (!one_interacting_var_present) {
            cli::cli_abort(c("In formula, the ranked regressor may occur only once, as a standalone term or interacting with a global, grouping variable.",
                "x" = "The ranked regressor interacts with multiple variables in a single term."
            ))
        }

        interacting_var <- is_variable_present_in_ranked_term
        interacting_var[ranked_regressor_index] <- FALSE
        interacring_var_present_in_every_term <- all(variables_table[interacting_var, ] != 0)

        if (!interacring_var_present_in_every_term) {
            cli::cli_abort(c("In formula, the ranked regressor may occur only once, as a standalone term or interacting with a global, grouping variable.",
                "x" = "The grouping variable does not interact with every other term in the formula."
            ))
        }

        # We need to replace intercept with the grouping factor
        if (attr(formula_terms, "intercept")) {
            new_terms <- attr(formula_terms, "term.labels")
            grouping_var_name <- rownames(variables_table)[interacting_var]

            if (!grouping_var_name %in% attr(formula_terms, "term.labels")) {
                new_terms <- c(new_terms, grouping_var_name)
            }

            # Reformulate without intercept, adding grouping variable explicitly
            modified_formula <- stats::reformulate(new_terms,
                response = formula_terms[[2]], intercept = FALSE
            )
            return(modified_formula)
        }
    }

    return(formula)
}

#' @title Formula Processor Factory
#' @description Creates a formula processor with specified interaction validation
#' @param interaction_validator Function to validate and process interactions
#' @return A formula processing function
#' @noRd
make_formula_processor <- function(interaction_validator) {
    function(formula, rank_env = NULL) {
        if (!inherits(formula, "formula")) {
            cli::cli_abort(c("{.var formula} must be a {.class formula} object.",
                "x" = "The passed {.var formula} is of {.cls {class(formula)}} class."
            ))
        }

        if (is.null(rank_env)) {
            rank_env <- environment(formula)
        }

        # Rest of the processing logic from original process_lmranks_formula
        # We need to re-parse the terms in case the formula was modified
        formula_terms <- stats::terms(formula,
            specials = "r",
            keep.order = TRUE,
            allowDotAsName = TRUE
        )

        rank_variables_indices <- attr(formula_terms, "specials")[["r"]]
        response_variable_index <- attr(formula_terms, "response")
        ranked_regressor_variable_index <- setdiff(rank_variables_indices, response_variable_index)

        if (length(ranked_regressor_variable_index) > 1) {
            cli::cli_abort(c("In formula there may be at most one term with ranked regressor.",
                "x" = "There are multiple terms with ranked regressors."
            ))
        }

        is_response_ranked <- response_variable_index %in% rank_variables_indices

        if (length(ranked_regressor_variable_index) == 0) {
            environment(formula) <- rank_env
            return(list(
                rank_terms_indices = integer(0),
                ranked_response = is_response_ranked,
                formula = formula
            ))
        }

        # Process interactions using the injected validator
        processed_formula <- interaction_validator(formula, ranked_regressor_variable_index)


        variables_terms_table <- attr(formula_terms, "factors")
        rank_regressor_present_in_term <- variables_terms_table[ranked_regressor_variable_index, ] != 0

        # Find rank terms for the final output
        rank_terms_names <- colnames(variables_terms_table)[rank_regressor_present_in_term]
        rearranged_formula_terms <- stats::terms(processed_formula,
            allowDotAsName = TRUE,
            keep.order = FALSE
        )
        rank_terms_indices <- which(attr(rearranged_formula_terms, "term.labels") %in% rank_terms_names)

        environment(processed_formula) <- rank_env
        return(list(
            rank_terms_indices = rank_terms_indices,
            ranked_response = is_response_ranked,
            formula = processed_formula
        ))
    }
}
