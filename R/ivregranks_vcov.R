#' @noRd
update_coefficients_by_dropping_exogenous_vars <- function(
  object, stage_1_projection_residual_coeffiecients
) {
  stage_1_coefficients <- coef(object, component = "stage1")
  regressor_dropped <- is.na(stage_1_coefficients)
  # We skip the instrumental variable Z.
  # Assumption: we only have 1 instrumental variable, which is also ranked.
  instrument_term <- object[["instruments"]][object[["rank_instruments_indices"]]]
  instrument_index <- which((1:length(stage_1_coefficients) == instrument_term)[!regressor_dropped])

  exogeneous_residual_matrix_without_Z_projection <- stage_1_projection_residual_coeffiecients[, -instrument_index]

  stage_1_coefficients_cleaned <- stage_1_coefficients[!regressor_dropped]

  substitute <- exogeneous_residual_matrix_without_Z_projection * stage_1_coefficients_cleaned[-instrument_index][col(exogeneous_residual_matrix_without_Z_projection)] * -1

  stage_1_coefficients_cleaned[row(substitute)] + substitute
}
