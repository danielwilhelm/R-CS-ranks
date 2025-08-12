#' Title
#'
#' @param object
#' @param vcov.
#' @param df
#' @param diagnostics
#' @param ...
#'
#' @return
#' @export
summary.ivregranks <- function(object, vcov. = NULL, df = NULL,
                               diagnostics = NULL, ...) {

}


#' Title
#'
#' @param object
#' @param component
#' @param complete
#' @param ...
#'
#' @return
#' @export
vcov.ivregranks <- function(object, component = c("stage2", "stage1"),
                            complete = TRUE, ...) {
  regressor_dropped_fs <- is.na(coef(object, component = "stage1"))
  Z <- stats::model.matrix(object, component = "instruments")
  if (any(regressor_dropped_fs)) {
    R <- qr.R(qr(Z[, !regressor_dropped_fs]))
  } else {
    R <- qr.R(Z)
  }
  projection_residual_matrix_fs <- calculate_projection_residual_matrix(
    R,
    regressor_dropped_fs,
    length(coef(object, component = "stage1"))
  )
  projection_residuals_fs <- Z %*% projection_residuals_fs

  H1 <- calculate_H1(object, projection_residuals_fs)
  H1_mean <- colMeans(H1)

  object_se <- object$object_se
  H2 <- calculate_H2(
    object, object_se, projection_residuals_fs,
    H1_mean
  )
  H3 <- calculate_H3(
    object, object_se, object_fs,
    projection_residual_matrix_fs, H1_mean
  )

  projection_residual_matrix_se <- get_projection_residual_matrix(object_se)
  X <- stats::model.matrix(object_se)
  projection_residuals_se <- X %*% projection_residual_matrix_se
  projection_variances <- colMeans(projection_residuals_se *
    projection_residuals_fs)
  psi <- t(t(H1 + H2 + H3) / projection_variances)

  sigmahat <- (t(psi) %*% psi) / (nrow(psi)^2)
  colnames(sigmahat) <- names(coef(object_se))
  rownames(sigmahat) <- colnames(sigmahat)

  if (!complete) {
    sigmahat <- sigmahat[!is.na(coef(object_se)), !is.na(coef(object_se))]
  }

  return(sigmahat)
}
