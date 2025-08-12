calculate_H1 <- function(object, projection_residuals) {
  UseMethod("calculate_H1")
}
#' Calculate H1 component for covariance estimation
#'
#' Originally defined as h_1(x,y) = (R_Y(Y)-rhoR_X(X)-Wbeta)(R_X(X) - Wgamma)
#'
#' @return n x p matrix
#' @noRd
calculate_H1.default <- function(object, projection_residuals) {
  original_resids <- resid(object)
  projection_residuals * original_resids
}
