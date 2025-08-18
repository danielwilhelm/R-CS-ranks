#' @describeIn ivregranks Summary and Inference Methods for \code{"ivregranks"}
#' Objects
#'
#' @param object An object of class \code{"ivregranks"}.
#'
#' @inheritParams ivreg::summary.ivreg
#' @export
summary.ivregranks <- function(object, vcov. = NULL, df = NULL,
                               diagnostics = NULL, ...) {
  if (!is.null(vcov.)) {
    cli::cli_abort("{.var vcov.} argument is not yet supported. ")
  }
  if (!is.null(df)) {
    cli::cli_abort("{.var df} argument is not yet supported. ")
  }

  object$df.residual <- stats::nobs(object) - length(coef(object))
  outcome <- NextMethod()
  object$df.residual <- NA


  outcome$sigma <- NA
  outcome$r.squared <- NA
  outcome$adj.r.squared <- NA
  outcome$waldtest[] <- NA
  outcome$diagnostics <- NULL

  cov_matrix <- outcome$vcov
  outcome$coefficients[, 2] <- sqrt(diag(cov_matrix))
  outcome$coefficients[, 3] <- outcome$coefficients[, 1] /
    outcome$coefficients[, 2]
  outcome$coefficients[, 4] <- 2 * stats::pnorm(-abs(outcome$coefficients[, 3]))

  colnames(outcome$coefficients)[3:4] <- c("z value", "Pr(>|z|)")

  class(outcome) <- c("summary.ivregranks", class(outcome))

  return(outcome)
}

#' @rdname ivregranks
#' @export
print.summary.ivregranks <- function(x, ...) {
  x$r.squared <- x$adj.r.squared <- 0
  text <- utils::capture.output(NextMethod())
  x$r.squared <- NA
  x$adj.r.squared <- NA

  text <- text[!grepl("^Residual standard error", text)]
  text <- text[!grepl("^Multiple R-Squared", text)]
  text <- text[!grepl("^Wald test", text)]
  text <- as.list(text)
  text$sep <- "\n"
  concatenated_text <- do.call(paste, text)
  cat(concatenated_text)
  return(invisible(x))
}

#' @rdname ivregranks
#'
#' @inheritParams ivreg::confint.ivreg
#' @export
confint.ivregranks <- function(
    object, parm, level = 0.95,
    component = c("stage2", "stage1"), complete = TRUE, vcov. = NULL,
    df = NULL, ...) {
  if (!is.null(vcov.)) {
    cli::cli_abort("{.var vcov.} argument is not yet supported. ")
  }
  if (!is.null(df)) {
    cli::cli_abort("{.var df} argument is not yet supported. ")
  }

  if (missing(parm)) {
    NextMethod(
      object = object, level = level, component = component,
      complete = complete, vcov. = vcov., ...
    )
  } else {
    NextMethod()
  }
}

#' @describeIn ivregranks Calculate Variance-Covariance Matrix for a Fitted
#' \code{ivregranks} object
#'
#' Returns the variance-covariance matrix of the regression coefficients
#' (main parameters) of a fitted \code{ivregranks} object. Its result is
#' theoretically valid and asymptotically consistent, in contrast to naively
#' running \code{vcov(ivreg(...))}.
#'
#' @param complete logical indicating if the full variance-covariance matrix
#' should be returned also in case of an over-determined system where
#' some coefficients are undefined and \code{coef(.)} contains NAs
#' correspondingly. When \code{complete = TRUE}, \code{vcov()} is compatible
#' with \code{coef()} also in this singular case.
#' @importFrom stats vcov
#' @export
vcov.ivregranks <- function(object, component = c("stage2", "stage1"),
                            complete = TRUE, ...) {
  component <- match.arg(component, c("stage2", "stage1"))
  ## default: stage 2
  if (component == "stage2") {
    regressor_dropped_fs <- is.na(coef(object, component = "stage1"))
    Z <- stats::model.matrix(object, component = "instruments")
    if (any(regressor_dropped_fs)) {
      R <- qr.R(qr(Z[, !regressor_dropped_fs]))
    } else {
      R <- qr.R(qr(Z))
    }
    projection_residual_matrix_fs <- calculate_projection_residual_matrix(
      R,
      regressor_dropped_fs,
      length(coef(object, component = "stage1"))
    )
    projection_residuals_fs <- Z %*% projection_residual_matrix_fs

    H1 <- calculate_H1(object, projection_residuals_fs)
    H1_mean <- colMeans(H1)

    object_seqn <- object$object_seqn
    H2 <- calculate_H2(object, projection_residuals_fs, H1_mean)
    H3 <- calculate_H3(object, projection_residual_matrix_fs, H1_mean)
    # print(H3)

    projection_residual_matrix_se <- get_projection_residual_matrix(object_seqn)
    X <- stats::model.matrix(object_seqn)
    projection_residuals_se <- X %*% projection_residual_matrix_se
    projection_variances <- colMeans(projection_residuals_se *
      projection_residuals_fs)
    psi <- t(t(H1 + H2 + H3) / projection_variances)

    sigmahat <- (t(psi) %*% psi) / (nrow(psi)^2)
    colnames(sigmahat) <- names(coef(object_seqn))
    rownames(sigmahat) <- colnames(sigmahat)

    if (!complete) {
      sigmahat <- sigmahat[!is.na(coef(object_seqn)), !is.na(coef(object_seqn))]
    }

    return(sigmahat)
  } else {
    return(vcov(object$object_fs, complete = complete, ...))
  }
}

calculate_H1.ivregranks <- function(object, projection_residuals) {
  NextMethod()
}

calculate_H2.ivregranks <- function(object, projection_residuals, H1_mean = NULL) {
  rank_column_index <- get_ranked_indices(
    object$object_seqn,
    component="regressors"
  )
  model_matrix_seqn <- stats::model.matrix(object$object_seqn)
  l <- get_and_separate_regressors(
    model_matrix_seqn,
    rank_column_index
  )
  RY <- stats::model.response(stats::model.frame(object$object_seqn))

  NextMethod(l = l, RY = RY)
}
calculate_H3.ivregranks <- function(object, projection_residual_matrix,
                                    H1_mean) {
  rank_column_index <- get_ranked_indices(object,
    component= "instruments"
  )
  model_matrix <- stats::model.matrix(object, component = "instruments")
  l <- get_and_separate_regressors(model_matrix, rank_column_index)

  NextMethod(l = l)
}
