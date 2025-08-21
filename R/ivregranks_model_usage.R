#' @describeIn ivregranks Predict method for IV Model for Ranks Fits
#' @param object \code{ivregrans} object.
#' @param newdata An optional data frame in which to look for variables with
#' which to predict. If omitted, the fitted values are used.
#' @export
predict.ivregranks <- function(object, newdata, ...) {
  # se.fit = FALSE, scale = NULL,
  # df = Inf, interval = "none", level = 0.95, type = "response",
  # terms = NULL, na.action = na.pass, pred.var = res.var/weights,
  # weights = 1,
  call <- match.call()
  illegal_argument_encountered <- !is.null(call$se.fit) ||
    !is.null(call$scale) || !is.null(call$df) || !is.null(call$interval) ||
    !is.null(call$type)
  # Disable not (yet) supported arguments
  if (illegal_argument_encountered) {
    cli::cli_abort(c("Only {.var object}, {.var newdata} and {.var na.action}
      arguments are currently supported.",
      "i" = "Currently, only basic prediction is supported, without calculation
      of standard error, confidence intervals, or per-term breakdown."
    ))
  }

  rank_env <- environment(object$terms)
  assign(".r_predict", TRUE, rank_env)
  out <- NextMethod()
  assign(".r_predict", FALSE, rank_env)
  out
}
