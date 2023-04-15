#' @export
simulate.lmranks <- function(object, nsim = 1, seed = NULL, ...){
  cli::cli_abort(c("This method does not return correct results.",
                 i = "The notion of prediction standard error is not theoretically developped yet."))
}

#' @describeIn lmranks Predict method for Linear Model for Ranks Fits
#' @param object \code{lmranks} object.
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @export
predict.lmranks <- function(object, newdata, ...){
  # se.fit = FALSE, scale = NULL, 
  # df = Inf, interval = "none", level = 0.95, type = "response",
  # terms = NULL, na.action = na.pass, pred.var = res.var/weights,
  # weights = 1,
  call <- match.call()
  illegal_argument_encountered <- !is.null(call$se.fit) || !is.null(call$scale) ||
    !is.null(call$df) || !is.null(call$interval) || !is.null(call$type)
  # Disable not (yet) supported arguments 
  if(illegal_argument_encountered){
    cli::cli_error(c("Only {.var object}, {.var newdata} and {.var na.action} arguments are currently suppoerted.",
                   "i" = "Currently, only basic prediction is supported, without calculation of standard error, confidence intervals, or per-term breakdown."))
  }
  
  rank_env <- environment(object$terms)
  assign(".r_predict", TRUE, rank_env)
  out <- NextMethod()
  assign(".r_predict", FALSE, rank_env)
  out
}