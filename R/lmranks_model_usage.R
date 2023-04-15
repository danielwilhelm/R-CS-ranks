#' @export
simulate.lmranks <- function(object, nsim = 1, seed = NULL, ...){
  cli::cli_abort(c("This method does not return correct results.",
                 i = "The notion of prediction standard error is not theoretically developped yet."))
}

#' @export
predict.lmranks <- function(object, newdata, se.fit = FALSE, scale = NULL, 
                            df = Inf, interval = "none", level = 0.95, type = "response",
                            terms = NULL, na.action = na.pass, pred.var = res.var/weights,
                            weights = 1, ...){
  # Disable not (yet) supported arguments 
  if(se.fit || !is.null(scale) || !is.infinite(df) || interval != "none" || 
     type != "response"){
    cli::cli_error(c("Only {.var object}, {.var newdata} and {.var na.action} arguments are currently suppoerted.",
                   "i" = "Currently, only basic prediction is supported, without calculation of standard error, confidence intervals, or per-term breakdown."))
  }
  
  rank_env <- environment(object$terms)
  assign(".r_predict", TRUE, rank_env)
  out <- NextMethod()
  assign(".r_predict", FALSE, rank_env)
  out
}

#' @export
plot.lmranks <- function(x,which = 1,...){
  if(length(which) != 1 || which != 1)
    cli::cli_abort('For now, only basic "residuals against fitted" plot is supported.')
  NextMethod(which = which)
}