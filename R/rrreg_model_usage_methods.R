#' @export
simulate.lmranks <- function(object, nsim = 1, seed = NULL, ...){
  cli::cli_warn("This method might not return correct results.")
  NextMethod()
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
                   "i" = "Currently, only basic prediction is supported, wihtout calculation of standard error, confidence intervals, or per-term breakdown."))
  }
  
  rank_env <- environment(object$terms)
  assign(".r_predict", TRUE, rank_env)
  out <- NextMethod()
  assign(".r_predict", FALSE, rank_env)
  out
}

#' @export
plot.lmranks <- function(x,...){
  cli::cli_warn("Some plots are not correct.")
  NextMethod()
}