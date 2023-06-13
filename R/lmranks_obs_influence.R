#' @importFrom stats cooks.distance dfbeta dfbetas influence rstandard rstudent
#' @export
rstandard.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                   "i" = "The notion of prediction standard error is not theoretically developped yet."))
  NextMethod()
}

#' @export
cooks.distance.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then."))
}

#' @export
dfbeta.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then."))
}

#' @export
dfbetas.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then."))
}

#' @export
influence.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then."))
}

#' @export
rstudent.lmranks <- function(model, ...){
  cli::cli_abort(c("This method returns incorrect results.",
                 "i" = "Deleting and observation changes values of ranks in (in)dependent values. Standard *smart* method breaks down then."))
}