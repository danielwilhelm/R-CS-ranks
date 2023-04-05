#' @export
rstandard.lmranks <- function(model, ...){
  cli::cli_warn("This method might not return correct results.")
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