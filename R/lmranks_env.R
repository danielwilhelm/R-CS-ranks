#' Create environment to interpret lmranks formula
#'
#' A formula in lmranks has the ranked variables (regressors and response) marked with
#' `r()`. A way to interpret this mark is needed. In R for this purpose we have
#' *environments* and *non-standard evaluation*.
#'
#' For all intents and purposes it is enough to know, which variables have been marked
#' (that's known from `process_lmranks_formula` output)
#' and values of ranks - that's done by treating `r()` as regular function and
#' evaluating it in an environment containing its correct definition.
#'
#' One difficulty is that the `r()` must return different values depending on whether
#' the linear model is being fitted or used for prediction. In *both* cases we want
#' to use fitting (training) data. For prediction, it is stored in cache.
#'
#' Every time the `lmranks` is called, a new environment of this kind is created
#' and "carried along" with an `lmranks` object (accessible with `environment(model$terms)`).
#'
#' The advantage of this solution is better reuse of `lm.fit` and `predict.lm`.
#'
#' @return an R environment, used later to call the `lm` function in it.
#' In total, it has 3 elements:
#' - `r`, the correct definition of `r` function, a call to frank_against
#' with preprocessed, correct arguments.
#' - `.r_cache`, a list with values saved for possible prediction later.
#' - `.r_predict`, a logical indicating whether we are in fitting (FALSE) or prediction mode.
#' In fitting mode, input data will be used for ranking and will be saved in cache;
#' in prediction mode, the data in cache will be used.
#'
#' Its parent environment is the parent frame of the caller of this function.
#' In the primary use case it is the environment of caller of `lmranks`.
#' In this way we ensure correct evaluation of other formula terms and variables.
#'
#' @seealso
#' [H. Wickham, Advanced R, Environments chapter](https://adv-r.hadley.nz/environments.html)
#' [environment()]
#' [csranks::frank_against()]
#' @noRd
create_env_to_interpret_r_mark <- function(omega, weights = NULL) {
  rank_env <- new.env(parent = parent.frame(2))
  r <- function(x, increasing = TRUE) x
  body(r) <- bquote({
    predict <- get(".r_predict", envir = environment(r), inherits = FALSE)
    cache <- get(".r_cache", envir = environment(r), inherits = FALSE)
    weights <- get(".weights", envir = environment(r), inherits = FALSE)
    was_na <- is.na(x)
    var_name <- paste0(as.character(substitute(x)), collapse = "")
    if (!predict) {
      cache[[var_name]] <- x
      assign(".r_cache", cache, envir = environment(r))
    } else if (is.null(cache[[var_name]])) {
      cli::cli_warn("New variable at predict time. Ranks will be calculated from scratch.")
    }
    v <- cache[[var_name]]
    out <- csranks::frank_against(x, v, increasing = increasing, omega = .(omega), na.rm = FALSE, weights = weights)
    out
  })
  environment(r) <- rank_env
  assign("r", r, envir = rank_env)
  assign(".r_cache", list(), envir = rank_env)
  assign(".r_predict", FALSE, envir = rank_env)
  assign(".weights", weights, envir = rank_env)
  return(rank_env)
}
