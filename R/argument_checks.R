process_csranks_args <- function(x, Sigma, indices, na.rm){
  if(!is.vector(x))
    cli::cli_abort(c("{.var x} must be a numeric vector.",
                     "x" = "{.var x} is of {.cls {class(x)}} class."))
  if(!is.numeric(x))
    cli::cli_abort(c("{.var x} must be a numeric vector.",
                     "x" = "{.var x} is of {.cls {typeof(x)}} type."))
  if(!is.matrix(Sigma)){
    msg <- c("{.var Sigma} must be a numeric matrix.",
             "x" = "{.var Sigma} if of {.cls {class(Sigma)}} class.")
    if(is.vector(Sigma) && is.numeric(Sigma))
      msg["i"] <- "Did you provide only a vector of variances? If so, use a diagonal covariance matrix."
    cli::cli_abort(msg)
  }
  if(nrow(Sigma) != length(x))
    cli::cli_abort(c("{.var Sigma} must be a matrix with number of rows and columns equal to length of {.var x}.",
                     "x" = "{.var Sigma} has {nrow(Sigma)} rows, but {.var x} has length of {length(x)}."))
  if(ncol(Sigma) != length(x))
    cli::cli_abort(c("{.var Sigma} must be a matrix with number of rows and columns equal to length of {.var x}.",
                     "x" = "{.var Sigma} has {ncol(Sigma)} columns, but {.var x} has length of {length(x)}."))
  indices = process_indices_argument(indices, length(x))
  # remove NAs
  if (na.rm) {
    ind <- !is.na(x) & apply(Sigma, 1, function(v) all(!is.na(v))) & apply(Sigma, 2, function(v) all(!is.na(v)))
    x <- x[ind]
    Sigma <- Sigma[ind, ind]
    indices <- adjust_indices_for_NAs(indices, ind)
  }
  if(any(is.na(x)))
    cli::cli_abort("NA values found in `x`.")
  if(any(is.na(Sigma)))
    cli::cli_abort("NA values found in `Sigma`.")
  list(x=x,Sigma=Sigma,indices=indices)
}

process_csranks_multinom_arguments <- function(x, indices, na.rm){
  if(!is.vector(x))
    cli::cli_abort(c("{.var x} must be an integer vector.",
                     "x" = "{.var x} is of {.cls {class(x)}} class."))
  if(!is.numeric(x))
    cli::cli_abort(c("{.var x} must be an integer vector.",
                     "x" = "{.var x} is of {.cls {typeof(x)}} type."))
  indices <- process_indices_argument(indices, length(x))
  
  if (na.rm) {
    ind <- !is.na(x)
    x <- x[ind]
    indices <- adjust_indices_for_NAs(indices, ind)
  }
  if(any(is.na(x)))
    cli::cli_abort("NA values found in `x`.")
  
  deviation_from_int <- which.max(abs(x - round(x)))
  max_deviation <- max(abs(x - round(x)))
  
  if(max_deviation > 1e-10)
    cli::cli_abort(c("{.var x} must be an integer vector.",
                     "x" = "{.var x[{deviation_from_int}]} = {x[deviation_from_int]} is not an integer."))
  if(any(x < 0)){
    wrong_index <- which(x < 0)[1]
    cli::cli_abort(c("{.var x} must be a vector of non negative integers.",
                     "x" = "{.var x[{wrong_index}]} == {x[wrong_index]} is negative."))
  }
  list(x=x, indices=indices)
}

process_indices_argument <- function(indices,p){
  if (any(is.na(indices))) indices <- 1:p
  else{
    if(!is.vector(indices))
      cli::cli_abort(c("{.var indices} must be an integer vector.",
                       "x" = "{.var indices} is of {.cls {class(indices)}} class."))
    if(!is.numeric(indices))
      cli::cli_abort(c("{.var indices} must be an integer vector.",
                       "x" = "{.var indices} is of {.cls {typeof(indices)}} type."))
    deviation_from_int <- which.max(abs(indices - round(indices)))
    max_deviation <- max(abs(indices - round(indices)))
    if(max_deviation > 1e-10)
      cli::cli_abort(c("{.var indices} must be an integer vector.",
                       "x" = "{.var indices[{deviation_from_int}]} == {indices[deviation_from_int]} is not an integer."))
    if(any(indices < 0) && any(indices > 0) || any(indices == 0) || any(abs(indices) > p)){
      msg <- c("{.var indices} must be a vector of integer indices of {.var x}.",
               "i" = "They must be either all strictly positive or strictly negative.",
               "i" = "They must be also between {.var -length(x)} and {.var length(x)}.")
      if(any(indices == 0)){
        zero_index <- which(indices == 0)[1]
        msg["x"] <- "{.var indices[{zero_index}]} == 0."
      } else if(any(indices < 0) && any(indices > 0)){
        positive_index <- which(indices > 0)[1]
        negative_index <- which(indices < 0)[1]
        msg["x"] <- "{.var indices[{positive_index}]} == {indices[positive_index]} is positive, but {.var indices[{negative_index}]} == {indices[negative_index]} is negative."
      } else {
        too_large_index <- which(abs(indices) > p)[1]
        msg["x"] <- "{.var indices[{too_large_index}]} == {indices[too_large_index]}, but {.var x} is of length {p}."
      }
      cli::cli_abort(msg)
    }
    # Correct for negative indices
    indices <- (1:p)[indices]
  }
  indices
}

adjust_indices_for_NAs <- function(original_indices, is_not_na){
  all_indices <- 1:length(is_not_na)
  filtered_indices <- (all_indices %in% original_indices) & is_not_na
  (all_indices - cumsum(!is_not_na))[filtered_indices]
}
