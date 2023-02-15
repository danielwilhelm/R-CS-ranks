check_csranks_args <- function(coverage, cstype, stepdown, R, simul, seed){
  assert_is_single_positive_integer(R, "R")
  assert_is_single_probability(coverage, "coverage")
  assert_is_single_positive_integer(seed, "seed", na_ok=TRUE)
  assert_is_single_boolean(stepdown, "stepdown")
  assert_is_single_boolean(simul, "simul")
  assert_is_one_of(cstype, "cstype", c("two-sided", "upper", "lower"))
}

process_csranks_args <- function(x, Sigma, indices, na.rm){
  # Do some checks and optionally handle NA values.
  assert_is_numeric_vector(x, "x", na_ok=TRUE)
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
  assert_has_no_NAs(Sigma, "Sigma")
  assert_has_no_NAs(x, "x")
  list(x=x,Sigma=Sigma,indices=indices)
}

check_csranks_multinom_args <- function(coverage, cstype, simul, multcorr){
  assert_is_single_probability(coverage, "coverage")
  assert_is_single_boolean(simul, "simul")
  assert_is_one_of(cstype, "cstype", c("two-sided", "upper", "lower"))
  assert_is_one_of(multcorr, "multcorr", c("Bonferroni", "Holm", "none"))
}

process_csranks_multinom_args <- function(x, indices, na.rm){
  # Do some checks and optionally handle NA values.
  assert_is_integer(x, "x", na_ok=TRUE)
  indices <- process_indices_argument(indices, length(x))
  
  if (na.rm) {
    ind <- !is.na(x)
    x <- x[ind]
    indices <- adjust_indices_for_NAs(indices, ind)
  }
  assert_has_no_NAs(x, "x")
  
  if(any(x < 0)){
    wrong_index <- which(x < 0)[1]
    cli::cli_abort(c("{.var x} must be a vector of non negative integers.",
                     "x" = "{.var x[{wrong_index}]} == {x[wrong_index]} is negative."))
  }
  list(x=x, indices=indices)
}

check_tau <- function(tau, p){
  assert_is_single_positive_integer(tau, "tau")
  if(tau > p){
    cli::cli_abort(c("{.var {tau}} must be a single integer between 0 and {.var length(x)}.",
                     "x" = "{.var tau} == {tau} is larger than .{var length(x)} == {length(x)}."))
  }
}

process_indices_argument <- function(indices,p){
  # Do some checks and optionally handle NA values.
  if (any(is.na(indices))) indices <- 1:p
  else{
    assert_is_integer(indices, "indices")
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

check_plotranking_args <- function(ranks, L, U, popnames, title, subtitle,
                                   caption, colorbins, horizontal){
  assert_is_numeric_vector(ranks, "ranks")
  assert_is_integer(L, "L")
  assert_is_positive(L, "L", na_ok = FALSE)
  assert_is_integer(U, "U")
  assert_is_positive(U, "U", na_ok = FALSE)
  assert_is_between(ranks, L, U, "ranks", "L", "U")
  if(!is.null(popnames)){
    assert_is_character(popnames, "popnames")
    assert_length(popnames, "popnames", length(ranks))
  }
  if(!is.null(title)){
    assert_is_character(title, "title")
    assert_is_single(title, "title")
  }
  if(!is.null(subtitle)){
    assert_is_character(subtitle, "subtitle")
    assert_is_single(subtitle, "subtitle")
  }
  if(!is.null(caption)){
    assert_is_character(caption, "caption")
    assert_is_single(caption, "caption")
  }
  assert_is_integer(colorbins, "colorbins")
  assert_is_between(colorbins, 1, length(ranks), "colorbins", "1", "length(ranks)")
  assert_is_single_boolean(horizontal)
}

assert_is_between <- function(middle, lower, upper, middle_name, lower_name, upper_name){
  if(!all(lower <= middle)){
    wrong_index <- which(lower > middle)[1]
    cli::cli_abort(c("{.var {middle_name}} must be between {.var {lower_name}} and {.var {upper_name}}.",
                     "x" = "{.var {lower_name}[{wrong_index}]} == {lower[wrong_index]} is larger than {.var {middle_name}[{wrong_index}]} == middle[wrong_index]}."))
  }
  if(!all(middle <= upper)){
    wrong_index <- which(middle > upper)[1]
    cli::cli_abort(c("{.var {middle_name}} must be between {.var {lower_name}} and {.var {upper_name}}.",
                     "x" = "{.var {upper_name}[{wrong_index}]} == {upper[wrong_index]} is smaller than {.var {middle_name}[{wrong_index}]} == middle[wrong_index]}."))
  }
}

adjust_indices_for_NAs <- function(original_indices, is_not_na){
  # If we are deleting NAs from x, the indices point to different entries (populations) in x than before
  # We have to correct for that
  all_indices <- 1:length(is_not_na)
  filtered_indices <- (all_indices %in% original_indices) & is_not_na
  (all_indices - cumsum(!is_not_na))[filtered_indices]
}

assert_is_single_positive_integer <- function(x, name, na_ok=FALSE){
  assert_is_integer(x, name, na_ok)
  assert_is_single(x, name)
  assert_is_positive(x, name, na_ok)
}

assert_is_positive <- function(x, name, na_ok){
  if(na_ok){
    if(all(is.na(x))) return()
    x <- x[!is.na(x)]
  }
  if(any(x <= 0)){
    msg <- c("{.var {name}} must be positive.")
    negative_index <- which(x <= 0)[1]
    msg["x"] <- "{.var {name}[{negative_index}]} == {x[negative_index]} <= 0."
    cli::cli_abort(msg)
  }
}

assert_is_single_probability <- function(x, name){
  assert_is_numeric_vector(x, name)
  assert_is_single(x, name)
  if(any(x < 0) || any(x > 1)){
    msg <- c("{.var {name}} must be a probability.",
             "i" = "It must be between zero and one.")
    if(any(x < 0)){
      negative_index <- which(x < 0)[1]
      msg["x"] <- "{.var {name}[{negative_index}]} == {x[negative_index]} is below 0."
    } else {
      too_large_index <- which(x > 1)[1]
      msg["x"] <- "{.var {name}[{too_large_index}]} == {x[too_large_index]} is above 1."
    }
    cli::cli_abort(msg)
  }
}

assert_is_integer <- function(x, name, na_ok=FALSE){
  assert_is_numeric_vector(x, name, na_ok)
  if(!na_ok)
    assert_has_no_NAs(x)
  else{
    if(all(is.na(x))) return()
    x <- x[!is.na(x)]
  }
  deviation_from_int <- which.max(abs(x - round(x)))
  max_deviation <- max(abs(x - round(x)))
  
  if(max_deviation > 1e-10)
    cli::cli_abort(c("{.var {name}} must be an integer vector.",
                     "x" = "{.var {name}[{deviation_from_int}]} = {x[deviation_from_int]} is not an integer."))
}

assert_is_numeric_vector <- function(x, name, na_ok=FALSE){
  assert_is_vector(x, name)
  if(!(is.numeric(x) || all(is.na(x)) && na_ok))
    cli::cli_abort(c("{.var {name}} must be a single number.",
                     "x" = "{.var {name}} is of {.cls {typeof(x)}} type."))
}

assert_is_single <- function(x, name){
  assert_length(x, name, 1)
}

assert_length <- function(x, name, length){
  if(length(x) != length)
    cli::cli_abort(c("{.var {name}} must be of length {length}.",
                     "x" = "{.var {name}} is of length {length(x)}."))
}

assert_is_single_boolean <- function(x, name){
  assert_is_vector(x, name)
  assert_is_single(x, name)
  if(!is.logical(x) || is.na(x)) # single NA is a boolean for R. Tri-valued logic
    cli::cli_abort(c("{.var {name}} must be a single boolean.",
                     "x" = "{.var {name}} is of {.cls {typeof(x)}} type."))
}

assert_is_one_of <- function(x, name, choices){
  assert_is_vector(x, name)
  assert_is_single(x, name)
  if(!is.character(x) || is.na(x))
    cli::cli_abort(c("{.var {name}} must be one of: {paste(choices, collapse=', ')}",
                     "x" = "{.var {name}} is of {.cls {typeof(x)}} type."))
  if(!x %in% choices){
    cli::cli_abort(c("{.var {name}} must be one of: {paste(choices, collapse=', ')}",
                     "x" = "{.var {name} is {x}}"))
  }
}

assert_is_character <- function(x, name){
  assert_is_vector(x, name)
  if(!is.character(x))
    cli::cli_abort(c("{.var {name}} must be a character.",
                     "x" = "{.var {name}} is of {.cls {typeof(x)}} type."))
  assert_has_no_NAs(x, name)
}

assert_has_no_NAs <- function(x, name){
  if(any(is.na(x))){
    na_indices <- utils::head(which(is.na(x)))
    cli::cli_abort("NA values found in {.var {name}} at positions {paste(na_indices, collapse=', ')}.")
  }
}

assert_is_vector <- function(x, name){
  if(!is.atomic(x))
    cli::cli_abort(c("{.var {name}} must be a vector.",
                     "x" = "{.var {name}} is of {.cls {class(x)}} class."))
}