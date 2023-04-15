### lmranks ###

test_that("lmranks and by-hand calculations provide same results",{
  df <- mtcars[1:10,]
  model <- lmranks(r(mpg) ~ r(hp) + disp + cyl + 1, data=df)
  
  expected_response <- c(0.5, 0.5, 0.2, 0.4, 0.8, 0.9, 1.0, 0.1, 0.2, 0.7)
  expected_model.matrix <- matrix(c(
    1, 0.4, 160.0, 6,
    1, 0.4, 160.0, 6,
    1, 0.9, 108.0, 4,
    1, 0.4, 258.0, 6,
    1, 0.2, 360.0, 8,
    1, 0.7, 225.0, 6,
    1, 0.1, 360.0, 8,
    1, 1.0, 146.7, 4,
    1, 0.8, 140.8, 4,
    1, 0.3, 167.6, 6
  ), byrow = TRUE, nrow = 10)
  expected_coef <- solve(t(expected_model.matrix) %*% expected_model.matrix) %*% 
    t(expected_model.matrix) %*% expected_response
  
  expect_equivalent(model.response(model.frame(model)),
                    expected_response)
  expect_equivalent(model.matrix(model),
                    expected_model.matrix)
  expect_equivalent(coef(model),
                    expected_coef)
})

test_that("lmranks and lm provide coherent results", {
  Y <- c(3,1,2,4,5)
  y_frank <- c(0.6, 1.0, 0.8, 0.4, 0.2)
  X <- 1:5
  omega <- 0.5
  x_frank <- c(1.0, 0.8, 0.6, 0.4, 0.2)
  W <- c(1,3,2,5,4)
  
  rank_m <- lmranks(r(Y) ~ r(X) + W)
  raw_rank_m <- unclass(rank_m)
  raw_rank_m$call <- as.character(raw_rank_m$call)
  raw_rank_m$terms <- NULL
  attr(raw_rank_m$model, "terms") <- NULL
  raw_rank_m$omega <- NULL
  raw_rank_m$rank_terms_indices <- NULL
  
  m <- lm(y_frank ~ x_frank + W)
  expected_m <- unclass(m)
  expected_m$df.residual <- NA
  expected_m$call <- as.character(str2lang("lmranks(r(Y) ~ r(X) + W)"))
  expected_m$terms <- NULL
  expected_m$ranked_response <- TRUE
  attr(expected_m$model, "terms") <- NULL 
  names(expected_m$coefficients)[2] <- "r(X)"
  names(expected_m$effects)[2] <- "r(X)"
  dimnames(expected_m$qr$qr)[[2]][2] <- "r(X)"
  colnames(expected_m$model)[1:2] <- c("r(Y)", "r(X)")
  
  expect_equal(raw_rank_m, expected_m)
})

test_that("process_lmranks_formula catches illegal formulas", {
  expect_error(process_lmranks_formula("y ~ x + w"))
  expect_error(process_lmranks_formula(r(y) ~ r(x) + r(w)))
  expect_error(process_lmranks_formula(r(y) ~ r(x) * w))
  expect_error(process_lmranks_formula(r(y) ~ r(x) + r(x):w + w))
  
  expect_silent(process_lmranks_formula(r(y) ~ r(x) + w))
})

test_that("process_lmranks_formula returns correct indices", {
  expect_equal(process_lmranks_formula(r(y) ~ r(x) + w),
               list(rank_terms_indices = 1,
                    ranked_response = TRUE))
  expect_equal(process_lmranks_formula(r(y) ~ w * z + r(x)),
               list(rank_terms_indices = 3,
                    ranked_response = TRUE))
  expect_equal(process_lmranks_formula(r(y) ~ w + z + w:z + r(x)),
               list(rank_terms_indices = 3,
                    ranked_response = TRUE))
  expect_equal(process_lmranks_formula(r(y) ~ w * z + r(x) - z),
               list(rank_terms_indices = 2,
                    ranked_response = TRUE))
  expect_equal(process_lmranks_formula(r(y) ~ w * z),
               list(rank_terms_indices = integer(0),
                    ranked_response = TRUE))
  expect_equal(process_lmranks_formula(y ~ w * z + r(x)),
               list(rank_terms_indices = 3,
                    ranked_response = FALSE))
  expect_equal(process_lmranks_formula(y ~ w * z),
               list(rank_terms_indices = integer(0),
                    ranked_response = FALSE))
})

test_that("process_lmranks_formula returns correct index for simplest fits", {
  expect_equal(process_lmranks_formula(r(y) ~ r(x) - 1),
               list(rank_terms_indices = 1,
                    ranked_response = TRUE))
  expect_equal(process_lmranks_formula(r(y) ~ r(x)),
               list(rank_terms_indices = 1,
                    ranked_response = TRUE))
})

test_that("prepare_lm_call works", {
  input_call <- str2lang("lmranks(r(y) ~ r(x) + W, data=data)")
  expected_call <- str2lang("stats::lm(r(y) ~ r(x) + W, data=data)")
  expect_equal(prepare_lm_call(input_call),
               expected_call)
  
  input_call <- str2lang("lmranks(r(y) ~ r(x) + W, data=data, omega=omega)")
  expected_call <- str2lang("stats::lm(r(y) ~ r(x) + W, data=data)")
  expect_equal(prepare_lm_call(input_call),
               expected_call)
  
  input_call <- str2lang("lmranks(r(y) ~ r(x) + W, data=data, na.rm=na.rm)")
  expected_call <- str2lang("stats::lm(r(y) ~ r(x) + W, data=data)")
  expect_equal(prepare_lm_call(input_call),
               expected_call)
})

test_that("create_env_to_interpret_r_mark has a correct parent env", {
  expected_parent_env <- new.env()
  caller_env <- new.env(parent = expected_parent_env)
  assign("wrapper", 
         function(omega, na.rm){create_env_to_interpret_r_mark(omega=omega,
                                                               na.rm=na.rm)},
         envir = expected_parent_env)
  created_env <- eval(quote(wrapper(0.4, FALSE)), envir = expected_parent_env)
  expect_reference(parent.env(created_env), expected_parent_env)
})

test_that("create_env_to_interpret_r_mark has correct contents", {
  expected_env_contents <- list(.r_predict = FALSE,
                                .r_cache = list(),
                                r = function(x, increasing=TRUE){})
  
  created_env <- create_env_to_interpret_r_mark(omega=0.4,na.rm=FALSE)
  actual_env_contents <- as.list(created_env, all.names = TRUE)
  expect_equal(names(actual_env_contents), 
               names(expected_env_contents))
  expect_equal(actual_env_contents[c(".r_predict", ".r_cache")],
               expected_env_contents[c(".r_predict", ".r_cache")])
  expect_true(is.function(actual_env_contents[['r']]))
})

test_that("create_env_to_interpret_r_mark's r function behaves correctly in fitting", {
  x_1 <- c(4,4,4,3,1,10,7,7)
  expected_r_output <- c(0.6, 0.6, 0.6, 0.875, 1, 0.125, 0.3, 0.3)
  wrapper <- function(omega, na.rm){create_env_to_interpret_r_mark(omega=omega,
                                                               na.rm=na.rm)}
  created_env <- wrapper(0.4, FALSE)
  actual_r_output <- eval(quote(r(x_1)), created_env)
  expect_equal(actual_r_output, expected_r_output)
})

test_that("create_env_to_interpret_r_mark's r function behaves correctly in prediction", {
  x_1 <- c(4,4,4,3,1,10,7,7)
  x_pred <- c(0,1,2,3,4,5,7,8,10,11)
  expected_r_output <- c(8.5, 8, 7.5, 7, 5, 3.5, 2.5, 1.5, 1, 0.5) / 8
  
  wrapper <- function(omega, na.rm){create_env_to_interpret_r_mark(omega=omega,
                                                                   na.rm=na.rm)}
  created_env <- wrapper(0.5, FALSE)
  
  assign(".r_cache", list(x_pred = x_1), envir = created_env)
  assign(".r_predict", TRUE, envir = created_env)
  
  actual_r_output <- eval(quote(r(x_pred)), created_env)
  expect_equal(actual_r_output, expected_r_output)
})

test_that("create_env_to_interpret_r_mark's r function handles NA", {
  x_1 <- c(4,4,NA,4,3,NA,1,10,7,7)
  expected_r_output <- c(0.6, 0.6, NA, 0.6, 0.875, NA, 1, 0.125, 0.3, 0.3)
  wrapper <- function(omega, na.rm){create_env_to_interpret_r_mark(omega=omega,
                                                                   na.rm=na.rm)}
  created_env <- wrapper(0.4, TRUE)
  actual_r_output <- eval(quote(r(x_1)), created_env)
  expect_equal(actual_r_output, expected_r_output)
})

test_that("omega argument works", {
  Y <- c(4,4,4,3,1,10,7,7)
  y_frank <- c(0.6, 0.6, 0.6, 0.875, 1, 0.125, 0.3, 0.3)
  X <- 1:8
  omega <- 0.5
  x_frank <- c(1.0, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125)
  W <- matrix(c(1,4,3,2,5,8,7,6), ncol = 1)
  
  rank_m <- lmranks(r(Y) ~ r(X) + W, omega = 0.4, y = TRUE)
  names(rank_m$y) <- NULL
  expect_equal(rank_m$y, y_frank)
})

