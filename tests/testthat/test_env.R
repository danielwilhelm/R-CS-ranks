test_that("create_env_to_interpret_r_mark has a correct parent env", {
  expected_parent_env <- new.env()
  caller_env <- new.env(parent = expected_parent_env)
  assign("wrapper",
    function(omega) {
      create_env_to_interpret_r_mark(omega = omega)
    },
    envir = expected_parent_env
  )
  created_env <- eval(quote(wrapper(0.4)), envir = expected_parent_env)
  expect_reference(parent.env(created_env), expected_parent_env)
})

test_that("create_env_to_interpret_r_mark's r function behaves correctly in fitting", {
  x_1 <- c(4, 4, 4, 3, 1, 10, 7, 7)
  expected_r_output <- c(0.475, 0.475, 0.475, 0.250, 0.125, 1.000, 0.800, 0.800)
  wrapper <- function(omega) {
    create_env_to_interpret_r_mark(omega = omega)
  }
  created_env <- wrapper(0.4)
  actual_r_output <- eval(quote(r(x_1)), created_env)
  expect_equal(actual_r_output, expected_r_output)

  expected_r_output_om0 <- c(0.375, 0.375, 0.375, 0.250, 0.125, 1.000, 0.750, 0.750)
  created_env <- wrapper(0)
  actual_r_output <- eval(quote(r(x_1)), created_env)
  expect_equal(actual_r_output, expected_r_output_om0)

  expected_r_output_om1 <- c(0.625, 0.625, 0.625, 0.25, 0.125, 1.0, 0.875, 0.875)
  created_env <- wrapper(1)
  actual_r_output <- eval(quote(r(x_1)), created_env)
  expect_equal(actual_r_output, expected_r_output_om1)
})

test_that("create_env_to_interpret_r_mark's r function behaves correctly in prediction", {
  x_1 <- c(4, 4, 4, 3, 1, 10, 7, 7)
  x_pred <- c(0, 1, 2, 3, 4, 5, 7, 8, 10, 11)
  expected_r_output <- ecdf(x_1)(x_pred)
  wrapper <- function(omega) {
    create_env_to_interpret_r_mark(omega = omega)
  }
  created_env <- wrapper(1.0)

  assign(".r_cache", list(x_pred = x_1), envir = created_env)
  assign(".r_predict", TRUE, envir = created_env)

  actual_r_output <- eval(quote(r(x_pred)), created_env)
  expect_equal(actual_r_output, expected_r_output)
})

test_that("create_env_to_interpret_r_mark's r function handles NA in prediction", {
  x_1 <- c(4, 4, 4, 3, 1, 10, 7, 7)
  x_pred <- c(0, 1, NA, 3, 4, 5, 7, NA, 10, 11)
  expected_r_output <- ecdf(x_1)(x_pred)

  wrapper <- function(omega, na.rm) {
    create_env_to_interpret_r_mark(omega = omega)
  }
  created_env <- wrapper(1.0)

  assign(".r_cache", list(x_pred = x_1), envir = created_env)
  assign(".r_predict", TRUE, envir = created_env)

  actual_r_output <- eval(quote(r(x_pred)), created_env)
  expect_equal(actual_r_output, expected_r_output)
})
