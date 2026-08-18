compare_for_tests <- function(i, j, v, omega = 1) {
  omega * (v[i] <= v[j]) + (1 - omega) * (v[i] < v[j])
}


test_that("get_projection_residual_matrix works", {
  X <- as.matrix(mtcars)
  expected <- diag(1, nrow = ncol(X), ncol = ncol(X))
  expected_resids <- matrix(nrow = nrow(X), ncol = ncol(X))
  for (j in 1:ncol(X)) {
    X_minusj <- X[, -j]
    Y <- X[, j]
    m <- lm(Y ~ X_minusj - 1)

    expected[-j, j] <- -coef(m)
    expected_resids[, j] <- resid(m)
  }

  object <- list(qr = qr(X))
  actual <- get_projection_residual_matrix(object)
  expect_equivalent(expected, actual)
  expect_equivalent(expected_resids, X %*% actual)
})

test_that("get_projection_residual_matrix works for singular matrix", {
  data(mtcars)
  X <- as.matrix(mtcars[, -1])
  X[, 2] <- X[, 1]
  expected <- diag(1, nrow = ncol(X), ncol = ncol(X))
  expected[2, -2] <- 0
  expected[, 2] <- NA
  expected_resids <- matrix(nrow = nrow(X), ncol = ncol(X))
  expected_resids[, 2] <- NA
  for (j in 1:ncol(X)) {
    if (j == 2) next
    X_minusj <- X[, -c(2, j)]
    Y <- X[, j]
    m <- lm(Y ~ X_minusj - 1)

    expected[-c(2, j), j] <- -coef(m)
    expected_resids[, j] <- resid(m)
  }

  Y <- mtcars[, 1]
  object <- lm(Y ~ X - 1)
  actual <- get_projection_residual_matrix(object)
  expect_equivalent(expected, actual)
  expect_equivalent(expected_resids, X %*% actual)
})

test_that("get_and_separate_regressors works", {
  data(mtcars)
  model_1 <- lm(mpg ~ disp + cyl + hp, data = mtcars)
  model_1$rank_terms_indices <- 1
  expected_RX <- mtcars$disp
  names(expected_RX) <- rownames(mtcars)
  expected_out_1 <- list(
    RX = expected_RX,
    rank_column_index = 2, # Intercept
    global_RX = expected_RX
  )
  rank_column_index <- get_ranked_indices(model_1, component = "regressors")
  model_matrix <- stats::model.matrix(model_1)
  expect_equal(
    get_and_separate_regressors(model_matrix, rank_column_index),
    expected_out_1
  )

  RX <- mtcars$disp
  Y <- mtcars$mpg
  W <- as.matrix(mtcars[, c("cyl", "hp")])
  model_2 <- lm(Y ~ W + RX)
  model_2$rank_terms_indices <- 2
  expected_out_2 <- expected_out_1
  names(expected_out_2$RX) <- 1:length(RX)
  names(expected_out_2$global_RX) <- 1:length(RX)
  expected_out_2$rank_column_index <- 4
  rank_column_index <- get_ranked_indices(model_2, component = "regressors")
  model_matrix <- stats::model.matrix(model_2)
  expect_equal(
    get_and_separate_regressors(model_matrix, rank_column_index),
    expected_out_2
  )
})
test_that("get_and_separate_regressors works for no ranked regressors", {
  data(mtcars)
  model <- lm(mpg ~ disp + cyl + hp, data = mtcars)
  model$rank_terms_indices <- numeric(0)

  expected_out <- list(
    RX = integer(0),
    rank_column_index = integer(0),
    global_RX = integer(0)
  )

  rank_column_index <- get_ranked_indices(model, component = "regressors")
  model_matrix <- stats::model.matrix(model)
  expect_equal(
    get_and_separate_regressors(model_matrix, rank_column_index),
    expected_out
  )
})

test_that("prepare_mat_om0 works without duplicates when mat has 1 column", {
  data(mtcars)
  v <- nrow(mtcars):1
  mat <- matrix(mtcars$mpg, ncol = 1)

  expected <- sapply(1:length(v), function(i) {
    I <- sapply(1:length(v), function(j) compare_for_tests(i, j, v, omega = 0))
    I %*% mat
  })
  expected <- matrix(expected, ncol = 1)

  actual <- apply(prepare_mat_om0(mat, v), 2, cumsum)
  expect_equivalent(actual, expected)
})

test_that("prepare_mat_om0 works with duplicates when mat has 1 column", {
  data(mtcars)
  v <- sort(mtcars$disp, decreasing = TRUE)
  mat <- matrix(mtcars$mpg, ncol = 1)

  expected <- sapply(1:length(v), function(i) {
    I <- sapply(1:length(v), function(j) compare_for_tests(i, j, v, omega = 0))
    I %*% mat
  })
  expected <- matrix(expected, ncol = 1)

  actual <- apply(prepare_mat_om0(mat, v), 2, cumsum)
  expect_equivalent(actual, expected)
})

test_that("prepare_mat_om1 works with duplicates when mat has 1 column", {
  data(mtcars)
  v <- nrow(mtcars):1
  mat <- matrix(mtcars$mpg, ncol = 1)

  expected <- sapply(1:length(v), function(i) {
    I <- sapply(1:length(v), function(j) compare_for_tests(i, j, v, omega = 1))
    I %*% mat
  })
  expected <- matrix(expected, ncol = 1)

  actual <- apply(prepare_mat_om1(mat, v), 2, cumsum)
  expect_equivalent(actual, expected)
})

test_that("prepare_mat_om0 works with duplicates when mat has many columns", {
  data(mtcars)
  v <- sort(mtcars$disp, decreasing = TRUE)
  mat <- as.matrix(mtcars)

  expected <- sapply(1:length(v), function(i) {
    I <- sapply(1:length(v), function(j) compare_for_tests(i, j, v, omega = 0))
    I %*% mat
  })
  expected <- t(expected)

  actual <- apply(prepare_mat_om0(mat, v), 2, cumsum)
  expect_equivalent(actual, expected)
})

test_that("prepare_mat_om1 works with duplicates when mat has many columns", {
  data(mtcars)
  v <- sort(mtcars$disp, decreasing = TRUE)
  mat <- as.matrix(mtcars)

  expected <- sapply(1:length(v), function(i) {
    I <- sapply(1:length(v), function(j) compare_for_tests(i, j, v, omega = 1))
    I %*% mat
  })
  expected <- t(expected)

  actual <- apply(prepare_mat_om1(mat, v), 2, cumsum)
  expect_equivalent(actual, expected)
})

test_that("ineq_indicator_matmult works for sorted input", {
  data(mtcars)
  v <- sort(mtcars$disp, TRUE)
  mat <- as.matrix(mtcars)

  expected <- sapply(1:length(v), function(i) {
    I <- sapply(1:length(v), function(j) compare_for_tests(i, j, v, omega = 0.4))
    I %*% mat
  })
  expected <- t(expected)

  actual <- ineq_indicator_matmult(v, mat, 0.4)
  expect_equivalent(actual, expected)
})

test_that("ineq_indicator_matmult works for unsorted input", {
  data(mtcars)
  v <- mtcars$disp
  mat <- as.matrix(mtcars)

  expected <- sapply(1:length(v), function(i) {
    I <- sapply(1:length(v), function(j) compare_for_tests(i, j, v, omega = 0.4))
    I %*% mat
  })
  expected <- t(expected)

  actual <- ineq_indicator_matmult(v, mat, 0.4)
  expect_equivalent(actual, expected)
})

test_that("findIntervalIncreasing works", {
  v <- c(10, 7, 7, 4, 4, 4, 3, 1)
  # v_{i_j} >= v_j > v_{i_{j+1}}
  expected <- c(1, 3, 3, 6, 6, 6, 7, 8)
  expect_equal(findIntervalIncreasing(v, left.open = FALSE), expected)
  # v_{i_j} > v_j >= v_{i_{j+1}}
  expected <- c(0, 1, 1, 3, 3, 3, 6, 7)
  expect_equal(findIntervalIncreasing(v, left.open = TRUE), expected)
})
