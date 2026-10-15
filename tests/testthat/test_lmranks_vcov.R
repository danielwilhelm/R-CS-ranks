test_that("summary does not raise errors", {
  mod <- lmranks(r(mpg) ~ r(cyl) + disp, data = mtcars)
  expect_silent(summary(mod))
})

test_that("print.summary.lmranks works", {
  mod <- lmranks(r(mpg) ~ r(cyl) + disp, data = mtcars)
  sumr <- summary(mod)
  testthat::expect_output(print(sumr))
})

test_that("vcov passes shallow checks", {
  model <- lmranks(r(mpg) ~ r(cyl) + disp, data = mtcars)
  V <- vcov(model)

  expect_true(isSymmetric(V))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals$values > 0))
  expect_equal(
    colnames(V),
    c("(Intercept)", "r(cyl)", "disp")
  )
  expect_equal(
    rownames(V),
    c("(Intercept)", "r(cyl)", "disp")
  )
})

test_that("vcov passes shallow checks in no ranked regressors case", {
  model <- lmranks(r(mpg) ~ cyl + disp, data = mtcars)
  V <- vcov(model)

  expect_true(isSymmetric(V))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals$values > 0))
})

test_that("vcov works for singular model matrix", {
  data(mtcars)
  W <- cbind(mtcars$disp, mtcars$disp)
  mod <- lmranks(r(mpg) ~ r(cyl) + W, data = mtcars)
  cov_mat <- vcov(mod)

  expect_true(all(!is.na(cov_mat[1:3, 1:3])))
  expect_true(all(is.na(cov_mat[4, ])))
  expect_true(all(is.na(cov_mat[, 4])))
})

test_that("vcov works for singular model matrix, complete=FALSE", {
  data(mtcars)
  W <- cbind(mtcars$disp, mtcars$disp)
  mod <- lmranks(r(mpg) ~ r(cyl) + W, data = mtcars)
  cov_mat <- vcov(mod, complete = FALSE)

  expect_true(all(!is.na(cov_mat)))
  expect_equal(
    c(nrow(cov_mat), ncol(cov_mat)),
    c(3, 3)
  )
})

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


test_that("get_projection_residual_matrix works for weighted case", {
  X <- as.matrix(mtcars)
  set.seed(2026)
  weights <- runif(nrow(mtcars), max = 2)
  expected <- diag(1, nrow = ncol(X), ncol = ncol(X))
  expected_resids <- matrix(nrow = nrow(X), ncol = ncol(X))
  for (j in 1:ncol(X)) {
    X_minusj <- X[, -j]
    Y <- X[, j]
    m <- lm(Y ~ X_minusj - 1, weights = weights)

    expected[-j, j] <- -coef(m)
    expected_resids[, j] <- resid(m)
  }

  object <- lmranks(r(weights) ~ . - 1, data = mtcars, weights = weights)
  actual <- get_projection_residual_matrix(object)
  expect_equivalent(expected, actual)
  expect_equivalent(expected_resids, X %*% actual)

  object$qr <- NULL
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

test_that("H2 handles no ranked regressor case without error", {
  data(mtcars)
  m <- lmranks(r(mpg) ~ ., data = mtcars)
  projection_residuals <- model.matrix(m) %*% get_projection_residual_matrix(m)
  expect_silent(calculate_H2(m, projection_residuals))
})

test_that("H2 handles not ranked response case without error", {
  data(mtcars)
  m <- lmranks(mpg ~ r(disp) + ., data = mtcars)
  projection_residuals <- model.matrix(m) %*% get_projection_residual_matrix(m)
  expect_silent(calculate_H2(m, projection_residuals))
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with covariates present", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- lmranks(r(Y) ~ r(X) + W - 1)
  sigma2hat.lmranks <- vcov(res)[1, 1] * n
  expect_equal(sigma2hat, sigma2hat.lmranks)
})

######################################################
### High-level checks against by-hand calculations ###
######################################################

test_that("h1 works for ranked regressor with no covariates", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- lmranks(r(Y) ~ r(X))
  proj_residuals <- model.matrix(res) %*% get_projection_residual_matrix(res)

  h1_lmranks <- calculate_H1(res, proj_residuals)
  expect_equivalent(h1_lmranks[, 2], h1)
})

test_that("h2 works for ranked regressor with no covariates", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- lmranks(r(Y) ~ r(X))
  proj_residuals <- model.matrix(res) %*% get_projection_residual_matrix(res)

  h2_lmranks <- calculate_H2(res, proj_residuals)
  expect_equal(h2_lmranks[, 2], h2)
})

test_that("h3 works for ranked regressor with no covariates", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- lmranks(r(Y) ~ r(X))
  proj_residual_matrix <- get_projection_residual_matrix(res)
  H1 <- calculate_H1(res, model.matrix(res) %*% proj_residual_matrix)
  H1_mean <- colMeans(H1)
  h3_lmranks <- calculate_H3(res, proj_residual_matrix, H1_mean)

  expect_equal(h3_lmranks[, 2], h3)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with no covariates", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- lmranks(r(Y) ~ r(X))
  sigma2hat.lmranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat.lmranks)
})

test_that("h1 works for ranked regressor with covariates", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- lmranks(r(Y) ~ r(X) + W)
  proj_residuals <- model.matrix(res) %*% get_projection_residual_matrix(res)

  h1_lmranks <- calculate_H1(res, proj_residuals)
  expect_equivalent(h1_lmranks[, 2], h1)
})

test_that("h2 works for ranked regressor with covariates", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- lmranks(r(Y) ~ r(X) + W)
  proj_residuals <- model.matrix(res) %*% get_projection_residual_matrix(res)

  h2_lmranks <- calculate_H2(res, proj_residuals)
  expect_equal(h2_lmranks[, 2], h2)
})

test_that("h3 works for ranked regressor with covariates", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- lmranks(r(Y) ~ r(X) + W)
  proj_residual_matrix <- get_projection_residual_matrix(res)
  H1 <- calculate_H1(res, model.matrix(res) %*% proj_residual_matrix)
  H1_mean <- colMeans(H1)
  h3_lmranks <- calculate_H3(res, proj_residual_matrix, H1_mean)

  expect_equal(h3_lmranks[, 2], h3)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with covariates", {
  load(test_path("testdata", "lmranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- lmranks(r(Y) ~ r(X) + W)
  sigma2hat.lmranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat.lmranks)
})

test_that("vcov produces asymptotic variance estimate of rank-rank slope equal to that of Hoeffding (1948)", {
  load(test_path("testdata", "lmranks_cov_sigmahat_Hoefding.rda"))
  res <- lmranks(r(Y) ~ r(X))
  sigma2hat.lmranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat.lmranks, tolerance = 1e-3)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with smaller datasets", {
  load(test_path("testdata", "lmranks_cov_sigmahat_n_10.rda"))
  res <- lmranks(r(Y) ~ r(X) + W)
  sigma2hat.lmranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat.lmranks)

  load(test_path("testdata", "lmranks_cov_sigmahat_n_50.rda"))
  res <- lmranks(r(Y) ~ r(X) + W)
  sigma2hat.lmranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat.lmranks)

  load(test_path("testdata", "lmranks_cov_sigmahat_n_100.rda"))
  res <- lmranks(r(Y) ~ r(X) + W)
  sigma2hat.lmranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat.lmranks)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with increasing=FALSE", {
  load(test_path("testdata", "lmranks_cov_sigmahat_increasing_FALSE.rda"))
  res <- lmranks(r(Y, increasing = FALSE) ~ r(X, increasing = FALSE) + W, omega = 1)
  sigma2hat.lmranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat.lmranks)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope with weights", {
  load(test_path("testdata", "lmranks_cov_sigmahat_weighted.rda"))
  res <- lmranks(r(Y) ~ r(X) + W, weights = weights)
  sigma2hat.lmranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat.lmranks)
})

test_that("h1 works for ranked regressor with weights", {
  load(test_path("testdata", "lmranks_cov_sigmahat_weighted.rda"))
  res <- lmranks(r(Y) ~ r(X) + W, weights = weights)
  proj_residuals <- model.matrix(res) %*% get_projection_residual_matrix(res)

  h1_lmranks <- calculate_H1(res, proj_residuals)
  expect_equivalent(h1_lmranks[, 2], h1)
})

test_that("h2 works for ranked regressor with weights", {
  load(test_path("testdata", "lmranks_cov_sigmahat_weighted.rda"))
  res <- lmranks(r(Y) ~ r(X) + W, weights = weights)
  proj_residuals <- model.matrix(res) %*% get_projection_residual_matrix(res)

  h2_lmranks <- calculate_H2(res, proj_residuals)
  expect_equal(h2_lmranks[, 2], h2)
})

test_that("h3 works for ranked regressor with weights", {
  load(test_path("testdata", "lmranks_cov_sigmahat_weighted.rda"))
  res <- lmranks(r(Y) ~ r(X) + W, weights = weights)
  proj_residual_matrix <- get_projection_residual_matrix(res)
  H1 <- calculate_H1(res, model.matrix(res) %*% proj_residual_matrix)
  H1_mean <- colMeans(H1)
  h3_lmranks <- calculate_H3(res, proj_residual_matrix, H1_mean)

  expect_equal(h3_lmranks[, 2], h3)
})
