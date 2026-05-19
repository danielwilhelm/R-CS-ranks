test_that("summary does not raise errors", {
  model <- ivregranks(r(mpg) ~ r(hp) + cyl | r(disp) + cyl, data = mtcars)
  expect_silent(summary(model))
})

test_that("print.summary.ivregranks works", {
  model <- ivregranks(r(mpg) ~ r(hp) + cyl | r(disp) + cyl, data = mtcars)
  sumr <- summary(model)
  expect_output(print(sumr))
})

test_that("vcov passes shallow checks", {
  model <- ivregranks(r(mpg) ~ r(hp) + cyl | r(disp) + cyl, data = mtcars)
  V <- vcov(model) # stage2
  expect_true(isSymmetric((V)))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals$values > 0))
  expect_equal(colnames(V), c("(Intercept)", "r(hp)", "cyl"))
  expect_equal(rownames(V), c("(Intercept)", "r(hp)", "cyl"))

  V <- vcov(model, component = "stage1") # stage1
  expect_true(isSymmetric((V)))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals$values > 0))
  expect_equal(colnames(V), c("(Intercept)", "r(disp)", "cyl"))
  expect_equal(rownames(V), c("(Intercept)", "r(disp)", "cyl"))
})

test_that("vcov passes shallow checks in no ranked  case", {
  testthat::skip("That's actually a separate use case from out theory. Need to consult.")
  model <- ivregranks(r(mpg) ~ hp + cyl | disp + cyl, data = mtcars)
  V <- vcov(model)

  expect_true(isSymmetric(V))
  vals <- eigen(V, only.values = TRUE)
  expect_true(all(vals$values > 0))
})

test_that("vcov works for singular model matrix", {
  # that XtX is singular
  w <- cbind(mtcars$qsec, mtcars$qsec)
  model <- expect_warning(ivregranks(r(mpg) ~ r(hp) + w | r(disp) + w, data = mtcars), "collinear")
  cov2 <- vcov(model, component = "stage2")
  cov1 <- vcov(model, component = "stage1")

  expect_true(all(!is.na(cov2[1:3, 1:3])))
  expect_true(all(is.na(cov2[4, ])))
  expect_true(all(is.na(cov2[, 4])))

  expect_true(all(!is.na(cov1[1:3, 1:3])))
  expect_true(all(is.na(cov1[4, ])))
  expect_true(all(is.na(cov1[, 4])))
})

test_that("vcov works for singular model matrix, complete=FALSE", {
  # that XtX is singular
  w <- cbind(mtcars$qsec, mtcars$qsec)
  model <- ivregranks(r(mpg) ~ r(hp) + w | r(disp) + w, data = mtcars)
  cov2 <- vcov(model, component = "stage2", complete = FALSE)
  cov1 <- vcov(model, component = "stage1", complete = FALSE)

  expect_true(all(!is.na(cov2)))
  expect_true(all(!is.na(cov1)))
  expect_equal(
    c(nrow(cov2), ncol(cov2)),
    c(3, 3)
  )
  expect_equal(
    c(nrow(cov1), ncol(cov1)),
    c(3, 3)
  )
})

#######################
### Low-level tests ###
#######################

test_that("update_coefficients_when_dropping_regressors works", {
  data(mtcars)
  model <- ivregranks(r(mpg) ~ r(hp) + cyl + drat | r(disp) + cyl + drat, data = mtcars)
  projection_matrix <- get_projection_residual_matrix_ivregranks(model)

  # expectation
  proj_1 <- lmranks(r(hp) ~ r(disp) + cyl + drat - 1, data = mtcars)
  proj_no_instrument <- lmranks(r(hp) ~ cyl + drat, data = mtcars)
  proj_2 <- lmranks(r(hp) ~ r(disp) + drat, data = mtcars)
  proj_3 <- lmranks(r(hp) ~ r(disp) + cyl, data = mtcars)

  expected <- matrix(c(
    0, coef(proj_1),
    coef(proj_no_instrument)[1], 0, coef(proj_no_instrument)[2:3],
    coef(proj_2)[1:2], 0, coef(proj_2)[3],
    coef(proj_3), 0
  ), ncol = 4, byrow = FALSE)

  actual <- update_coefficients_when_dropping_regressors(model, projection_matrix)

  expect_equal(actual, expected)
})

test_that("calculate_projection_residual_matrix_ivregranks works", {
  data(mtcars)
  model <- ivregranks(r(mpg) ~ r(hp) + cyl | r(disp) + cyl, data = mtcars)

  # expectation

  X_proj_1 <- lmranks(r(hp) ~ r(disp) + cyl - 1, data = mtcars)
  X_proj_2 <- lmranks(r(hp) ~ r(disp), data = mtcars)

  Z_proj <- lmranks(r(disp) ~ cyl, data = mtcars)
  intercept <- rep(1, nrow(mtcars))
  cyl_proj <- lm(cyl ~ X_proj_2$fitted.values, data = mtcars)
  intercept_proj <- lm(intercept ~ X_proj_1$fitted.values + cyl - 1, data = mtcars)

  expected <- matrix(c(
    1, -coef(intercept_proj),
    -coef(Z_proj)[1], 1, -coef(Z_proj)[2],
    -coef(cyl_proj), 1
  ), ncol = 3, byrow = FALSE)

  actual <- calculate_projection_residual_matrix_ivregranks(model)

  expect_equal(actual, expected)
})

######################################################
### High-level checks against by-hand calculations ###
######################################################

test_that("h1 works for ranked regressor with no covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- ivregranks(r(Y) ~ r(X) | r(Z))
  regressor_dropped_fs <- is.na(coef(res, component = "stage1"))
  U <- stats::model.matrix(res, component = "instruments")
  R <- qr.R(qr(U[, !regressor_dropped_fs]))
  proj_resid_matrix <- calculate_projection_residual_matrix(
    R,
    regressor_dropped_fs, length(coef(res, component = "stage1"))
  )
  proj_residuals <- U %*% proj_resid_matrix
  h1_ivregranks <- calculate_H1(res, proj_residuals)
  expect_equivalent(h1_ivregranks[, 2], h1)
})

test_that("h2 works for ranked regressor with no covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- ivregranks(r(Y) ~ r(X) | r(Z))
  regressor_dropped_fs <- is.na(coef(res, component = "stage1"))
  U <- stats::model.matrix(res, component = "instruments")
  R <- qr.R(qr(U[, !regressor_dropped_fs]))
  proj_resid_matrix <- calculate_projection_residual_matrix(
    R,
    regressor_dropped_fs, length(coef(res, component = "stage1"))
  )
  proj_residuals <- U %*% proj_resid_matrix
  h1_ivregranks <- calculate_H1(res, proj_residuals)
  h1_mean <- colMeans(h1_ivregranks)
  h2_ivregranks <- calculate_H2(res, proj_residuals, h1_mean)
  expect_equivalent(h2_ivregranks[, 2], h2)
})

test_that("h3 works for ranked regressor with no covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- ivregranks(r(Y) ~ r(X) | r(Z))
  regressor_dropped_fs <- is.na(coef(res, component = "stage1"))
  U <- stats::model.matrix(res, component = "instruments")
  R <- qr.R(qr(U[, !regressor_dropped_fs]))
  proj_resid_matrix <- calculate_projection_residual_matrix(
    R,
    regressor_dropped_fs, length(coef(res, component = "stage1"))
  )
  proj_residuals <- U %*% proj_resid_matrix
  h1_ivregranks <- calculate_H1(res, proj_residuals)
  h1_mean <- colMeans(h1_ivregranks)
  h3_ivregranks <- calculate_H3(res, proj_resid_matrix, h1_mean)
  expect_equivalent(h3_ivregranks[, 2], h3)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slow
  with no covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_FALSE.rda"))
  res <- ivregranks(r(Y) ~ r(X) | r(Z))
  sigma2hat_ivregranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)
})

test_that("h1 works for ranked regressor with covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)
  regressor_dropped_fs <- is.na(coef(res, component = "stage1"))
  U <- stats::model.matrix(res, component = "instruments")
  R <- qr.R(qr(U[, !regressor_dropped_fs]))
  proj_resid_matrix <- calculate_projection_residual_matrix(
    R,
    regressor_dropped_fs, length(coef(res, component = "stage1"))
  )
  proj_residuals <- U %*% proj_resid_matrix
  h1_ivregranks <- calculate_H1(res, proj_residuals)
  expect_equivalent(h1_ivregranks[, 2], h1)
})

test_that("h2 works for ranked regressor with covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)
  regressor_dropped_fs <- is.na(coef(res, component = "stage1"))
  U <- stats::model.matrix(res, component = "instruments")
  R <- qr.R(qr(U[, !regressor_dropped_fs]))
  proj_resid_matrix <- calculate_projection_residual_matrix(
    R,
    regressor_dropped_fs, length(coef(res, component = "stage1"))
  )
  proj_residuals <- U %*% proj_resid_matrix
  h1_ivregranks <- calculate_H1(res, proj_residuals)
  h1_mean <- colMeans(h1_ivregranks)
  h2_ivregranks <- calculate_H2(res, proj_residuals, h1_mean)
  expect_equivalent(h2_ivregranks[, 2], h2)
})

test_that("h3 works for ranked regressor with covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)
  regressor_dropped_fs <- is.na(coef(res, component = "stage1"))
  U <- stats::model.matrix(res, component = "instruments")
  R <- qr.R(qr(U[, !regressor_dropped_fs]))
  proj_resid_matrix <- calculate_projection_residual_matrix(
    R,
    regressor_dropped_fs, length(coef(res, component = "stage1"))
  )
  proj_residuals <- U %*% proj_resid_matrix
  h1_ivregranks <- calculate_H1(res, proj_residuals)
  h1_mean <- colMeans(h1_ivregranks)
  h3_ivregranks <- calculate_H3(res, proj_resid_matrix, h1_mean)
  expect_equivalent(h3_ivregranks[, 2], h3)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope
  with covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)
  sigma2hat_ivregranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope
  with smaller datasets", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_n_10.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)
  sigma2hat_ivregranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)

  load(test_path("testdata", "ivregranks_cov_sigmahat_n_50.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)
  sigma2hat_ivregranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)

  load(test_path("testdata", "ivregranks_cov_sigmahat_n_100.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)
  sigma2hat_ivregranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope
  with increasing=FALSE", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_increasing_FALSE.rda"))
  res <- ivregranks(r(Y, increasing = FALSE) ~ r(X, increasing = FALSE) + W |
    r(Z, increasing = FALSE) + W, omega = 1)
  sigma2hat_ivregranks <- vcov(res)[2, 2] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)
})
