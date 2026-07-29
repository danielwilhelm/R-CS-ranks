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
  model <- expect_warning(ivregranks(r(mpg) ~ r(hp) + w | r(disp) + w, data = mtcars), "collinear")
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

test_that("get_projection_residual_matrix_ivregranks works", {
  data(mtcars)
  model <- ivregranks(r(mpg) ~ r(hp) + cyl + drat | r(disp) + cyl + drat, data = mtcars)
  intercept <- rep(1, nrow(mtcars))

  # expectation
  proj_1 <- lmranks(intercept ~ r(disp) + cyl + drat - 1, data = mtcars)
  proj_of_instrument <- lmranks(r(disp) ~ cyl + drat, data = mtcars)
  proj_2 <- lmranks(cyl ~ r(disp) + drat, data = mtcars)
  proj_3 <- lmranks(drat ~ r(disp) + cyl, data = mtcars)

  expected <- matrix(c(
    1, -coef(proj_1),
    -coef(proj_of_instrument)[1], 1, -coef(proj_of_instrument)[2:3],
    -coef(proj_2)[1:2], 1, -coef(proj_2)[3],
    -coef(proj_3), 1
  ), ncol = 4, byrow = FALSE)

  actual <- get_projection_residual_matrix_ivregranks(model, "stage1")

  expect_equal(actual, expected)
})

test_that("update_endogenous_coefficients_when_dropping_instrument works", {
  data(mtcars)
  model <- ivregranks(r(mpg) ~ r(hp) + cyl + drat | r(disp) + cyl + drat, data = mtcars)
  projection_matrix <- get_projection_residual_matrix_ivregranks(model, "stage1")

  # expectation
  proj_no_instrument <- lmranks(r(hp) ~ cyl + drat, data = mtcars)

  expected <- coef(proj_no_instrument)

  actual <- update_endogenous_coefficients_when_dropping_instrument(model, projection_matrix)

  expect_equal(actual, expected)
})

test_that("substitute_coefs_change_base_to_stage_1 works", {
  data(mtcars)
  model <- ivregranks(r(mpg) ~ cyl + r(hp) | cyl + r(disp), data = mtcars)
  proj_resid <- get_projection_residual_matrix_ivregranks(model, "stage1")
  instrument_index <- 3

  # The actual inserted matrix into substitute_coefs_change_base_to_stage_1 is different,
  # but the behaviour should be the same

  substituted_proj_resid <- substitute_coefs_change_base_to_stage_1(model, proj_resid)
  expected <- model.matrix(model, "projected") %*% proj_resid[, -instrument_index]
  actual <- model.matrix(model, "instruments") %*% substituted_proj_resid[, -instrument_index]

  expect_equal(actual, expected)
  expect_equivalent(substituted_proj_resid[, instrument_index], proj_resid[, instrument_index])
})


test_that("get_projection_variances works", {
  data(mtcars)
  model <- ivregranks(r(mpg) ~ r(hp) + cyl | cyl + r(disp), data = mtcars)

  X_on_W <- 1:5
  residuals_rest <- matrix(seq(0, 1, length.out = 15), ncol = 3)

  expected <- c(
    mean(residuals_rest[, 1]^2),
    mean(residuals_rest[, 3] * X_on_W),
    mean(residuals_rest[, 2]^2)
  )
  actual <- get_projection_variances(model, X_on_W, residuals_rest)

  expect_equivalent(actual, expected)
})

test_that("get_instrument_index_after_dropping_NAs works for simplest case", {
  data(mtcars)
  model <- ivregranks(r(mpg) ~ r(hp) + cyl | cyl + r(disp), data = mtcars)
  expected <- 3
  actual <- get_instrument_index_after_dropping_NAs(model)
  expect_equal(actual, expected)
})

test_that("get_instrument_index_after_dropping_NAs works for case without intercept", {
  model_2 <- ivregranks(r(mpg) ~ r(hp) + cyl - 1 | cyl + r(disp) - 1, data = mtcars)
  expected <- 2
  actual <- get_instrument_index_after_dropping_NAs(model_2)
  expect_equal(actual, expected)
})

test_that("get_instrument_index_after_dropping_NAs works for reordered terms", {
  model_3 <- ivregranks(r(mpg) ~ r(hp) + cyl | r(disp) + cyl, data = mtcars)
  expected <- 2
  actual <- get_instrument_index_after_dropping_NAs(model_3)
  expect_equal(actual, expected)
})

test_that("get_instrument_index_after_dropping_NAs works for * term", {
  model_4 <- ivregranks(r(mpg) ~ r(hp) + cyl * wt | cyl * wt + r(disp), data = mtcars)
  expected <- 4
  actual <- get_instrument_index_after_dropping_NAs(model_4)
  expect_equal(actual, expected)
})

test_that("get_instrument_index_after_dropping_NAs works for interaction with factor", {
  cyl_f <- as.factor(mtcars$cyl)
  model_5 <- ivregranks(r(mpg) ~ r(hp) + cyl_f * wt | cyl_f * wt + r(disp), data = mtcars)
  expected <- 5
  actual <- get_instrument_index_after_dropping_NAs(model_5)
  expect_equal(actual, expected)
})

test_that("get_instrument_index_after_dropping_NAs works for colinear case", {
  W <- mtcars$cyl
  expect_warning(model_6 <- ivregranks(r(mpg) ~ r(hp) + cyl + W | cyl + W + r(disp), data = mtcars), "collinear")
  expected <- 3
  actual <- get_instrument_index_after_dropping_NAs(model_6)
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
  res <- ivregranks(r(Y) ~ r(X) + W - 1 | r(Z) + W - 1)
  regressor_dropped_fs <- is.na(coef(res, component = "stage1"))
  U <- stats::model.matrix(res, component = "instruments")
  R <- qr.R(qr(U[, !regressor_dropped_fs]))
  proj_resid_matrix <- calculate_projection_residual_matrix(
    R,
    regressor_dropped_fs, length(coef(res, component = "stage1"))
  )
  proj_residuals <- U %*% proj_resid_matrix
  h1_ivregranks <- calculate_H1(res, proj_residuals)
  expect_equivalent(h1_ivregranks[, 1], h1)
})

test_that("h2 works for ranked regressor with covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W - 1 | r(Z) + W - 1)
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
  expect_equivalent(h2_ivregranks[, 1], h2)
})

test_that("h3 works for ranked regressor with covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W - 1 | r(Z) + W - 1)
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
  expect_equivalent(h3_ivregranks[, 1], h3)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope
  with covariates", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_covariates_TRUE.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W - 1 | r(Z) + W - 1)
  sigma2hat_ivregranks <- vcov(res)[1, 1] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope
  with smaller datasets", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_n_10.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W - 1 | r(Z) + W - 1)
  sigma2hat_ivregranks <- vcov(res)[1, 1] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)

  load(test_path("testdata", "ivregranks_cov_sigmahat_n_50.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W - 1 | r(Z) + W - 1)
  sigma2hat_ivregranks <- vcov(res)[1, 1] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)

  load(test_path("testdata", "ivregranks_cov_sigmahat_n_100.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W - 1 | r(Z) + W - 1)
  sigma2hat_ivregranks <- vcov(res)[1, 1] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)
})

test_that("vcov produces correct asymptotic variance estimate of rank-rank slope
  with increasing=FALSE", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_increasing_FALSE.rda"))
  res <- ivregranks(r(Y, increasing = FALSE) ~ r(X, increasing = FALSE) + W - 1 |
    r(Z, increasing = FALSE) + W - 1, omega = 1)
  sigma2hat_ivregranks <- vcov(res)[1, 1] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)
})

test_that("vcov produces correct asymptotic variance estimate of regressor variance", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_regressor_1.rda"))
  res <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)
  sigma2hat_ivregranks <- vcov(res)[3, 3] * n
  expect_equal(sigma2hat, sigma2hat_ivregranks)
})

test_that("h1 works for regressor variance estimation", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_regressor_1.rda"))
  object <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)

  regressor_dropped_fs <- is.na(coef(res, component = "stage1"))
  projection_residual_matrix_stage_2 <- get_projection_residual_matrix_ivregranks(object, "stage2")
  projection_residual_matrix_stage_1 <- get_projection_residual_matrix_ivregranks(object, "stage1")
  instrument_index <- get_instrument_index_after_dropping_NAs(object)
  projection_residual_matrix_stage_2[, instrument_index] <- projection_residual_matrix_stage_1[, instrument_index]
  projection_residual_matrix_stage_2_in_terms_stage_1 <- substitute_coefs_change_base_to_stage_1(
    object,
    projection_residual_matrix_stage_2
  )

  Z <- stats::model.matrix(object, component = "instruments")[, !regressor_dropped_fs, drop = FALSE]
  projection_residuals_fs <- Z %*% projection_residual_matrix_stage_2_in_terms_stage_1

  h1_ivregranks <- calculate_H1(object, projection_residuals_fs)
  expect_equivalent(h1_ivregranks[, 3], h1)
})

test_that("h2 works for regressor variance estimation", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_regressor_1.rda"))
  object <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)

  regressor_dropped_fs <- is.na(coef(object, component = "stage1"))
  projection_residual_matrix_stage_2 <- get_projection_residual_matrix_ivregranks(object, "stage2")
  projection_residual_matrix_stage_1 <- get_projection_residual_matrix_ivregranks(object, "stage1")
  instrument_index <- get_instrument_index_after_dropping_NAs(object)
  projection_residual_matrix_stage_2[, instrument_index] <- projection_residual_matrix_stage_1[, instrument_index]
  projection_residual_matrix_stage_2_in_terms_stage_1 <- substitute_coefs_change_base_to_stage_1(
    object,
    projection_residual_matrix_stage_2
  )
  Z <- stats::model.matrix(object, component = "instruments")[, !regressor_dropped_fs]
  projection_residuals_fs <- Z %*% projection_residual_matrix_stage_2_in_terms_stage_1

  H1 <- calculate_H1(object, projection_residuals_fs)
  H1_mean <- colMeans(H1)
  H2 <- calculate_H2(object, projection_residuals_fs, H1_mean)
  expect_equivalent(H2[, 3], h2)
})

test_that("h3 works for regressor variance estimation", {
  load(test_path("testdata", "ivregranks_cov_sigmahat_regressor_1.rda"))
  object <- ivregranks(r(Y) ~ r(X) + W | r(Z) + W)

  regressor_dropped_fs <- is.na(coef(object, component = "stage1"))
  projection_residual_matrix_stage_2 <- get_projection_residual_matrix_ivregranks(object, "stage2")
  projection_residual_matrix_stage_1 <- get_projection_residual_matrix_ivregranks(object, "stage1")
  instrument_index <- get_instrument_index_after_dropping_NAs(object)
  projection_residual_matrix_stage_2[, instrument_index] <- projection_residual_matrix_stage_1[, instrument_index]
  projection_residual_matrix_stage_2_in_terms_stage_1 <- substitute_coefs_change_base_to_stage_1(
    object,
    projection_residual_matrix_stage_2
  )
  Z <- stats::model.matrix(object, component = "instruments")[, !regressor_dropped_fs]
  projection_residuals_fs <- Z %*% projection_residual_matrix_stage_2_in_terms_stage_1

  H1 <- calculate_H1(object, projection_residuals_fs)
  H1_mean <- colMeans(H1)
  H3 <- calculate_H3(object, projection_residual_matrix_stage_2_in_terms_stage_1, H1_mean)
  expect_equivalent(H3[, 3], h3)
})

test_that("vcov produces identical asymptotic variance estimate to simulations from dwilhelm's guys", {
  load(test_path("testdata", "ivregranks_vcov_sims.rda"))
  res <- ivregranks(r(Y) ~ r(X) | r(Z), data = df)
  sigma2hat_ivregranks <- vcov(res)
  expect_equal(var_est[2, 2], sigma2hat_ivregranks[2, 2])
})
