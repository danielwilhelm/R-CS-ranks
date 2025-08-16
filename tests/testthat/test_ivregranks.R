test_that("ivregranks and by-hand calculation provide same results", {
  df <- mtcars[1:10, ]
  model <- ivregranks(r(mpg) ~ cyl | r(hp) | r(disp), data = df)
  expected_response <- c(0.6, 0.6, 0.9, 0.7, 0.3, 0.2, 0.1, 1.0, 0.9, 0.4)
  expected_model_matrix_regressors <- matrix(c(
    1, 0.7, 6,
    1, 0.7, 6,
    1, 0.2, 4,
    1, 0.7, 6,
    1, 0.9, 8,
    1, 0.4, 6,
    1, 1.0, 8,
    1, 0.1, 4,
    1, 0.3, 4,
    1, 0.8, 6
  ), byrow = TRUE, nrow = 10)
  expected_model_matrix_instruments <- matrix(c(
    1, 0.5, 6,
    1, 0.5, 6,
    1, 0.1, 4,
    1, 0.8, 6,
    1, 1.0, 8,
    1, 0.7, 6,
    1, 1.0, 8,
    1, 0.3, 4,
    1, 0.2, 4,
    1, 0.6, 6
  ), byrow = TRUE, nrow = 10)

  expected_coef_stage1 <- solve(t(expected_model_matrix_instruments) %*%
    expected_model_matrix_instruments) %*%
    t(expected_model_matrix_instruments) %*%
    expected_model_matrix_regressors[, 2]

  PX <- expected_model_matrix_instruments %*% expected_coef_stage1
  expected_model_matrix_projected <- expected_model_matrix_regressors
  expected_model_matrix_projected[, 2] <- PX
  expected_coef_stage2 <- solve(t(expected_model_matrix_projected) %*%
    expected_model_matrix_projected) %*%
    t(expected_model_matrix_projected) %*% expected_response

  expect_equivalent(model.response(model.frame(model)), expected_response)
  expect_equivalent(
    model.matrix(model, component = "instruments"),
    expected_model_matrix_instruments
  )
  expect_equivalent(
    model.matrix(model, component = "projected"),
    expected_model_matrix_projected
  )
  expect_equivalent(
    model.matrix(model, component = "regressors"),
    expected_model_matrix_regressors
  )

  expect_equivalent(coef(model, component = "stage1"), expected_coef_stage1)
  expect_equivalent(coef(model, component = "stage2"), expected_coef_stage2)
})

test_that("ivregranks and ivreg provide coherent results", {
  Y <- c(3, 1, 2, 4, 5)
  y_frank <- c(0.6, 0.2, 0.4, 0.8, 1.0)
  X <- 1:5
  x_frank <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  Z <- c(7, 6, 4, 5, 8)
  z_frank <- c(0.8, 0.6, 0.2, 0.4, 1.0)
  W <- c(1, 3, 2, 5, 4)
  omega <- 0.5

  rank_m <- ivregranks(r(Y) ~ W | r(X) | r(Z))
  raw_rank_m <- unclass(rank_m)
  raw_rank_m$call <- as.character(raw_rank_m$call)
  raw_rank_m$terms <- NULL
  attr(raw_rank_m$model, "terms") <- NULL
  raw_rank_m$omega <- NULL
  raw_rank_m$rank_terms_indices <- NULL
  raw_rank_m$ranked_instruments_indices <- NULL
  raw_rank_m$object_seqn <- NULL
  raw_rank_m$object_fs <- NULL

  m <- ivreg::ivreg(y_frank ~ W | x_frank | z_frank)
  expected_m <- unclass(m)
  expected_m$df.residual <- NA
  expected_m$formula <- formula("r(Y) ~ r(X) + W | r(Z) + W")
  expected_m$call <- as.character(str2lang(
    "ivregranks(r(Y) ~ W | r(X) | r(Z))"
  ))
  expected_m$ranked_response <- TRUE
  expected_m$terms <- NULL
  attr(expected_m$model, "terms") <- NULL
  names(expected_m$coefficients)[2] <- colnames(expected_m$residuals1)[2] <-
    colnames(expected_m$coefficients1)[2] <-
    dimnames(expected_m$qr$qr)[[2]][2] <- names(expected_m$endogenous) <-
    colnames(expected_m$cov.unscaled)[2] <-
    rownames(expected_m$cov.unscaled)[2] <- "r(X)"
  names(expected_m$instruments) <- rownames(expected_m$coefficients1)[2] <-
    dimnames(expected_m$qr1$qr)[[2]][2] <- "r(Z)"
  colnames(expected_m$model)[c(1, 2, 4)] <- c("r(Y)", "r(X)", "r(Z)")

  expect_equal(raw_rank_m, expected_m)
})

test_that("ivregranks falls back to ivreg in no rank case", {
  expect_warning(
    m <- ivregranks(mpg ~ cyl | hp | disp, data = mtcars),
    "no ranked terms"
  )
  m2 <- ivreg::ivreg(mpg ~ cyl | hp | disp, data = mtcars)

  expect_equivalent(m, m2)
})

test_that("ivregranks raises error if NA is encountered in data", {
  df1 <- mtcars
  df1[5, "disp"] <- NA
  expect_error(ivregranks(r(mpg) ~ r(hp) | disp, data = df1), "missing values")

  df2 <- mtcars
  df2[5, "cyl"] <- NA
  expect_error(ivregranks(r(mpg) ~ r(hp) + cyl | disp + cyl,
    data = df2
  ), "missing values")
})

test_that("process_ivregranks_formula catches illegal formulas", {
  expect_error(process_ivregranks_formula("y ~ x + w | z + w"))
  expect_error(process_ivregranks_formula(y ~ r(x) + r(w) | z + w))
  expect_error(process_ivregranks_formula(y ~ x + w | r(z) + r(w)))
  expect_error(process_ivregranks_formula(r(y) ~ r(x) * w | z))
  expect_error(process_ivregranks_formula(r(y) ~ r(x) | r(z) * w))
  expect_error(process_ivregranks_formula(r(y) ~ r(x) | r(z) + r(z):w))
  expect_error(process_ivregranks_formula(r(y) ~ r(x) | r(z):v:w))
  expect_error(process_ivregranks_formula(r(y) ~ r(x) | r(z):v + w))
  expect_error(process_ivregranks_formula(r(y) ~ w | r(x) + v | r(z1) + r(z2)))

  expect_silent(process_ivregranks_formula(r(y) ~ r(x) | r(z)))
  expect_silent(process_ivregranks_formula(r(y) ~ r(x) + w | r(z) + w))
  expect_silent(process_ivregranks_formula(r(y) ~ r(x) | r(z) + w))
  expect_silent(process_ivregranks_formula(r(y) ~ w | r(x) | r(z)))
  expect_silent(process_ivregranks_formula(r(y) ~ w | r(x) + v | r(z1) + z2))
  expect_silent(process_ivregranks_formula(r(y) ~ r(x) | r(z):G))
  expect_silent(process_ivregranks_formula(r(y) ~ r(x) | (r(z) + w):G))
  expect_silent(process_ivregranks_formula(r(y) ~ r(x) | r(z):G + G))
  expect_silent(process_ivregranks_formula(r(y) ~ r(x) | r(z):G - 1))
})

test_that("process_ivregranks_formula env to formula", {
  env <- new.env()

  actual <- process_ivregranks_formula(r(y) ~ w | r(x) | r(z), env)$formula
  expect_equal(environment(actual), env)

  actual <- process_ivregranks_formula(r(y) ~ r(x) | r(z):G, env)$formula
  expect_equal(environment(actual), env)

  actual <- process_ivregranks_formula(
    r(y) ~ r(x) | r(z):G - 1,
    env
  )$formula
  expect_equal(environment(actual), env)
})

test_that("process_ivregranks_formula returns correct index for simplest
  fits", {
  expect_equal(
    process_ivregranks_formula(r(y) ~ r(x) | r(z)),
    list(
      rank_terms_indices = 2, ranked_instruments_indices = 2,
      ranked_response = TRUE, formula = Formula::as.Formula(r(y) ~ r(x) | r(z))
    )
  )
  expect_equal(
    process_ivregranks_formula(r(y) ~ r(x) - 1 | r(z)),
    list(
      rank_terms_indices = 1, ranked_instruments_indices = 2,
      ranked_response = TRUE,
      formula = Formula::as.Formula(r(y) ~ r(x) - 1 | r(z))
    )
  )
  expect_equal(
    process_ivregranks_formula(r(y) ~ r(x) - 1 | r(z) - 1),
    list(
      rank_terms_indices = 1, ranked_instruments_indices = 1,
      ranked_response = TRUE,
      formula = Formula::as.Formula(r(y) ~ r(x) - 1 | r(z) - 1)
    )
  )
  expect_equal(
    process_ivregranks_formula(y ~ r(x) - 1 | r(z) - 1),
    list(
      rank_terms_indices = 1, ranked_instruments_indices = 1,
      ranked_response = FALSE,
      formula = Formula::as.Formula(y ~ r(x) - 1 | r(z) - 1)
    )
  )
})

test_that("prepare_ivreg_call works", {
  input_call <- str2lang("ivregranks(r(y) ~ r(x) + w | r(z) + w, data=data)")
  expected_call <-
    str2lang("ivreg::ivreg(r(y) ~ r(x) + w | r(z) + w, data=data,
      na.action=stats::na.fail)")
  expect_equal(prepare_ivreg_call(input_call), expected_call)

  input_call <- str2lang("ivregranks(r(y) ~ r(x) + w | r(z) + w, data=data,
    omega=omega)")
  expected_call <-
    str2lang("ivreg::ivreg(r(y) ~ r(x) + w | r(z) + w, data=data,
      na.action=stats::na.fail)")
  expect_equal(prepare_ivreg_call(input_call), expected_call)

  input_call <- str2lang("ivregranks(r(y) ~ r(x) + w | r(z) + w, data=data,
    na.rm=na.rm)")
  expected_call <-
    str2lang("ivreg::ivreg(r(y) ~ r(x) + w | r(z) + w, data=data,
      na.action=stats::na.fail)")
  expect_equal(prepare_ivreg_call(input_call), expected_call)
})

test_that("prepare_ivreg_call catches unsupported arguments", {
  input_call <- str2lang("ivregranks(r(y) ~ r(x) + w | r(z) + w, data=data,
    weights=weight)")
  expect_error(prepare_ivreg_call(input_call), "weights")

  input_call <- str2lang("ivregranks(r(y) ~ r(x) + w | r(z) + w, data=data,
    subset=x>0)")
  expect_error(prepare_ivreg_call(input_call), "subset")

  input_call <- str2lang("ivregranks(r(y) ~ r(x) + w | r(z) + w, data=data,
    na.action=na.action)")
  expect_error(prepare_ivreg_call(input_call), "na.action")
})
