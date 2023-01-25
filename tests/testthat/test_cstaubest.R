context("CS for tau-best")
V <- diag(rep(1, 10))
test_that("return value is of correct class and size", {
  for (stepdown in c(TRUE, FALSE)) {
    res1 <- cstaubest(1:10, V, tau = 2, coverage = 0.95, stepdown = stepdown, R = 100)
    res2 <- cstaubest(1:10, V, tau = 5, coverage = 0.95, stepdown = stepdown, R = 100)

    expect_is(res1, "logical")
    expect_is(res2, "logical")

    expect_equal(length(res1), 10)
    expect_equal(length(res2), 10)

    expect_false(any(is.na(res1)))
    expect_false(any(is.na(res2)))
  }
})


test_that("NAs are handled correctly", {
  expect_error(cstaubest(c(1:8, NA, 2), V, tau = 2, coverage = 0.95, stepdown = FALSE, R = 100))
  expect_error(cstaubest(c(1:8, NA, NA), V, tau = 2, coverage = 0.95, stepdown = FALSE, R = 100))

  res <- cstaubest(c(1:8, NA, 2), V, tau = 2, coverage = 0.95, stepdown = FALSE, R = 100, na.rm = TRUE)
  expect_is(res, "logical")
  expect_equal(length(res), 9)
  expect_false(any(is.na(res)))

  V_with_NAs <- diag(c(rep(1, 8), NA, 1))
  res <- cstaubest(1:10, V_with_NAs, tau = 2, coverage = 0.95, stepdown = FALSE, R = 100, na.rm = TRUE)
  expect_is(res, "logical")
  expect_equal(length(res), 9)
  expect_false(any(is.na(res)))
  
  res <- cstaubest(1:4, matrix(c(1,NA,0.1,0.1,
                                  0.1,1,0.1,0.1,
                                  0.1,0.1,1,0.1,
                                  0.1,0.1,0.1,1),4), tau = 1,
                      coverage = 0.95, stepdown = FALSE, R = 100, na.rm = TRUE)
  expect_is(res, "logical")
  expect_equal(length(res), 2)
  expect_false(any(is.na(res)))
})


test_that("correct behavior for varying tau", {
  for (stepdown in c(TRUE, FALSE)) {
    res1 <- cstaubest(1:10, V, tau = 2, coverage = 0.95, stepdown = stepdown, R = 100, seed = 1)
    res2 <- cstaubest(1:10, V, tau = 5, coverage = 0.95, stepdown = stepdown, R = 100, seed = 1)

    expect_true(identical(res1 & res2, res1))
    expect_true(sum(res2) >= sum(res1))
  }
})
