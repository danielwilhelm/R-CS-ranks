context("CS for tau-worst")
V <- diag(rep(1, 10))
test_that("return value is of correct class and size", {
  for (stepdown in c(TRUE, FALSE)) {
    res1 <- cstauworst(1:10, V, tau = 2, coverage = 0.95, stepdown = stepdown, R = 100)
    res2 <- cstauworst(1:10, V, tau = 5, coverage = 0.95, stepdown = stepdown, R = 100)

    expect_is(res1, "logical")
    expect_is(res2, "logical")

    expect_equal(length(res1), 10)
    expect_equal(length(res2), 10)

    expect_false(any(is.na(res1)))
    expect_false(any(is.na(res2)))
  }
})

test_that("correct behavior for varying tau", {
  for (stepdown in c(TRUE, FALSE)) {
    res1 <- cstauworst(1:10, V, tau = 2, coverage = 0.95, stepdown = stepdown, R = 100, seed = 1)
    res2 <- cstauworst(1:10, V, tau = 5, coverage = 0.95, stepdown = stepdown, R = 100, seed = 1)

    expect_true(identical(res1 & res2, res1))
    expect_true(sum(res2) >= sum(res1))
  }
})
