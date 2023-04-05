test_that("summary does not raise errors", {
  mod <- lmranks(r(mpg) ~ r(cyl) + disp, data=mtcars)
  expect_silent(summary(mod))
})