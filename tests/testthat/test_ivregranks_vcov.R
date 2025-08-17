test_that("summary does not raise errors", {
  model <- ivregranks(r(mpg) ~ r(hp) + cyl | r(disp) + cyl, data = mtcars)
  expect_silent(summary(model))
})
