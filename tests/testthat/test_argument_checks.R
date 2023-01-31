V <- diag(1:10)
V_with_NAs <- diag(c(1:8, NA, 10))

expected_args <- list(x=c(1:8, 10),
                      Sigma = diag(c(1:8, 10)))
test_that("NAs are handled correctly by process_csranks_args", {
  expect_error(process_csranks_args(c(1:8, NA, 10),
                                    V, na.rm=FALSE))
  expect_error(process_csranks_args(c(1:10),
                                    V_with_NAs, na.rm=FALSE))
  expect_equal(process_csranks_args(c(1:8, NA, 10),
                                    V, na.rm=TRUE),
               expected_args)
  expect_equal(process_csranks_args(1:10,
                                    V_with_NAs, na.rm=TRUE),
               expected_args)
})

test_that("NAs are handled correctly at the higher level", {
  expect_error(csranks(c(1:8, NA, 10), V, na.rm=FALSE))
  expect_error(csranks(c(1:10), V_with_NAs, na.rm=FALSE))
  out_handled_NAs_internally <- csranks(c(1:8, NA, 10), V, na.rm=TRUE, seed=2023)
  out_handled_NAs_manually <- csranks(c(1:8, 10), diag(c(1:8, 10)), seed=2023)
  expect_equal(out_handled_NAs_internally, out_handled_NAs_manually)
  out_handled_NAs_internally <- csranks(1:10, V_with_NAs, na.rm=TRUE, seed=2023)
  expect_equal(out_handled_NAs_internally, out_handled_NAs_manually)
})