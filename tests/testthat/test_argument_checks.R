V <- diag(1:10)
V_with_NAs <- diag(c(1:8, NA, 10))

expected_args <- list(x = c(1:8, 10),
                      Sigma = diag(c(1:8, 10)),
                      indices = 1:9)

test_that("`x` argument is handled correctly by process_csranks_args",{
  expect_equal(process_csranks_args(1:10, V, NA, na.rm=TRUE),
               list(x=1:10, Sigma=V, indices=1:10))
  expect_error(process_csranks_args(data.frame(x=1:10),
                                    V, NA, na.rm=TRUE))
  expect_error(process_csranks_args(paste0(1:10),
                                    V, NA, na.rm=TRUE))
})

test_that("`Sigma` argument is handled correctly by process_csranks_args",{
  expect_error(process_csranks_args(1:10, data.frame(x=1:10), NA, na.rm=TRUE))
  expect_error(process_csranks_args(1:10, 1:10 / 10, NA, na.rm=TRUE),
               "variances")
  expect_error(process_csranks_args(1:10, V[-1,], NA, na.rm=TRUE))
  expect_error(process_csranks_args(1:10, V[,-1], NA, na.rm=TRUE))
})

test_that("NAs are handled correctly by process_csranks_args", {
  expect_error(process_csranks_args(c(1:8, NA, 10),
                                    V, NA, na.rm=FALSE))
  expect_error(process_csranks_args(c(1:10),
                                    V_with_NAs, NA, na.rm=FALSE))
  expect_equal(process_csranks_args(c(1:8, NA, 10),
                                    V, NA, na.rm=TRUE),
               expected_args)
  expect_equal(process_csranks_args(1:10,
                                    V_with_NAs, NA, na.rm=TRUE),
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

test_that("adjust_indices_for_NAs work correctly",{
  expect_equal(adjust_indices_for_NAs(1:10, rep(TRUE, 10)),
               1:10)
  expect_equal(adjust_indices_for_NAs(1:10, c(rep(TRUE, 4),
                                              rep(FALSE, 3), 
                                              rep(TRUE, 3))),
               1:7)
  expect_equal(adjust_indices_for_NAs(seq(2,10,2), c(TRUE, FALSE, TRUE, FALSE, FALSE,
                                              TRUE, TRUE, TRUE, FALSE, TRUE)),
               c(3,5,6))
})

test_that("`x` argument is handled correctly by process_csranks_multinom_argss",{
  expect_equal(process_csranks_multinom_args(1:10, NA, na.rm=TRUE)[["x"]],
               1:10)
  expect_equal(process_csranks_multinom_args(c(NA, 2:10), NA, na.rm=TRUE)[["x"]],
               2:10)
  expect_error(process_csranks_multinom_args(data.frame(x=1:10), NA, na.rm=TRUE))
  expect_error(process_csranks_multinom_args(paste0(1:10), NA, na.rm=TRUE))
  expect_error(process_csranks_multinom_args(c(-1,2:10), NA, na.rm=TRUE))
  expect_error(process_csranks_multinom_args(c(3.14, 2:10), NA, na.rm=TRUE))
  expect_error(process_csranks_multinom_args(c(NA, 2:10), NA, na.rm=FALSE))
})

test_that("`indices` argument is handled correctly by process_csranks_multinom_args",{
  expect_equal(process_csranks_multinom_args(1:10, NA, na.rm=TRUE)[["indices"]],
               1:10)
  expect_equal(process_csranks_multinom_args(c(NA, 2:10), NA, na.rm=TRUE)[["indices"]],
               1:9)
  expect_equal(process_csranks_multinom_args(c(NA, 2:10), 1:10, na.rm=TRUE)[["indices"]],
               1:9)
  expect_equal(process_csranks_multinom_args(c(NA, 2:10), 2:10, na.rm=TRUE)[["indices"]],
               1:9)
  expect_equal(process_csranks_multinom_args(c(NA, 2,3,NA,NA,6:10), 
                                         c(1,3,5,6), na.rm=TRUE)[["indices"]],
               c(2,3))
})

test_that("`indices` argument is handled correctly by process_indices_argument",{
  expect_equal(process_indices_argument(1:10, 10),
               1:10)
  expect_equal(process_indices_argument(c(NA, 2:10), 10),
               1:10)
  expect_equal(process_indices_argument(c(-2,-5,-6), 10),
               c(1,3,4,7,8,9,10))
  expect_error(process_indices_argument(data.frame(x=1:10), 10))
  expect_error(process_indices_argument(paste0(1:10), 10))
  expect_error(process_indices_argument(c(0,2:10), 10))
  expect_error(process_indices_argument(c(1:9, 11), 10))
  expect_error(process_indices_argument(c(3.14, 2:10), 10))
  expect_error(process_indices_argument(c(0, 2:10), 10))
  expect_error(process_indices_argument(c(-1, 2:10), 10))
  expect_error(process_indices_argument(2:11, 10))
})
