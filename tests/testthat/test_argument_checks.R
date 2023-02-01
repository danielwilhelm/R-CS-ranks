V <- diag(1:10)
V_with_NAs <- diag(c(1:8, NA, 10))

expected_args <- list(x=c(1:8, 10),
                      Sigma = diag(c(1:8, 10)))

test_that("`x` argument is handled correctly by process_csranks_args",{
  expect_equal(process_csranks_args(1:10, V, na.rm=TRUE),
               list(x=1:10, Sigma=V))
  expect_error(process_csranks_args(data.frame(x=1:10),
                                    V, na.rm=TRUE))
  expect_error(process_csranks_args(paste0(1:10),
                                    V, na.rm=TRUE))
})

test_that("`Sigma` argument is handled correctly by process_csranks_args",{
  expect_error(process_csranks_args(1:10, data.frame(x=1:10),
                                    V))
  expect_error(process_csranks_args(1:10, diag(paste0(1:10))))
  expect_error(process_csranks_args(1:10, V[-1,]))
  expect_error(process_csranks_args(1:10, V[,-1]))
})

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

test_that("`x` argument is handled correctly by process_x_counts_argument",{
  expect_equal(process_x_counts_argument(1:10, na.rm=TRUE),
               1:10)
  expect_equal(process_x_counts_argument(c(NA, 2:10), na.rm=TRUE),
               2:10)
  expect_error(process_x_counts_argument(data.frame(x=1:10), na.rm=TRUE))
  expect_error(process_x_counts_argument(paste0(1:10), na.rm=TRUE))
  expect_error(process_x_counts_argument(c(-1,2:10), na.rm=TRUE))
  expect_error(process_x_counts_argument(c(3.14, 2:10), na.rm=TRUE))
  expect_error(process_x_counts_argument(c(NA, 2:10), na.rm=FALSE))
})

test_that("`indices` argument is handled correctly by process_indices_argument",{
  expect_equal(process_indices_argument(1:10, 10),
               1:10)
  expect_equal(process_indices_argument(c(NA, 2:10), 10),
               1:10)
  expect_error(process_indices_argument(data.frame(x=1:10), 10))
  expect_error(process_indices_argument(paste0(1:10), 10))
  expect_error(process_indices_argument(c(0,2:10), 10))
  expect_error(process_indices_argument(c(1:9, 11), 10))
  expect_error(process_indices_argument(c(3.14, 2:10), 10))
})