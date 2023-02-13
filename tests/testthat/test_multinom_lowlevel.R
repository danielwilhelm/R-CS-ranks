### which_pairs_are_compared ###
p <- 5
indices <- c(2,4)
expected_partial_pairs <- matrix(c(FALSE, FALSE, FALSE, FALSE, FALSE,
                                TRUE, FALSE, TRUE, TRUE, TRUE,
                                FALSE, FALSE, FALSE, FALSE, FALSE,
                                TRUE, TRUE, TRUE, FALSE, TRUE,
                                FALSE, FALSE, FALSE, FALSE, FALSE),
                              byrow = TRUE, nrow = 5)

expected_symmetric_pairs <- matrix(c(FALSE, TRUE, FALSE, TRUE, FALSE,
                                  TRUE, FALSE, TRUE, TRUE, TRUE,
                                  FALSE, TRUE, FALSE, TRUE, FALSE,
                                  TRUE, TRUE, TRUE, FALSE, TRUE,
                                  FALSE, TRUE, FALSE, TRUE, FALSE),
                                byrow = TRUE, nrow = 5)

test_that("which_pairs_are_compared works", {
  expect_equal(which_pairs_are_compared(p, indices, "upper"),
               expected_partial_pairs)
  expect_equal(which_pairs_are_compared(p, indices, "symmetric"),
               expected_symmetric_pairs)
  expect_equal(which_pairs_are_compared(p, indices, "lower"),
               t(expected_partial_pairs))
})

### calculate_pairwise_p_values ###
set.seed(2023)
x <- rmultinom(p, 100, 1:p)
x_j <- x[c(2,4,4,2,4,2,2,4)]
S <- c(x[1] + x[2],
       x[1] + x[4],
       x[2] + x[4],
       x[3] + x[2],
       x[3] + x[4],
       x[4] + x[2],
       x[5] + x[2],
       x[5] + x[4])
expected_pvalues <- calculate_p_value(x_j, S)
expected_df <- data.frame(
  j = c(2,4,4,2,4,2,2,4),
  k = c(1,1,2,3,3,4,5,5),
  pval = expected_pvalues
)
test_that("calculate_pairwise_p_values works", {
  out <- calculate_pairwise_p_values(x, expected_partial_pairs)
  expect_equal(out, expected_df)
})













