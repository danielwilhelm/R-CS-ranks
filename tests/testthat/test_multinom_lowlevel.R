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
  expect_equal(which_pairs_are_compared(p, indices, "two-sided"),
               expected_symmetric_pairs)
  expect_equal(which_pairs_are_compared(p, indices, "lower"),
               t(expected_partial_pairs))
})

### calculate_pairwise_p_values ###
set.seed(2023)
x <- as.vector(rmultinom(1, 100, 1:p))
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

test_that("calculate_p_value returns sensible result", {
  p_values <- calculate_p_value(8:2, 10)
  expect_true(is.numeric(p_values))
  expect_length(p_values, 7)
  expect_true(all(p_values >= 0 & p_values <= 1))
  expect_true(!is.unsorted(p_values))
})

### reject_or_accept ###
df <- data.frame(j = 1:10, pval = c(1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2,
                                   1e-1, 3e-1))
unordered_df <- df[c(5, 1, 6, 4, 8, 2, 10, 3, 7, 9),]
coverage <- 0.95
expected_df_B <- df # effective beta == 0.005 == 5e-3
expected_df_B$rej <- c(rep(TRUE, 6), rep(FALSE, 4))
expected_df_H <- df # effective beta == 0.05 / (10:1)
# 0.005 0.006 0.006 0.007 0.008 0.010 0.015 0.016 0.025 0.05
expected_df_H$rej <- c(rep(TRUE, 7), rep(FALSE, 3))

test_that("reject_or_accept_works", {
  expect_equal(reject_or_accept(df, "Bonferroni", coverage),
               expected_df_B)
  expect_equal(reject_or_accept(unordered_df, "Bonferroni", coverage),
               expected_df_B)
  expect_equal(reject_or_accept(df, "Holm", coverage),
               expected_df_H)
  expect_equal(reject_or_accept(unordered_df, "Holm", coverage),
               expected_df_H)
})

# Edge cases:
# 1) reject all
df_edge_H <- df
df_edge_H$pval <- rep(1e-4, 10)
expected_df_edge_H <- df_edge_H
expected_df_edge_H$rej <- rep(TRUE, 10)
# 2) accept second from end despite p-value lower than 0.05 / 2
# (cause we accepted the third from end)
df_edge_2_H <- data.frame(j = 1:4,
                          pval = c(1e-5, 1.7e-2, 2e-2, 1e-1))
expected_df_edge_2_H <- df_edge_2_H
expected_df_edge_2_H$rej <- c(TRUE, FALSE, FALSE, FALSE)

test_that("reject_or_accept works for Holm in edge cases", {
  expect_equal(reject_or_accept(df_edge_H, "Holm", coverage),
               expected_df_edge_H)
  expect_equal(reject_or_accept(df_edge_2_H, "Holm", coverage),
               expected_df_edge_2_H)
})

### calculate_N_plus_minus ###

# Feature values:
# X1 < X2 == X3 < X4
df <- data.frame(
  j = c(1,1,1,2,2,2,3,3,3,4,4,4),
  k = c(2,3,4,1,3,4,1,2,4,1,2,3),
  rej = c(rep(FALSE, 3), 
          TRUE, FALSE, FALSE,
          TRUE, FALSE, FALSE, 
          TRUE, TRUE, TRUE)
)

unordered_df <- df[sample(nrow(df)),]

expected_symmetic <- list(Nminus = c(3,1,1,0),
                          Nplus = c(0,1,1,3))
expected_upper <- list(Nminus = c(0,0,0,0),
                       Nplus = c(0,1,1,3))
expected_lower <- list(Nminus = c(3,1,1,0),
                       Nplus = c(0,0,0,0))

test_that("calculate_N_plus_minus works", {
  expect_equal(calculate_N_plus_minus(df, cstype = "two-sided", indices = 1:4),
               expected_symmetic)
  expect_equal(calculate_N_plus_minus(unordered_df, cstype = "two-sided", indices = 1:4),
               expected_symmetic)
  expect_equal(calculate_N_plus_minus(df, cstype = "upper", indices = 1:4),
               expected_upper)
  expect_equal(calculate_N_plus_minus(unordered_df, cstype = "upper", indices = 1:4),
               expected_upper)
  expect_equal(calculate_N_plus_minus(df, cstype = "lower", indices = 1:4),
               expected_lower)
  expect_equal(calculate_N_plus_minus(unordered_df, cstype = "lower", indices = 1:4),
               expected_lower)
})

indices <- c(2,4)
expected_symmetic <- list(Nminus = c(1,0),
                          Nplus = c(1,3))
expected_lower <- list(Nminus = c(1,0),
                       Nplus = c(0,0))
expected_upper <- list(Nminus = c(0,0),
                       Nplus = c(1,3))

test_that("calculate_N_plus_minus works with indices argument", {
  expect_equal(calculate_N_plus_minus(df, cstype = "two-sided", indices = indices),
               expected_symmetic)
  expect_equal(calculate_N_plus_minus(df, cstype = "upper", indices = indices),
               expected_upper)
  expect_equal(calculate_N_plus_minus(df, cstype = "lower", indices = indices),
               expected_lower)
})








