# setup
n <- 10 # not larger than 23
x <- seq(1,3,length=n)
sd <- rep(0.2,n)
ranks <- irank(x)
CS <- csranks(x, sd)
popnames <- rev(LETTERS[1:n])

test_that("default plotranking returns a ggplot object",{
  expect_s3_class(plotranking(ranks, CS$L, CS$U), "ggplot")
})

test_that("custom plotranking returns a ggplot object",{
  expect_s3_class(plotranking(ranks, CS$L, CS$U, popnames = popnames), "ggplot")
  expect_s3_class(plotranking(ranks, CS$L, CS$U, title = "title"), "ggplot")
  expect_s3_class(plotranking(ranks, CS$L, CS$U, subtitle = "subtitle"), "ggplot")
  expect_s3_class(plotranking(ranks, CS$L, CS$U, caption = "caption"), "ggplot")
  expect_s3_class(plotranking(ranks, CS$L, CS$U, colorbins = 2), "ggplot")
  expect_s3_class(plotranking(ranks, CS$L, CS$U, horizontal = FALSE), "ggplot")
})

