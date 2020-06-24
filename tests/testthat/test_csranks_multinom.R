context("CS for ranks based on multinomial data")

test_that("return value is of correct class and size", {

	for (multcorr in c("Bonferroni","Holm")) {
			res1 <- csranks_multinom(1:10, coverage=0.95, multcorr=multcorr, indices=NA)
			res2 <- csranks_multinom(1:10, coverage=0.95, multcorr=multcorr, indices=1)
			res3 <- csranks_multinom(1:10, coverage=0.95, multcorr=multcorr, indices=1:3)

			expect_is(res1$L, "integer")
			expect_is(res1$U, "integer")
			expect_is(res2$L, "integer")
			expect_is(res2$U, "integer")
			expect_is(res3$L, "integer")
			expect_is(res3$U, "integer")

			expect_equal(length(res1$L), 10)
			expect_equal(length(res1$U), 10)
			expect_equal(length(res2$L), 1)
			expect_equal(length(res2$U), 1)
			expect_equal(length(res3$L), 3)
			expect_equal(length(res3$U), 3)

			expect_false(any(is.na(res1$L)))
			expect_false(any(is.na(res1$U)))
			expect_false(any(is.na(res2$L)))
			expect_false(any(is.na(res2$U)))
			expect_false(any(is.na(res3$L)))
			expect_false(any(is.na(res3$U)))
	}
})


test_that("NAs are handled correctly", {
	expect_error(csranks_multinom(c(1:8,NA,2), coverage=0.95, indices=NA))
	expect_error(csranks_multinom(c(1:8,NA,NA), coverage=0.95, indices=NA))

	res <- csranks_multinom(c(1:8,NA,2), coverage=0.95, indices=NA, na.rm=TRUE)
	expect_is(res$L, "integer")
	expect_is(res$U, "integer")
	expect_equal(length(res$L), 9)
	expect_equal(length(res$U), 9)
	expect_false(any(is.na(res$L)))
	expect_false(any(is.na(res$U)))
})


test_that("lower and upper bounds are in the correct range of values", {
	x <- c(rmultinom(1,1000,1:20))

	for (cstype in c("two-sided","lower","upper")) {
		for (multcorr in c("Bonferroni","Holm")) {
				res1 <- csranks_multinom(x, coverage=0.95, cstype=cstype, multcorr=multcorr, indices=NA)
				res2 <- csranks_multinom(x, coverage=0.95, cstype=cstype, multcorr=multcorr, indices=1)
				res3 <- csranks_multinom(x, coverage=0.95, cstype=cstype, multcorr=multcorr, indices=1:3)
	
				expect_true(all(res1$L<=res1$U))
				expect_true(all(res2$L<=res2$U))
				expect_true(all(res3$L<=res3$U))
	
				expect_true(all(res1$L<=20 & res1$U<=20 & res1$L>=1 & res1$U>=1))
				expect_true(all(res2$L<=20 & res2$U<=20 & res2$L>=1 & res2$U>=1))			
				expect_true(all(res3$L<=20 & res3$U<=20 & res3$L>=1 & res3$U>=1))
	
				if (cstype=="lower") {
					expect_equal(res1$U, rep(20,20))
					expect_equal(res2$U, 20)
					expect_equal(res3$U, rep(20,3))
				}
				if (cstype=="upper") {
					expect_equal(res1$L, rep(1,20))
					expect_equal(res2$L, 1)
					expect_equal(res3$L, rep(1,3))
				}
			
		}
	}
})


test_that("simultaneous CS is not wider than marginal CS", {
	x <- c(rmultinom(1,1000,1:20))

	for (multcorr in c("Bonferroni","Holm")) {
			res1S <- csranks_multinom(x, coverage=0.95, simul=TRUE, multcorr=multcorr, indices=NA)
			res2S <- csranks_multinom(x, coverage=0.95, simul=TRUE, multcorr=multcorr, indices=1)
			res3S <- csranks_multinom(x, coverage=0.95, simul=TRUE, multcorr=multcorr, indices=1:3)
			res1M <- csranks_multinom(x, coverage=0.95, simul=FALSE, multcorr=multcorr, indices=NA)
			res2M <- csranks_multinom(x, coverage=0.95, simul=FALSE, multcorr=multcorr, indices=1)
			res3M <- csranks_multinom(x, coverage=0.95, simul=FALSE, multcorr=multcorr, indices=1:3)

			expect_true(all(res1M$L>=res1S$L & res1M$U<=res1S$U))
			expect_true(all(res2M$L>=res2S$L & res2M$U<=res2S$U))
			expect_true(all(res3M$L>=res3S$L & res3M$U<=res3S$U))
	}
})