test_that("lmranks correctly estimates OLS slope and intercept in grouped case", {
  set.seed(100)
  
  for (dist in c("discr", "cont")) {

  	if (dist == "discr") {
  		X <- rbinom(200, 5, 0.5)
  		Y <- c(rep(1,100),rep(2,100))*X + rbinom(200, 2, 0.5)
  	} else {
  		X <- rnorm(200)
		Y <- c(rep(1,100),rep(2,100))*X + rnorm(200)
  	}
  	G <- factor(c(rep(1,100),rep(2,100)))
  	dat <- data.frame(Y=Y, X=X)
	  	  
	for (omega in c(0, 0.5, 1)) {
	    RY <- frank(Y, omega=omega, increasing=TRUE)
	    RX <- frank(X, omega=omega, increasing=TRUE)
	    RY1 <- RY[1:100]; RX1 <- RX[1:100]
	    RY2 <- RY[101:200]; RX2 <- RX[101:200]
	    coeffs.lm <- as.numeric(cbind(coef(lm(RY1~RX1)), coef(lm(RY2~RX2))))

	    res <- lmranks(r(Y) ~ r(X):G, data=dat, omega=omega)
	    coeffs.grouped_lmranks <- as.numeric(coef(res)[c(1,3,2,4)])
	    expect_equal(coeffs.lm, coeffs.grouped_lmranks)     
	}

  }
 
})

test_that("lmranks catches non-factor groupings", {
  X <- rnorm(200)
  Y <- c(rep(1,100),rep(2,100))*X + rnorm(200)
  dat <- data.frame(Y=Y, X=X)
  G <- c(rep(1,100),rep(2,100))
  
  expect_error(lmranks(r(Y) ~ r(X):G, data=dat), "factor")
  
  G <- matrix(c(rep(1,100),rep(2,100),
                rep(1,100),rep(2,100)), ncol=2)
  expect_error(lmranks(r(Y) ~ r(X):G, data=dat), "factor")
  
  G <- factor(c(rep(1,100),rep(2,100)))
  expect_silent(lmranks(r(Y) ~ r(X):G, data=dat))
})

test_that("lmranks catches wrong contrasts", {
  X <- rnorm(200)
  Y <- c(rep(1,100),rep(2,100))*X + rnorm(200)
  dat <- data.frame(Y=Y, X=X)
  G <- factor(c(rep(1,100),rep(2,100)))
  
  expect_error(lmranks(r(Y) ~ r(X):G, data=dat, contrasts = list(G = "contr.sum")), 
               "contrasts")
})