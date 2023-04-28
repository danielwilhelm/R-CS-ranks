test_that("vcov produces correct asymptotic variance estimate of rank-rank slope", {
  set.seed(100)

  for (covariates in c(TRUE,FALSE)) {
    
    # draw data
    n <- 10000
    G <- c(rep(1,n/2),rep(2,n/2))
    beta <- G
    if (covariates) {
      X <- rnorm(n)
      W <- matrix(rnorm(n*2), n, 2)
      Y <- beta/2 + beta*X + rowSums(W) + rnorm(n,0,0.5)  
      W <- cbind(1,W)  
    } else {
      X <- rnorm(n)
      Y <- beta/2 + beta*X + rnorm(n,0,0.5) 
      W <- matrix(1,n,1)
    }

    # compute ranks
    RY <- frank(Y, increasing=TRUE)
    RX <- frank(X, increasing=TRUE)
    RY1 <- RY[1:(n/2)]; RX1 <- RX[1:(n/2)]; W1 <- W[1:(n/2)]
    RY2 <- RY[(n/2+1):n]; RX2 <- RX[(n/2+1):n]; W2 <- W[(n/2+1):n]


    # ------- compute asymptotic variance "by hand"
      
      Ifn <- function(u, v) return( u<=v )

      # first stage
      res1 <- lm(RX1 ~ W1-1)
      Wgammahat1 <- predict(res1)
      nuhat1 <- resid(res1)
      gammahat1 <- coef(res1)

      # outcome equation
      res <- lm(RY~RX+W-1)
      rhohat <- coef(res)[1]
      betahat <- coef(res)[-1]
      epsilonhat <- resid(res)

      # construct h1
      h1 <- epsilonhat * nuhat

      # construct h2
      h2fn <- function(xy) mean((Ifn(xy[2],Y)-rhohat*Ifn(xy[1],X)-c(W%*%betahat)) * nuhat)
      h2 <- apply(cbind(X,Y), 1, h2fn)

      # construct h3
      h3fn <- function(x) mean(epsilonhat * (Ifn(x,X)-Wgammahat))
      h3 <- sapply(X, h3fn)

      # compute asymptotic variance
      sigma2hat <- mean((h1+h2+h3)^2) / var(nuhat)^2


    # ------- compute asymptotic variance using lmranks

      if (covariates) { 
        res <- lmranks(r(Y) ~ r(X) + W - 1)
        sigma2hat.lmranks <- vcov(res)[1,1]*n
      } else {
        res <- lmranks(r(Y) ~ r(X))
        sigma2hat.lmranks <- vcov(res)[2,2]*n
      }

      
    # ------- test equality
    
      expect_equal(sigma2hat, sigma2hat.lmranks, tolerance=1e-5)
  }
})