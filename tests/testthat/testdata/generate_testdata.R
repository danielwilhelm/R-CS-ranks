library(lmranks)
set.seed(100)

for (covariates in c(TRUE,FALSE)) {
  
  # draw data
  n <- 10000
  if (covariates) {
    X <- rnorm(n)
    W <- matrix(rnorm(n*2), n, 2)
    Y <- X + rowSums(W) + rnorm(n,0,0.5)  
    W <- cbind(1,W)  
  } else {
    X <- rnorm(n)
    Y <- X + rnorm(n,0,0.5) 
    W <- matrix(1,n,1)
  }
  
  # compute ranks
  RY <- frank(Y, increasing=TRUE)
  RX <- frank(X, increasing=TRUE)
  
  
  # ------- compute asymptotic variance "by hand"
  
  Ifn <- function(u, v) return( u<=v )
  
  # first stage
  res <- lm(RX ~ W-1)
  Wgammahat <- predict(res)
  nuhat <- resid(res)
  gammahat <- coef(res)
  
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
  
  # save the result
  save(sigma2hat, Y, W, X, n, file = file.path("tests", "testthat", "testdata",
                                   paste0("lmranks_cov_sigmahat_covariates_", covariates, ".rda")))
}

############################
### Hoeffding's variance ###
############################

set.seed(100)

# draw data
n <- 1000
X <- rnorm(n)
Y <- X + rnorm(n,0,0.5) 

# ------- compute Hoeffding's variance "by hand"

n <- length(Y)
rhohat <- cor(frank(Y, omega=1/2, increasing=TRUE),frank(X, omega=1/2, increasing=TRUE))
oX <- outer(X,X,'<=')
oY <- outer(Y,Y,'<=')
FhatX <- colMeans(oX)
FhatY <- colMeans(oY)
FhatAvgX <- function(x) mean(colMeans(outer(X,rep(x,n),'<=')*oY))
FhatAvgY <- function(y) mean(colMeans(outer(Y,rep(y,n),'<=')*oX))
psihatX <- sapply(X, FhatAvgX) - FhatX*mean(FhatY)
psihatY <- sapply(Y, FhatAvgY) - mean(FhatX)*FhatY

W <- 3*( (2*FhatX-1)*(2*FhatY-1) + 4*psihatX + 4*psihatY )

# compute asymptotic variance
sigma2hat <- var(W)

save(sigma2hat, Y, W, X, n, file = file.path("tests", "testthat", "testdata",
                                             "lmranks_cov_sigmahat_Hoefding.rda"))
