library(csranks)
set.seed(100)

####################
### vcov.lmranks ###
####################

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

############################
### vcov.grouped_lmranks ###
############################

set.seed(100)

for (covariates in c(TRUE,FALSE)) {
  
  # draw data
  n <- 1000
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
  RY <- frank(Y, increasing=TRUE, omega=1)
  RX <- frank(X, increasing=TRUE, omega=1)
  RY1 <- RY[1:(n/2)]; RX1 <- RX[1:(n/2)]; W1 <- W[1:(n/2),]
  RY2 <- RY[(n/2+1):n]; RX2 <- RX[(n/2+1):n]; W2 <- W[(n/2+1):n,]
  
  
  # ------- compute asymptotic variance "by hand"
  
  Ifn <- function(u, v) return( u<=v )
  
  # first stage
  res1 <- lm(RX1 ~ W1-1)
  Wgammahat1 <- predict(res1)
  nuhat1 <- resid(res1)
  gammahat1 <- coef(res1)
  res2 <- lm(RX2 ~ W2-1)
  Wgammahat2 <- predict(res2)
  nuhat2 <- resid(res2)
  gammahat2 <- coef(res2)
  nuhat <- c(nuhat1,nuhat2)
  Wgammahat <- c(Wgammahat1,Wgammahat2)
  
  # outcome equation
  res1 <- lm(RY1~RX1+W1-1)
  rhohat1 <- coef(res1)[1]
  betahat1 <- coef(res1)[-1]
  epsilonhat1 <- resid(res1)
  res2 <- lm(RY2~RX2+W2-1)
  rhohat2 <- coef(res2)[1]
  betahat2 <- coef(res2)[-1]
  epsilonhat2 <- resid(res2)
  epsilonhat <- c(epsilonhat1,epsilonhat2)
  
  # construct H1
  H11 <- c(epsilonhat1 * nuhat1, rep(0,n/2))
  H12 <- c(rep(0,n/2), epsilonhat2 * nuhat2)
  
  # construct H2
  H21fn <- function(xy) mean((G==1) * (Ifn(xy[2],Y)-rhohat1*Ifn(xy[1],X)-c(W%*%betahat1)) * nuhat)
  H21 <- apply(cbind(X,Y), 1, H21fn)
  H22fn <- function(xy) mean((G==2) * (Ifn(xy[2],Y)-rhohat2*Ifn(xy[1],X)-c(W%*%betahat2)) * nuhat)
  H22 <- apply(cbind(X,Y), 1, H22fn)
  
  # construct H3
  H31fn <- function(x) mean((G==1) * epsilonhat * (Ifn(x,X)-Wgammahat))
  H31 <- sapply(X, H31fn)
  H32fn <- function(x) mean((G==2) * epsilonhat * (Ifn(x,X)-Wgammahat))
  H32 <- sapply(X, H32fn)
  
  # compute asymptotic variance
  sigma2hat <- c(mean((H11+H21+H31)^2) / mean((G==1)*nuhat^2)^2,
                 mean((H12+H22+H32)^2) / mean((G==2)*nuhat^2)^2)
  
  # save the result
  save(sigma2hat, Y, W, X, G, n, file = file.path("tests", "testthat", "testdata",
                                               paste0("grouped_lmranks_cov_sigmahat_covariates_", covariates, ".rda")))
}
