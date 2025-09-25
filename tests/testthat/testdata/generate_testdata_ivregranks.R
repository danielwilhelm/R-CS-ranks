library(csranks)
library(ivreg)
set.seed(100)

#######################
### vcov.ivregranks ###
#######################

for (covariates in c(TRUE, FALSE)) {
  n <- 10000
  if (covariates) {
    Z <- rnorm(n)
    W <- matrix(rnorm(n * 2), n, 2)
    X <- Z + rowSums(W) + rnorm(n, 0, 0.5)
    Y <- X + rowSums(W) + rnorm(n, 0, 1)
    W <- cbind(1, W)
  } else {
    Z <- rnorm(n)
    X <- Z + rnorm(n, 0, 0.5)
    Y <- X + rnorm(n, 0, 1)
    W <- matrix(1, n, 1)
  }

  RY <- frank(Y, increasing = TRUE)
  RX <- frank(X, increasing = TRUE)
  RZ <- frank(Z, increasing = TRUE)

  # ------- compute asymptotic variance "by hand"
  Ifn <- function(u, v) u <= v

  res1 <- lm(RZ ~ W - 1)
  Wgammahat <- predict(res1)
  nuhat <- resid(res1)
  gammahat <- coef(res1)

  res2 <- lm(RX ~ W - 1)
  zetahat <- resid(res2)

  res3 <- ivreg(RY ~ RX + W - 1 | RZ + W - 1)
  rhohat <- coef(res3)[1]
  betahat <- coef(res3)[-1]
  epsilonhat <- resid(res3)

  # construct h1
  h1 <- epsilonhat * nuhat

  # construct h2
  h2fn <- function(xy) {
    mean((Ifn(xy[2], Y) - rhohat * Ifn(xy[1], X) - c(W %*% betahat)) * nuhat)
  }
  h2 <- apply(cbind(X, Y), 1, h2fn)

  # construct h3
  h3fn <- function(z) mean(epsilonhat * (Ifn(z, Z) - Wgammahat))
  h3 <- sapply(Z, h3fn)

  # compute asymptotic variance
  sigma2hat <- mean((h1 + h2 + h3)^2) / mean(nuhat * zetahat)^2

  save(sigma2hat, Y, X, W, Z, n, h1, h2, h3, file = file.path(
    "tests", "testthat",
    "testdata", paste0(
      "ivregranks_cov_sigmahat_covariates_",
      covariates, ".rda"
    )
  ))
}

for (n in c(10, 50, 100)) {
  Z <- rnorm(n)
  W <- matrix(rnorm(n * 2), n, 2)
  X <- Z + rowSums(W) + rnorm(n, 0, 0.5)
  Y <- X + rowSums(W) + rnorm(n, 0, 1)
  W <- cbind(1, W)

  RY <- frank(Y, increasing = TRUE)
  RX <- frank(X, increasing = TRUE)
  RZ <- frank(Z, increasing = TRUE)

  Ifn <- function(u, v) u <= v

  res1 <- lm(RZ ~ W - 1)
  Wgammahat <- predict(res1)
  nuhat <- resid(res1)
  gammahat <- coef(res1)

  res2 <- lm(RX ~ W - 1)
  zetahat <- resid(res2)

  res3 <- ivreg(RY ~ RX + W - 1 | RZ + W - 1)
  rhohat <- coef(res3)[1]
  betahat <- coef(res3)[-1]
  epsilonhat <- resid(res3)

  # construct h1
  h1 <- epsilonhat * nuhat

  # construct h2
  h2fn <- function(xy) {
    mean((Ifn(xy[2], Y) - rhohat * Ifn(xy[1], X) - c(W %*% betahat)) * nuhat)
  }
  h2 <- apply(cbind(X, Y), 1, h2fn)

  # construct h3
  h3fn <- function(z) mean(epsilonhat * (Ifn(z, Z) - Wgammahat))
  h3 <- sapply(Z, h3fn)

  # compute asymptotic variance
  sigma2hat <- mean((h1 + h2 + h3)^2) / mean(nuhat * zetahat)^2

  save(sigma2hat, Y, X, W, Z, n, h1, h2, h3, file = file.path(
    "tests", "testthat",
    "testdata", paste0(
      "ivregranks_cov_sigmahat_n_",
      n, ".rda"
    )
  ))
}

########################
### increasing=FALSE ###
########################

n <- 50

Z <- rnorm(n)
W <- matrix(rnorm(n * 2), n, 2)
X <- Z + rowSums(W) + rnorm(n, 0, 0.5)
Y <- X + rowSums(W) + rnorm(n, 0, 1)
W <- cbind(1, W)

RY <- frank(Y, omega = 1, increasing = FALSE)
RX <- frank(X, omega = 1, increasing = FALSE)
RZ <- frank(Z, omega = 1, increasing = FALSE)

Ifn <- function(u, v) u >= v

res1 <- lm(RZ ~ W - 1)
Wgammahat <- predict(res1)
nuhat <- resid(res1)
gammahat <- coef(res1)

res2 <- lm(RX ~ W - 1)
zetahat <- resid(res2)

res3 <- ivreg(RY ~ RX + W - 1 | RZ + W - 1)
rhohat <- coef(res3)[1]
betahat <- coef(res3)[-1]
epsilonhat <- resid(res3)

# construct h1
h1 <- epsilonhat * nuhat

# construct h2
h2fn <- function(xy) {
  mean((Ifn(xy[2], Y) - rhohat * Ifn(xy[1], X) - c(W %*% betahat)) * nuhat)
}
h2 <- apply(cbind(X, Y), 1, h2fn)

# construct h3
h3fn <- function(z) mean(epsilonhat * (Ifn(z, Z) - Wgammahat))
h3 <- sapply(Z, h3fn)

# compute asymptotic variance
sigma2hat <- mean((h1 + h2 + h3)^2) / mean(nuhat * zetahat)^2

save(sigma2hat, Y, X, W, Z, n, h1, h2, h3, file = file.path(
  "tests", "testthat",
  "testdata", "ivregranks_cov_sigmahat_increasing_FALSE.rda"
))
