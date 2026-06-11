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

n <- 1000
Z <- rnorm(n)
W <- rnorm(n)
X <- Z + W + rnorm(n, 0, 0.5)
Y <- X + W + rnorm(n, 0, 1)

RY <- frank(Y, increasing = TRUE)
RX <- frank(X, increasing = TRUE)
RZ <- frank(Z, increasing = TRUE)

Ifn <- function(u, v) u <= v

main_model <- ivreg(RY ~ RX + W | RZ + W)
RX_fitted_values <- RX - resid(main_model, component = "stage1")

stage_1_model <- lm(RX ~ RZ + W)

projection_model_1 <- ivreg(W ~ RX | RZ)
theta_1 <- coef(projection_model_1)["RX"]
v_hat_1 <- W - theta_1 * RX_fitted_values - coef(projection_model_1)["(Intercept)"]

rhohat <- coef(main_model)["RX"]
betahat_W <- coef(main_model)["W"]
intercept <- coef(main_model)["(Intercept)"]
epsilonhat <- resid(main_model)

# construct h1
h1 <- epsilonhat * v_hat_1

# construct h2
h2fn <- function(xy) {
  mean((Ifn(xy[2], Y) - rhohat * Ifn(xy[1], X) - W * betahat_W - intercept) * v_hat_1)
}
h2 <- apply(cbind(X, Y), 1, h2fn)

# construct h3
h3fn <- function(z) {
  prediction_data <- data.frame(W = W, RZ = as.numeric(Ifn(z, Z)))
  X_fitted_using_ifn <- predict(stage_1_model, prediction_data)
  new_v_hat <- (W - theta_1 * X_fitted_using_ifn - coef(projection_model_1)["(Intercept)"])
  mean(epsilonhat * new_v_hat)
}
h3 <- sapply(Z, h3fn)

# compute asymptotic variance
sigma2hat <- mean((h1 + h2 + h3)^2) / mean(v_hat_1^2)^2

save(sigma2hat, Y, X, W, Z, n, h1, h2, h3, file = file.path(
  "tests", "testthat",
  "testdata", "ivregranks_cov_sigmahat_regressor_1.rda"
))

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
