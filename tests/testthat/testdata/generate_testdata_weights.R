library(csranks)
set.seed(100)

n <- 1000
X <- rnorm(n)
W <- matrix(rnorm(n * 2), n, 2)
Y <- X + rowSums(W) + rnorm(n, 1, 0.5)

weights <- runif(n, max = 2)
normalized_weights <- weights / mean(weights)

# compute ranks
RY <- frank(Y, increasing = TRUE, weights = weights, omega = 1)
RX <- frank(X, increasing = TRUE, weights = weights, omega = 1)

# ------- compute asymptotic variance "by hand"

Ifn <- function(u, v) {
  return(u <= v)
}

# first stage
res <- lm(RX ~ W, weights = weights)
Wgammahat <- predict(res)
nuhat <- resid(res)
gammahat <- coef(res)

# outcome equation
main_fit <- lm(RY ~ RX + W, weights = weights)
rhohat <- coef(main_fit)["RX"]
betahat <- coef(main_fit)[3:4]
epsilonhat <- resid(main_fit)
intercept <- coef(main_fit)[1]

# construct h1
# h1 <- epsilonhat * nuhat
h1 <- normalized_weights * epsilonhat * nuhat

# construct h2
h2fn <- function(xy) mean((Ifn(xy[2], Y) - rhohat * Ifn(xy[1], X) - c(W %*% betahat) - intercept) * nuhat * normalized_weights)
h2 <- apply(cbind(X, Y), 1, h2fn) * normalized_weights

# construct h3
h3fn <- function(x) mean(epsilonhat * (Ifn(x, X) - Wgammahat) * normalized_weights)
h3 <- sapply(X, h3fn) * normalized_weights

# compute asymptotic variance
sigma2hat <- mean((h1 + h2 + h3)^2) / mean(nuhat^2 * normalized_weights)^2

# so actually, the weights cancel each other out in sigma2hat?

# save the result
save(sigma2hat, Y, W, X, n, h1, h2, h3, weights,
  file = file.path("tests", "testthat", "testdata", "lmranks_cov_sigmahat_weighted.rda")
)
