#' Linear model for ranks
#' 
#' Fit a linear model with a single rank response and
#' a single rank covariate (and possibly other usual covariates).
#'
#' @param Y Numeric vector. Based on it, ranks for response will be constructed
#' and used as final response in the linear model.
#' @param X Numeric vector. Based on it, ranks for this covariate will be constructed
#' and used as final covariate in the linear model.
#' @param W Matrix. Other covariates to include in the model. Should not include an intercept column.
#' @param omega numeric; numeric value in [0,1], each corresponding to a different definition of the rank; default is \code{0}. See \code{\link{frank}} for details.
#' @param increasing logical; if \code{TRUE}, then large elements in \code{X} and \code{Y} receive a large rank. Otherwise, large elements receive small ranks. 
#' @param na.rm logical; if \code{TRUE}, then \code{NA}'s are removed from \code{X} and \code{Y} (if any). 

#' @return A list with two items: 
#' \itemize{
#' \item{\code{rhohat} - the linear coefficient for rank of X.}
#' \item{\code{se} - Estimated standard error of the coefficient rhohat.}
#' }
#' 
#' @examples
#' Y <- c(3,1,2,4,5)
#' y_frank <- c(0.6, 1.0, 0.8, 0.4, 0.2)
#' X <- 1:5
#' omega <- 0.5
#' x_frank <- c(1.0, 0.8, 0.6, 0.4, 0.2)
#' W <- matrix(y_frank * 0.1 + 5 + rnorm(5, sd = 0.1), ncol = 1)
#'
#' simple_lmranks(Y, X, W)
#' 
#' @importFrom stats lm predict resid coef var
#' @export
simple_lmranks <- function(Y, X, W=NULL, omega=0, increasing=FALSE, na.rm=FALSE) {
  I_Y <- compare(Y, omega=omega, increasing=increasing, na.rm=na.rm)
  I_X <- compare(X, omega=omega, increasing=increasing, na.rm=na.rm)
  RY <- (rowSums(I_Y) + 1 - omega) / length(Y)
  RX <- (rowSums(I_X) + 1 - omega) / length(X)
  
	if (any(is.null(W))) {
		W <- matrix(rep(1,length(Y)), ncol = 1)
	} else {
	  W <- cbind(rep(1,length(Y)),W)
	}
  
  RX_by_Ws <- lm(RX ~ W-1)
  Wgammahat <- predict(RX_by_Ws)
  nuhat <- resid(RX_by_Ws)
  gammahat <- coef(RX_by_Ws)
  
  res <- lm(RY~RX+W-1)
  rhohat <- coef(res)[1]
  betahat <- coef(res)[-1]
  epsilonhat <- resid(res)
  
  h1 <- epsilonhat * nuhat
  
  # vectorized `-` and `*` goes over rows
  h2 <- colMeans((t(I_Y - rhohat * I_X) - c(W %*% betahat)) * nuhat)
  
  h3 <- colMeans(epsilonhat * (t(I_X) - Wgammahat))
  
  sigmahat <- mean((h1+h2+h3)^2) / var(nuhat)^2
  se <- sqrt(sigmahat / length(Y))
  names(rhohat) <- NULL
	return(list(rhohat=rhohat, se=se))
}

#' @describeIn simple_lmranks Calculate the standard error of rank regression coefficient

simple_lmranks_rho_se <- function(Y, X, W=NULL, omega=0, increasing=FALSE, na.rm=FALSE){
  lmranks_outcome <- simple_lmranks(Y, X, W=W, omega=omega, increasing=increasing, na.rm=na.rm)
  return(lmranks_outcome$se)
}