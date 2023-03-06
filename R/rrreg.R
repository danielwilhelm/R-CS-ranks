rrreg <- function(Y, X, W=NULL, omega=0, increasing=FALSE, na.rm=FALSE) {
  I_Y <- compare(Y, omega=omega, increasing=increasing, na.rm=na.rm)
  I_X <- compare(X, omega=omega, increasing=increasing, na.rm=na.rm)
  RY <- rowMeans(I_Y)
  RX <- rowMeans(I_X)
  
	if (any(is.null(W))) {
		W <- rep(1,length(Y))
	} else {
	  W <- cbind(rep(1,length(Y)),W)
	}
  
  RX_by_W <- lm(RX ~ W-1)
  Wgammahat <- predict(RX_by_Ws)
  nuhat <- resid(RX_by_Ws)
  gammahat <- coef(RX_by_Ws)
  
  res <- lm(RY~RX+W-1)
  rhohat <- coef(res)[1]
  betahat <- coef(res)[-1]
  epsilonhat <- resid(res)
  
  h1 <- epsilonhat * nuhat
  
  h2 <- colMeans((I_Y - rhohat * I_X - c(W %*% betahat)) * nuhat)
  
  h3 <- colMeans(epsilonhat * (I_X - Wgammahat))
  
  sigmahat <- mean((h1+h2+h3)^2) / var(nuhat)^2
  se <- sqrt(sigmahat / length(Y))
	return(list(rhohat=rhohat, se=se))
}

rrregSE <- function(Y, X, W=NULL, omega, increasing, na.rm){
  rrreg_outcome <- rrreg(Y, X, W=W, omega=omega, increasing=increasing, na.rm=na.rm)
  return(rrreg_outcome$se)
}