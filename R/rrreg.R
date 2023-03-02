rrreg <- function(Y, X, W=NULL, ...) {
	RY <- frank(Y, ...)		
	RX <- frank(X, ...)	 
	
	if (any(is.null(W))) {
		res <- lm(RY ~ RX)
	} else {
		res <- lm(RY ~ RX + W)
	}
	rhohat <- coef(summary(res))[2,1]
	se <- rrregSE(Y, X, W, ...)
	return(list(rhohat=rhohat, se=se))
}

rrregSE <- function(Y, X, W=NULL, ...) return( sqrt(rrregAVar(Y, X, W, ...) / length(Y)) )

rrregAVar <- function(Y, X, W=NULL, ...) {
	RYhat <- function(y) mean(Ifn(y, Y, ...))
	RY <- frank(Y, ...)	
	RXhat <- function(x) mean(Ifn(x, X, ...))
	RX <- frank(X, ...)

	W <- cbind(rep(1,length(Y)),W)

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
	h2fn <- function(xy) mean((Ifn(xy[2],Y,...)-rhohat*Ifn(xy[1],X,...)-c(W%*%betahat)) * nuhat)
	h2 <- apply(cbind(X,Y), 1, h2fn)

	# construct h3
	h3fn <- function(x) mean(epsilonhat * (Ifn(x,X,...)-Wgammahat))
	h3 <- sapply(X, h3fn)

	# compute asymptotic variance
	sigmahat <- mean((h1+h2+h3)^2) / var(nuhat)^2

	return(sigmahat)
}

Ifn <- function(u, v, omega=0, increasing=FALSE, na.rm=FALSE) {
	if (increasing) {
		return( omega*(u<=v) + (1-omega)*(u<v) )
	} else {
		return( omega*(u>=v) + (1-omega)*(u>v) )
	}
}