rrreg <- function(Y, X, W=NULL, omega=0, increasing=FALSE, na.rm=FALSE) {
	RY <- frank(Y, omega=omega, increasing=increasing, na.rm=na.rm)		
	RX <- frank(X, omega=omega, increasing=increasing, na.rm=na.rm)	 
	
	if (any(is.null(W))) {
		res <- lm(RY ~ RX)
	} else {
		res <- lm(RY ~ RX + W)
	}
	rhohat <- coef(summary(res))[2,1]
	se <- rrregSE(Y, X, W, omega=omega, increasing=increasing, na.rm=na.rm)
	return(list(rhohat=rhohat, se=se))
}

rrregSE <- function(Y, X, W=NULL, omega, increasing, na.rm) 
  return( sqrt(rrregAVar(Y, X, W, omega=omega, increasing=increasing, na.rm=na.rm) / length(Y)) )

#rrreg_se_private <- function()

rrregAVar <- function(Y, X, W=NULL, omega=omega, increasing=increasing, na.rm=na.rm) {
	RYhat <- function(y) mean(Ifn(y, Y, omega=omega, increasing=increasing, na.rm=na.rm))
	RY <- frank(Y, omega=omega, increasing=increasing, na.rm=na.rm)	
	RXhat <- function(x) mean(Ifn(x, X, omega=omega, increasing=increasing, na.rm=na.rm))
	RX <- frank(X, omega=omega, increasing=increasing, na.rm=na.rm)

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
	h2fn <- function(xy) mean((Ifn(xy[2],Y,omega=omega, increasing=increasing, na.rm=na.rm)-rhohat*Ifn(xy[1],X,omega=omega, increasing=increasing, na.rm=na.rm)-c(W%*%betahat)) * nuhat)
	h2 <- apply(cbind(X,Y), 1, h2fn)

	# construct h3
	h3fn <- function(x) mean(epsilonhat * (Ifn(x,X,omega=omega, increasing=increasing, na.rm=na.rm)-Wgammahat))
	h3 <- sapply(X, h3fn)

	# compute asymptotic variance
	sigmahat <- mean((h1+h2+h3)^2) / var(nuhat)^2

	return(sigmahat)
}

# u: single number; v - vector
# Returns a vector of length v with following values (for increasing = FALSE):
# 1 if v[i] > u
# 0 if v[i] < u
# omega if v[i] == u

Ifn <- function(u, v, omega=0, increasing=FALSE, na.rm=FALSE) {
	if (increasing) {
		return( omega*(u<=v) + (1-omega)*(u<v) )
	} else {
		return( omega*(u>=v) + (1-omega)*(u>v) )
	}
}

Ifn_vectorized <- function(x, omega, increasing, na.rm){
  # It is used only in self-context (u is always in vector v)
  # And it is used for every u in v
  # Return: matrix M. For increasing = FALSE M[i,j] = 
  # 1 if v[i] < v[j]
  # 0 if v[i] > v[j]
  # omega if v[i] == v[j]
  # TODO: move to utils and change irank so it uses the Ifn internally
  # Order matters!
  if(!increasing){
    x <- -x
  }
  if (na.rm){
    was_NA <- rep(FALSE, sum(is.na(x)))
  } else {
    was_NA <- is.na(x)
  }
  
  x <- x[!is.na(x)]
  n_NAs <- sum(was_NA)
  
  ranking <- order(x)
  sorted <- x[ranking]
  equal_to_next <- c(diff(sorted) == 0, FALSE)
  # block is a sequence of equal values
  # In the outcome matrix for sorted values, those are blocks on diagonal
  block_ends <- which(!equal_to_next)
  block_sizes <- diff(c(0, block_ends))
  block_starts <- block_ends - block_sizes + 1
  
  out_for_sorted <- matrix(0, nrow = length(x), ncol = length(x))
  out_for_sorted[lower.tri(out_for_sorted)] <- 1
  for(i in 1:length(block_ends)){
    out_for_sorted[block_starts[i]:block_ends[i],
                   block_starts[i]:block_ends[i]] <- omega
  }
  
  # return in order of original x
  original_order <- order(ranking)
  out_without_nas <- out_for_sorted[original_order,][,original_order]
  
  # correct for NAs
  out <- matrix(NA, nrow = length(x) + n_NAs, ncol = length(x) + n_NAs)
  out[!was_NA, !was_NA] <- out_without_nas
  out
}