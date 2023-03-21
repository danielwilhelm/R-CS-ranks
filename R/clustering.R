library(MASS)

computeT <- function(X, Sigma=NA, studentize=TRUE) return(max(computeTSj(X, Sigma=Sigma, studentize=studentize)))

computeTSj <- function(X, Sigma=NA, studentize=TRUE) {
	stopifnot( !(any(is.na(Sigma)) & studentize) )

	X.sorted <- sort(X, decreasing=TRUE)
	p <- length(X)	
	
	if (studentize) {
		fn <- function(j) {
			sigmaij <- sqrt( outer(diag(Sigma)[1:j], diag(Sigma)[(j+1):p], '+') - 2*Sigma[1:j,(j+1):p] )
			return( min(outer(X.sorted[1:j], X.sorted[(j+1):p], '-') / sigmaij) )
		}
		TSj <- sapply(1:(p-1), fn)
	} else {
		TSj <- X.sorted[1:(p-1)] - X.sorted[2:p]
	}
	return(TSj)	
}

clusterRanks <- function(X, Sigma, R=1000, alpha=0.05, studentize=TRUE) {

	p <- length(X)
	stopifnot(dim(Sigma)==c(p,p))

	# compute T
	T <- computeT(X, Sigma=Sigma, studentize=studentize)

	# compute critical value
	Z <- mvrnorm(n=R, mu=rep(0,p), Sigma=Sigma)
	Tstar <- apply(Z, 1, computeT, Sigma=Sigma, studentize=studentize)
	cval <- quantile(Tstar, prob=1-alpha)

	# collect rejections
	TSj <- computeTSj(X, Sigma=Sigma, studentize=studentize)
	rej <- TSj > cval

	# create clusters
	if (any(rej)) {
		j1 <- 1
		SR <- (1:(p-1))[rej]
		if (max(SR)!=p) SR <- c(SR,p)
		clusters <- rep(NA, p)
		ind <- order(X, decreasing=TRUE)
		for (i in 1:length(SR)) {
			clusters[ind[j1:SR[i]]] <- i
			j1 <- SR[i]+1
		}
		return(clusters)
	} else {
		return(rep(1,p))
	}
}