
computeT <- function(X) return(max(computeTSj(X)))

computeTSj <- function(X) {
	X.sorted <- sort(X, decreasing=TRUE)
	p <- length(X)
	diff <- X.sorted[1:(p-1)] - X.sorted[2:p]
	return(diff)
}

clusterRanks <- function(X, Sigma, R=1000, alpha=0.05) {

	p <- length(X)
	stopifnot(dim(Sigma)==c(p,p))

	# compute T
	T <- computeT(X)

	# compute critical value
	Z <- mvrnorm(n=R, mu=rep(0,p), Sigma=Sigma)
	Tstar <- apply(Z, 1, computeT)
	cval <- quantile(Tstar, prob=1-alpha)

	# collect rejections
	TSj <- computeTSj(X)
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