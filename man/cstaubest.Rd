% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranks.R
\name{cstaubest}
\alias{cstaubest}
\title{Projection confidence sets for the tau-best}
\usage{
cstaubest(x, sd, tau = 2, coverage = 0.95, stepdown = TRUE,
  R = 1000, na.rm = FALSE, seed = NA)
}
\arguments{
\item{x}{vector of estimates}

\item{sd}{vector of standard errors of \code{x}}

\item{tau}{the confidence set contains indicators for the elements in \code{x} whose rank is less than or equal to \code{tau}.}

\item{coverage}{nominal coverage of the confidence set. Default is 0.95.}

\item{stepdown}{logical; if \code{TRUE} (default), stepwise procedure is used, otherwise single step procedure is used.}

\item{R}{number of bootstrap replications. Default is 1000.}

\item{na.rm}{logical; if \code{TRUE}, then \code{NA}'s are removed from \code{x} and \code{sd} (if any).}

\item{seed}{seed for bootstrap random variable draws. If set to \code{NA} (default), then seed is not set.}
}
\value{
logical vector indicating which of the elements of \code{x} are in the confidence set for the tau-best.
}
\description{
Projection confidence sets for the tau-best
}
\section{Details}{

Implentation of the projection confidence sets for the tau-best proposed in Mogstad, Romano, Shaikh, and Wilhelm ("Inference for Ranks with Applications to Mobility across Neighborhoods and Academic Achievements across Countries", \href{https://www.ucl.ac.uk/~uctpdwi/papers/cwp1020.pdf}{CeMMAP Working Paper CWP10/20}). The confidence set contains indicators for the elements in \code{x} whose rank is less than or equal to \code{tau} with probability approximately equal to the coverage indicated in \code{coverage}. 
Parametric bootstrap based on the normal distribution with independent populations.
}

\examples{
x <- seq(1,3,length=10)
sd <- rep(0.2,10)
cstaubest(x, sd, tau=3)
}
