\name{utilities}

\alias{crit.fun}
\alias{bwmc}

\title{Basic utilities for the EBcoexpress package}

\description{
At present there are two utilties: crit.fun() and bwmc(). The former is used
to compute soft thresholds for FDR control, the latter is like cor()
but uses biwieght midcorrelation instead of the usual Pearson's
correlation coefficient.
}

\usage{
crit.fun(ecPostProbs, targetFDR)
bwmc(X)
}

\arguments{
  \item{ecPostProbs}{
An array of posterior probabilities of equivalent coexpression for all pairs
}
  \item{targetFDR}{
A target FDR rate
}
  \item{X}{
An expression matrix in one condition where the rows correspond to genes
}
}

\details{
crit.fun() returns a soft threshold for FDR control. It is similar
to the function of the same name in the package EBarrays.
bwmc() computes the biweight midcorrelation for an expression matrix;
it is used internally to generate the D correlations matrix by makeMyD()
when useBWMC is TRUE. It is also a handy little function so we made it
visible at the top level. The guts of this function are in C for speed
}

\value{
crit.fun returns a single value; under a soft thresholding approach,
any pair with total posterior probability of differential co-expression
(i.e., 1 - posterior probability of equivalent co-expression) greater
than this value is deemed to be DC

If X has 1st dimension m, bwmc(t(X)) returns an m-by-m matrix of
pairwise biweight midcorrelations as a matrix, in a manner similar to cor().
}

\references{
Dawson JA and Kendziorski C. An empirical Bayesian approach for
identifying differential co-expression in high-throughput experiments.
(2011) Biometrics. E-publication before print:
http://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2011.01688.x/abstract
}

\author{
John A. Dawson <jadawson@wisc.edu>
}

\examples{
set.seed(1)
ecs <- c(runif(950),runif(50,0,0.01))
thresh <- crit.fun(ecs, 0.05)

set.seed(1)
X <- matrix(runif(10*100),10,100)
print(cor(t(X)))
print(bwmc(t(X)))
}

\keyword{ utilities }
\keyword{ manip }
