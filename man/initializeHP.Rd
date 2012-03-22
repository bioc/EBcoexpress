\name{initializeHP}

\alias{initializeHP}

\title{Initializing Hyperparameters for EM}

\description{
A function for initializing the EM hyperparameters. While the user
is free to do this in any manner s/he deems fit, we use
the excellent Mclust approach of package R/mclust
}

\usage{
initializeHP(D, conditions, seed = NULL, plottingOn = FALSE, colx = "red",
applyTransform = TRUE, verbose = 0, subsize = NULL, ...)
}

\arguments{
  \item{D}{
The correlation matrix output of makeMyD()
}
  \item{conditions}{
The conditions array
}
  \item{seed}{
A seed for making this procedure deterministic
}
  \item{plottingOn}{
Should the weighted vs. unweighted comparison plots be shown to the user?
Default is FALSE (no plotting)
}
  \item{colx}{
The color of the fitted empirical density curve, if plottingOn is TRUE
}
  \item{applyTransform}{
Should Fisher's Z-transformation be applied to the correlations?
}
  \item{verbose}{
An option to control comments as initialization proceeds. Set to 1 to
see comments, default is 0
}
  \item{subsize}{
A value less than the 1st dimension of D, if it is desired that only
some of the pairs be used in the estimation process (for computational reasons)
}
  \item{\dots}{
Other options to be passed to plot()
}
}

\details{
initializeHP() initializes the hyperparameters by asking Mclust to find the
1-, 2- or 3- component Normal mixture model that best fits the
(transformed) correlations as a whole. Mclust directly returns
estimates for G and the MUS; TAUS are estimated using Mclust's estimates
and sample sizes, per our model. WEIGHTS are estimated
using the mixture component classifications Mclust provides;
however, it is unclear whether those classifications should
we weighted by how confident Mclust is in their accuracy.
We have tried both approaches and neither is superior to
the other in all cases, so we compute both, compare the model
fits and return the WEIGHTS that best empirically fit the data.
This process is shown visually if plottingOn is TRUE
and the comments (if enabled) will describe the comparison.
In the event that the TAUS are estimated to be less than 0,
half of the (TAUS+SS_VARIANCE) estimate provided by Mclust
is used for TAUS; we especially do not recommend the use of
ebCoexpressZeroStep when this occurs (as opposed to generally
favoring the one-step version over the zero-step)
}

\value{
A list with four components, describing the hyperparameters:
 \item{G}{The number of mixture components (1, 2 or 3)}
 \item{MUS}{Means across mixture components}
 \item{TAUS}{Standard deviations across mixture components}
 \item{WEIGHTS}{Weights across mixture components (these sum to 1)}
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
data(fiftyGenes)
tinyCond <- c(rep(1,100),rep(2,25))
tinyPat <- ebPatterns(c("1,1","1,2"))
D <- makeMyD(fiftyGenes, tinyCond, useBWMC=TRUE)
set.seed(3)
initHP <- initializeHP(D, tinyCond)
}

\keyword{ smooth }
