\name{showPair}
\alias{showPair}
\title{
A function for looking at a pair's differential co-expression
}
\description{
This function plots the expression data for a given pair, coloring by
condition. A regression line may be added; this line may be made robust
so that it only uses those data points used by biweight midcorrelation
}
\usage{
showPair(pair, X, conditions, gsep = "~", regLine = TRUE, useBWMC = TRUE, colors = NULL, ...)
}
\arguments{
  \item{pair}{
A pair name, such as ABC~XYZ, where ~ is the separator given by gsep
}
  \item{X}{
An m-by-n expression matrix, where rows are genes and columns are chips
(subjects); include all chips in X, indicate condition in the
conditions array
}
  \item{conditions}{The conditions array
}
  \item{gsep}{
A separator that indicates a gene-pair, such as P53~MAPK1. The separator
should not appear in any of the gene names
}
  \item{regLine}{
Should a regression line be drawn for each condition?
}
  \item{useBWMC}{
Should the regression line be robust, a la biweight midcorrelation?
}
  \item{colors}{Colors for the different conditions. Defaults to palette()
}
  \item{\dots}{
Other options to be passed to plot(), with three exceptions.
The lty= and lwd= options will be passed to abline()
and will have an effect on the plot when regLine=TRUE;
col= is overwritten by the colors= array and may not be specified.
All other ... options will be passed to the main plot.
}
}
\value{Returns invisble(NULL)
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
showPair("X1~X2",fiftyGenes, conditions=tinyCond, pch=20,
  xlim=c(-4,4), ylim=c(-4,4))
#
showPair("X26~X35",fiftyGenes, conditions=tinyCond, pch=20,
  xlim=c(-4,4), ylim=c(-4,4))
#
showPair("X1~X35",fiftyGenes, conditions=tinyCond, pch=20,
  xlim=c(-4,4), ylim=c(-4,4))
}
\keyword{ hplot }