\name{makeMyD}
\alias{makeMyD}
\title{A function to convert the X expression matrix into the D correlation matrix}

\description{
A function to convert the X expression matrix into the D correlation matrix;
uses either Pearson's correlation coefficient or biweight midcorrelation
}

\usage{
makeMyD(X, conditions, useBWMC = FALSE, gpsep = "~")
}

\arguments{
  \item{X}{
An m-by-n expression matrix, where rows are genes and columns are chips
(subjects); include all chips in X, indicate condition in the
conditions array
%%     ~~Describe \code{X} here~~
}
  \item{conditions}{
The conditions array
}
  \item{useBWMC}{
Should biweight midcorrelation be used instead of Pearson's correlation
coefficient?
}
  \item{gpsep}{
A separator that indicates a gene-pair, such as P53~MAPK1. The separator
should not appear in any of the gene names
}
}

\value{
A p-by-K matrix of observed correlations for all p gene-pairs,
where p is choose(m,2), m is the 1st dimension of X and K is the
number of conditions specified by the conditions array
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
}

\keyword{ models }
