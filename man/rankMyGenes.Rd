\name{rankMyGenes}
\alias{rankMyGenes}
\title{A function to rank the genes by the number of DC pairs
in which they appear
}
\description{
This function uses a threshold to determine the names of the DC pairs.
It then splits those pairs into their constituent genes and tables them.
A sorted version of that table is then returned. This information
may be useful for those investigating `differential hubbing' -- see
the Hudson et al. reference for more information
}
\usage{
rankMyGenes(emOut, thresh = 0.95, sep = "~")
}
\arguments{
  \item{emOut}{The output of an ebCoexpressSeries function call
}
  \item{thresh}{A threshold for determining whether a pair is DC. This
may be set as a hard threshold (default is hard 5% FDR) or a soft
threshold, as returned by crit.fun
}
  \item{sep}{The separator used in the pair names
}
}
\value{A sorted, named array of gene counts}
\references{
Dawson JA and Kendziorski C. An empirical Bayesian approach for
identifying differential co-expression in high-throughput experiments.
(2011) Biometrics. E-publication before print:
http://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2011.01688.x/abstract

Hudson NJ, Reverter A, Dalrymple BP (2009) A Differential Wiring Analysis
of Expression Data Correctly Identifies the Gene Containing the Causal
Mutation. PLoS Comput Biol 5(5): e1000382. doi:10.1371/journal.pcbi.1000382
}
\author{
John A. Dawson <jadawson@wisc.edu>
}
\seealso{ebCoexpressSeries, crit.fun
}

\examples{
data(fiftyGenes)
tinyCond <- c(rep(1,100),rep(2,25))
tinyPat <- ebPatterns(c("1,1","1,2"))
D <- makeMyD(fiftyGenes, tinyCond, useBWMC=TRUE)
set.seed(3)
initHP <- initializeHP(D, tinyCond)

zout <- ebCoexpressZeroStep(D, tinyCond, tinyPat, initHP)
rankMyGenes(zout)
}

\keyword{ univar }