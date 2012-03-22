\name{ebCoexpressSeries}

\alias{ebCoexpressFullTCAECM}
\alias{ebCoexpressOneStep}
\alias{ebCoexpressZeroStep}

\title{Functions to run single-study analyses}

\description{Core functions that run single-study DC analyses in EBcoexpress}

\usage{
ebCoexpressFullTCAECM(D, conditions, pattern, hpEsts, controlOptions = list())
ebCoexpressOneStep(D, conditions, pattern, hpEsts, controlOptions = list())
ebCoexpressZeroStep(D, conditions, pattern, hpEsts, controlOptions = list())
}

\arguments{
  \item{D}{
The correlation matrix output of makeMyD()
}
  \item{conditions}{
The conditions array
}
  \item{pattern}{
An appropriate output from ebPatterns()
}
  \item{hpEsts}{
The initial hyperparameter estimates from initializeHP() or some other method
}
  \item{controlOptions}{

A list with many options for controlling execution:

** These two are common to all members of the series:

 applyTransform: Should Fisher's Z-transformation be applied? Defaults to TRUE

 verbose: Controls auto-commenting; set to 0 to turn off comments

** These are used only by the OneStep and FullTCAECM versions:

 convtol: Convergence tolerance for the EM; default is 5e-04

 enforceFloor: Should EC proportion never drop below 0.8? Default is TRUE

 subsize: If non-NULL, a value less than the 1st dimension of D (p).
  The EM will use subsize randomly chosen pairs in its computations
  instead of all p pairs; by default, all pairs are used. We suggest
  use of this option when the number of pairs is very large

 m2MaxIter: Upper limit of the number of itertions the optimizer
  from minqa will use in the M2 step; defaults to 100. For unrestricted
  iterations, set this to NULL
}
}
\details{
These three functions represent different flavors of the TCAECM.
The FullTCAECM version will run a full TCAECM. The OneStep version
will perform a single iteration of the TCAECM and return the results.
The ZeroStep version does not perform any EM calculations and instead
uses the initial estimates of the hyperparameters to generate
posterior probabilities of DC.

We recommend using the OneStep version in most cases.
}
\value{

The output is a list with two members, MODEL and POSTPROBS:

MODEL is a list containing an array MIX and a list HPS.
MIX contains estimated mixing proportions for EC/DC classes.
HPS contains model specifications:
 \item{G}{The number of mixture components (1, 2 or 3)}
 \item{MUS}{Means across mixture components}
 \item{TAUS}{Standard deviations across mixture components}
 \item{WEIGHTS}{Weights across mixture components (these sum to 1)}

The list of lists required by ebCoexpressMeta() is obtained
by listing the separate analyses' HPS lists

POSTPROBS is a p-by-L matrix containing posterior probabilities
of EC and DC over all L EC/DC classes. The EC posterior probabilities
will always be in the first column (which should be fed into crit.fun()
if using the soft threshold). Total posterior probabilities of DC
for each gene pair are found by summing over the other L-1 columns
(or taking 1 minus the first (EC) column)

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

zout <- ebCoexpressZeroStep(D, tinyCond, tinyPat, initHP)
\dontrun{oout <- ebCoexpressOneStep(D, tinyCond, tinyPat, initHP)}
\dontrun{fout <- ebCoexpressFullTCAECM(D, tinyCond, tinyPat, initHP)}

softThresh <- crit.fun(zout$POSTPROB[,1], 0.05)
}

\keyword{ models }