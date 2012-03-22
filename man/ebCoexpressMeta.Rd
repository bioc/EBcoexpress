\name{ebCoexpressMeta}

\alias{ebCoexpressMeta}

\title{A function for DC meta-analysis that combines individual study analyses}

\description{
This function performs a DC meta-analysis, using the hyperparameter estimates
obtained from individual-study DC analyses via one of the ebCoexpress series.
}

\usage{
ebCoexpressMeta(DList, conditionsList, pattern, hpEstsList, controlOptions = list())
}

\arguments{
  \item{DList}{
A list of the individual D matrices from the studies under consideration.
They should have the same dimensions, corresponding to the same gene-pairs
and the same conditions
}
  \item{conditionsList}{
A list of the individual conditions arrays from the studies under consideration
}
  \item{pattern}{
The appropriate output from ebPatterns()
}
  \item{hpEstsList}{
A list of the individual hyperparameter estimate objects from the studies
under consideration, as outputted from the ebCoexpress series
(or initializeHP(), the format is the same)
}
  \item{controlOptions}{

A list with many options for controlling execution:

 applyTransform: Should Fisher's Z-transformation be applied? Defaults to TRUE

 verbose: Controls auto-commenting; set to 0 to turn off comments

 convtol: Convergence tolerance for the EM; default is 5e-04

 enforceFloor: Should EC proportion never drop below 0.8? Default is TRUE
}
}

\details{
Since the meta-analysis model assumes that each study has its own
study-specific parameters, those parameters should be estimated
using a single-study DC function (one of the other members of the
ebCoexpress series); their outputs are used by the hpEstsList option.
The EM is then run to determine the system-wide mixing proportions,
which are used to compute meta posterior probabilities for all
EC/DC classes
}

\value{
The output is a list with two members, MODEL and POSTPROBS:

MODEL is a list containing an array MIX and a list HPS.
MIX contains estimated mixing proportions for EC/DC classes.
HPS is the inputed list of lists

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
D <- makeMyD(fiftyGenes, tinyCond, useBWMC=TRUE)
set.seed(3)
initHP <- initializeHP(D, tinyCond)

D1 <- D
D2 <- D
DList <- list(D1, D2)
cond1 <- tinyCond
cond2 <- tinyCond
conditionsList <- list(cond1, cond2)
pattern <- ebPatterns(c("1,1","1,2"))
initHP1 <- initHP
initHP2 <- initHP
out1 <- ebCoexpressZeroStep(D1, cond1, pattern, initHP1)
out2 <- ebCoexpressZeroStep(D2, cond2, pattern, initHP2)
hpEstsList <- list(out1$MODEL$HPS, out2$MODEL$HPS)

metaResults <- ebCoexpressMeta(
    DList, conditionsList, pattern, hpEstsList)
}

\keyword{ models }
