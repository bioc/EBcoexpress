\name{showNetwork}
\alias{showNetwork}
\title{
A function for looking at the co-expression among a small group of genes
}
\description{
This function draws a network for a selected group of genes using igraph.
The edges are colored in accordance with the correlation strength
indicated by the inputted D matrix, ranging from red (strong negative
correlation) to blue (strong positive correlation)
}
\usage{
showNetwork(geneSet, D, condFocus, gsep = "~", layout = "kamada.kawai", seed = NULL, hidingThreshold=NULL, ...)
}
\arguments{
  \item{geneSet}{
An array of genes of interest; should not be larger than a dozen or so
}
  \item{D}{
The correlation matrix output of makeMyD()
}
  \item{condFocus}{
The condition of interest for this network. Should be one of the
integers in the conditions array
}
  \item{gsep}{
A separator that indicates a gene-pair, such as P53~MAPK1. The separator
should not appear in any of the gene names
}
  \item{layout}{
A layout to be parsed and used by igraph. Examples include circle (the default)
and kamada.kawai; see the documentation for igraph for more information.
At this time it is not possible to specify parameters specific to
particular layouts
}
  \item{seed}{
A seed to be set before invoking igraph's layout generation. This is useful
for layouts such as random, where node postion is not deterministic
}
  \item{hidingThreshold}{
A threshold which we will shorthand by 'h'. If this value is non-NULL, all
correlations in [-h, h] will not be plotted in the network. This is
useful for removing clutter in busy networks will relatively high (say, 20+)
numbers of genes
}
  \item{\dots}{
Other options to be passed to plot.igraph(). Networks generated
by igraph require quite a bit of formatting, and it is up to the
user to do so by specifying appropriate options from the following:

vertex.shape=, vertex.label.cex=, vertex.color=, vertex.frame.color=,
vertex.size=, vertex.label.color=, vertex.label.family=, and edge.width=

The following options are hard-coded and may not be overwritten:

vertex.label=geneSet, edge.arrow.mode=0, edge.color=[red/blue colors]

where [red/blue colors] is determined by the correlation information
contained in D, possibly overwritten in some cases if hidingThreshold
is non-NULL
}
}
\value{Returns invisible(NULL)
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

\seealso{igraph, igraph.layout
}

\examples{
data(fiftyGenes)
tinyCond <- c(rep(1,100),rep(2,25))
tinyPat <- ebPatterns(c("1,1","1,2"))
D <- makeMyD(fiftyGenes, tinyCond, useBWMC=TRUE)
twentyGeneNames <- dimnames(fiftyGenes)[[1]][c(1:10,26:35)]

showNetwork(twentyGeneNames, D, condFocus = 1, gsep = "~",
  layout = "kamada.kawai", seed = 5, vertex.shape="circle",
  vertex.label.cex=1, vertex.color="white", edge.width=2,
  vertex.frame.color="black", vertex.size=20,
  vertex.label.color="black", vertex.label.family="sans",
  hidingThreshold=0.3)
#
showNetwork(twentyGeneNames, D, condFocus = 2, gsep = "~",
  layout = "kamada.kawai", seed = 5, vertex.shape="circle",
  vertex.label.cex=1, vertex.color="white", edge.width=2,
  vertex.frame.color="black", vertex.size=20,
  vertex.label.color="black", vertex.label.family="sans",
  hidingThreshold=0.3)
#
}
\keyword{ hplot }
