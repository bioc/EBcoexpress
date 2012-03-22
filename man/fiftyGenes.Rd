\name{fiftyGenes}
\alias{fiftyGenes}
\docType{data}
\title{
The fiftyGenes expression matrix
}
\description{
A simulated expression matrix for fifty genes in two conditions, with
one hundred chips and twenty-five chips in the two coniditions.
Most gene pairs are uncorrelated, but all gene pairs involving only genes X1
to X25, or only X26 to X50, are DC. So for instance, X1~X2 is DC but X1~X30
is not
}
\usage{data(fiftyGenes)}
\format{
 num [1:50, 1:125] 0.655 0.0188 1.0786 1.6856 0.4814 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:50] "X1" "X2" "X3" "X4" ...
  ..$ : chr [1:125] "C1x1" "C1x2" "C1x3" "C1x4" ...
}
\examples{
data(fiftyGenes)
}
\keyword{datasets}
