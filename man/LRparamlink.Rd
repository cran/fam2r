\name{LRparamlink}
\alias{LRparamlink}
\title{
Calculates likelihoods and likelihood ratios using 'paramlink'
}
\description{
Typically export from  windows 'Familias' is converted to 
\code{\link[paramlink]{linkdat}}
using \code{\link[paramlink]{Familias2linkdat}} prior to calling the function.
}
\usage{
LRparamlink(x, ref, markers)
}

\arguments{
  \item{x}{
A list of pedigrees, in the form of  \code{\link{linkdat}} objects.
}
  \item{ref}{
Integer. Index of the pedigree to be used in the denominator of LR.
}

  \item{markers}{
Integer vector. Index for markers for which calculations should be done.
}
}

\value{
\item{LR }{Likelihood ratios.}
\item{LRperMarker }{Likelihood ratios for each marker.}
\item{likelihoodsPerSystem }{Likelihoods for each marker.}
\item{time }{User, system and elapsed time.}
}

\author{
Magnus Dehli Vigeland and Thore Egeland
}



\examples{
data(adoption)
x = Familias2linkdat(adoption$pedigrees, adoption$datamatrix, adoption$loci)
result = LRparamlink(x, ref=2)
# Only marker 11 and 33
result33 = LRparamlink(x, ref=2, marker=c(11,33))
}
