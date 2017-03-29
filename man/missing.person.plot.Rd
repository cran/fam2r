\name{missing.person.plot}
\alias{missing.person.plot}
\alias{internalID}
\alias{getSex}
\title{
Plots hypotheses for family reunion
}
\description{
A panel with two plots are generated. To the left, the one corresponding to
the person of interest (POI) being identical to the Missing Person (MP).
To the right these persons are urelated.
}
\usage{
missing.person.plot(ped_related, missing, id.labels=NULL,
  available="shaded", marker=NULL, width=c(4,4,1), newdev=TRUE,
	frametitles=c("H1: POI related", "H2:POI unrelated"), ...)
	
internalID(x, orig.ids)
getSex(x, orig.ids) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped_related}{
A \code{\link[paramlink]{linkdat}} object.
}
  \item{missing}{
An integer giving the MP.
}
  \item{id.labels}{
Character vector naming individuals in plot.
}
  \item{available}{
Character determining how genotyped individuals are displayed.
}
  \item{marker}{
Integer vector indicating markers to be displayed.
}
  \item{width}{
A numeric vector of relative widths of the subplots. Recycled to \code{length(plot.arg.list)} if necessary, before passed on to layout. Note that the vector does not need to sum to 1.
}
\item{newdev}{
logical. If \code{TRUE}, new plot window is generated.
}
  \item{frametitles}{
A character vector of length two giving titles.
}
  \item{\dots}{
Further arguments passed on to each call to \code{\link[paramlink]{plotPedList}}
}
\item{x}{
\code{\link[paramlink]{linkdat}} object.
}
\item{orig.ids}{
Character or integer vector,
}

}
\details{
See  example in the Examples section below.
}


\author{
Magnus Dehli Vigeland and Thore Egeland
}

\seealso{
\code{\link{plotPedList}}
}
\examples{
data(grandmother)
pedigrees = grandmother$pedigrees
datamatrix = grandmother$datamatrix
datamatrix[3,]=1
loci = grandmother$loci
x1 = Familias2linkdat(pedigrees, datamatrix, loci)
missing.person.plot(x1[[1]],2, available="shaded", frametitles=c("",""), 
  newdev=TRUE, width=c(1,2), marker=1)
}
