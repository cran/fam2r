\name{PE}
\alias{PE}

\title{
Exclusion probabilities
}
\description{
This is a wrapper for \code{\link[paramlink]{exclusionPower}}.
The probability of exclusion (PE) is calculated for each marker and combined
for all markers.
}
\usage{
PE(pedigrees, datamatrix, loci, claim = 1, true = 2, available = NULL, 
  file = NULL, ignore = FALSE)
}

\arguments{
  \item{pedigrees}{
A \code{\link[Familias]{FamiliasPedigree}} object or a list of such.
}
  \item{datamatrix}{
A data frame, each row gives the genotype of one individual,
a pair of columns represents a marker.
}
  \item{loci}{
A \code{\link[Familias]{FamiliasLocus}} object or a list of such.
}
  \item{claim}{
Integer. Index of the claimed pedigree, typically the one where \code{MP=POI}.
}
  \item{true}{
Integer. Index of the assumed true pedigree for PE calculation, 
typically the one where POI is assumed to be an unrelated individual.
}
  \item{available}{
Character or integer identifying person to be simulated.
}
  \item{file}{
Charachter. If not \code{NULL}, a file is written.
}
\item{ignore}{
Logical. If \code{TRUE}, mutations are ignored.
}
}
\details{
Wrapper for \code{\link[paramlink]{exclusionPower}}
}
\value{
A data frame. The first column gives the name of the markers,
the second the PE. The last line gives the combined result for all markers.
}
\references{
Egeland, Pinto and Vigeland (2014).
}
\author{
Thore Egeland <Thore.Egeland@gmail.com>
}


\seealso{
 \code{\link[paramlink]{exclusionPower}}
}
\examples{
data(F21)
pedigrees = F21$pedigrees
datamatrix = F21$datamatrix
loci = F21$loci
available = "Missing Person"
PE(pedigrees, datamatrix, loci, claim = 1, true = 2, 
  available = available, file = NULL)
}

