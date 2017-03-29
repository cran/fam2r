\name{fam2r-package}
\alias{fam2r-package}
\alias{fam2r}
\docType{package}
\title{
\packageTitle{fam2r}
}
\description{
Conditional simulation of marker data on pedigrees, plotting and likelihood calculations, 
using \code{\link{Familias}} or \code{\link{paramlink}}, on pedigrees. The package is 
specifically tailored to work on export of data from windows \code{\link{Familias}}.
}
\details{
The main function is \code{\link{conditionalLR}}. Marker data for a specified individual is simulated conditionally on 
pedigrees and previously typed individuals using \code{\link{markerSim}}. Then likelihoods are calculated using
\code{\link{Familias}} or \code{\link{paramlink}}. 
}
\author{
Thore Egeland <thore.egeland@gmail.com> and Magnus Dehli Vigeland.

Maintainer: Thore Egeland 
}
\references{
Kling et al (2017)
}

\examples{
# Example.There is one marker. A grand mother (GM=1) is
# is genotyped and we simulate the genotype of the grand son (GS=POI).
data(grandmother)
pedigrees = grandmother$pedigrees
datamatrix = grandmother$datamatrix
loci = grandmother$loci
x = Familias2linkdat(pedigrees, datamatrix, loci)
missing.person.plot(x[[1]], missing=3, marker = 1, newdev=TRUE)
Nsim = 5 #Increase to 1000
# Calculations can be done using 'Familias' or 'paramlink':
res1 = FamiliasConditional(Nsim = Nsim, datamatrix, loci, pedigrees, 
  truePed = 1, available = "GS", ref=2, seed=17)
res2 = paramlinkConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
  truePed = 1, available = "GS", ref=2, seed=17)
LR1 = data.frame(LR.H1.Familias=res1[[1]][,1], LR.H1.paramlink=res2[[1]][,1])
# Checks that 'Familias' and 'paramlink' gives the same result:
stopifnot(all(apply(LR1,1, function(x) abs(x[1]-x[2])<1e-10)))

# Altenatively, there is a wrapper function
res3 = conditionalLR(Nsim = Nsim, datamatrix, loci, pedigrees, 
  available = "GS", ref=2, seed=17, verbose=FALSE, 
  program="Familias", simplify=TRUE)
	 

}
