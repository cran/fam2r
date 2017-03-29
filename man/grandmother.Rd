\name{grandmother}
\alias{grandmother}
\alias{Demo3Markers}
\alias{symmetric}
\alias{F21}
\alias{E004}
\alias{dc}
\alias{adoption}
\docType{data}
\title{
Dataset for missing grandchildren paper
}
\description{
The datasets are used for examples in Kling et al. (2017) or to validate
conditional simulation of marker data followed by likelihood calculation for alternative pedigrees.
Further brief description is given for the examples below,
see \code{\link{FamiliasConditional}} or \code{\link{paramlinkConditional}} for examples involving computations.
}
\usage{
data(grandmother)
data(Demo3Markers)
data(symmetric)
data(F21)
data(dc)
data(adoption)
}
\format{
A list with components
pedigrees, datamatrix, loci, see \code{\link{FamiliasLocus}}
}

\examples{
data(grandmother) #One marker
data(Demo3Markers) #Three markers, one with mutation
#Three pedigrees, Half-sibs, avuncular and grandparent:
data(symmetric) 
data(E004) # E zero zero four
pedigrees = E004$pedigrees
datamatrix = E004$datamatrix
loci = E004$loci
avail = 7
persons = rownames(datamatrix)
x = Familias2linkdat(pedigrees, datamatrix, loci) 
missing.person.plot(x[[2]][[1]], missing=avail, marker=1:2, newdev=TRUE,
                    frametitle=c("H1", "H2"), dev.height=3.3, dev.width=9)
x=Familias2linkdat(pedigrees, datamatrix, loci)
missing.person.plot(x[[2]][[1]], missing=avail, marker=1:2, newdev=TRUE,
                    frametitle=c("H1", "H2"))
\dontrun{ #Takes a few minutes
Nsim=1000
res.new = conditionalLR(Nsim=Nsim, datamatrix, loci, pedigrees, available=avail,
                     seed=1483118516, program = "Familias", ref=1 )
LRnew.mean =cbind(res.new[,,1][,2],res.new[,,2][,2])
apply(LRnew.mean,2,quantile)
}

}

