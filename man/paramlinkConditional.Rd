\name{paramlinkConditional}
\alias{paramlinkConditional}
\alias{paramlinkConditionalOne}
\title{
Conditional simulation of marker data on pedigrees and \code{\link[paramlink]{paramlink}} likelihood
}
\description{
Marker data is simulated for a specified person and several markers using \code{paramlinkConditionalOne} for each marker.
The only difference between \code{\link{FamiliasConditional}} and \code{paramlinkConditional} is that the former uses
the C implementation of \code{\link[Familias]{FamiliasPosterior}} for likelihood calculation while the latter is based on the \code{R} 
implementation in \code{\link{paramlink}}.
}
\usage{
paramlinkConditional(Nsim = 5, datamatrix, loci, pedigrees, 
	truePed = 1, ref=NULL, available = NULL, prior=NULL, seed = NULL)
paramlinkConditionalOne(Nsim = 5, mark = 1, ref = 2, datamatrix, loci, 
    pedigrees, truePed = 1, available = NULL, prior = NULL, seed = NULL)
}


\arguments{
  \item{Nsim}{
Integer. Number of simulations.
}
  \item{ref}{
Integer Denominator of LR.
}
  \item{datamatrix}{
A data frame or a matrix. The row names must be the names of the persons you have data for. The columns contain the alleles, two columns for each marker, in the same order used in the loci list.
}

  \item{loci}{
A list of \code{\link[Familias]{FamiliasLocus}} objects. 
}
  \item{pedigrees}{
An list with elements of type \code{\link[Familias]{FamiliasPedigree}}.
}
  \item{truePed}{
Integer. Index of pedigree from which marker data are simulated.
}
  \item{available}{
A character giving the name of the person to be simulated or the integer ID.
}
  \item{prior}{
Double vector. The prior on pedigrees.
}
  \item{seed}{
Integer used to fix simulations.
}

  \item{mark}{
Integer. Index of \code{\link[paramlink]{marker}} used in \code{paramlinkConditionalOne}.
}
}
\details{
Marker data is simulated for a specified person and one specified marker 
using \code{paramlinkConditional} which calls \code{\link[paramlink]{markerSim}}. 
The marker data is then loaded into a datamatrix and likelihoods calculated using 'Familias'.
}

\value{
\item{LR.All.Markers }{One LR for each simulation for each pedigree.}
\item{lik.All.Markers }{One likehood for each simulation for each pedigree.}
\item{LR.Per.Marker}{One LR for each simulation and marker for each pedigree.}
\item{lik.Per.Marker}{Onelikelihood for each simulation and marker for each pedigree.}
\item{first.Sim}{NULL}
}

\references{
Kling et al. (2017)
}
\author{
Thore Egeland <Thore.Egeland@gmail.com> and Magnus Dehli Vigeland
}

\seealso{
See also \code{\link{FamiliasConditional}}
}
\examples{

# Example 
data(grandmother)
pedigrees = grandmother$pedigrees
datamatrix = grandmother$datamatrix
loci = grandmother$loci
persons = rownames(datamatrix)
Nsim = 5 #Increase to 1000
res1 = paramlinkConditional(Nsim = Nsim, datamatrix, loci, pedigrees, 
                           truePed = 1, available = "GS", ref=2, seed=17)
LR.H1 = table(res1$LR.All.Markers[,1])/Nsim 
# LR(H1) distribution, agrees well with theory:
# Pr(LR=0.5|H1) = 0.49005, Pr(LR=25.5|H1) = 0.50490, Pr(LR=50.5|H1) = 0.00505
# Simulate from unrelated alternative: 
res2 = paramlinkConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                           truePed = 2, available = "GS", ref=2, seed=17)
LR.H2 = table(res2$LR.All.Markers[,1])/Nsim
# Try mutation 
p = as.double(loci[[1]]$alleles)
loci = list(FamiliasLocus(p, 1:3, "L1", MutationModel = "Proportional",
                          MutationRate=0.005))
x = Familias2linkdat(pedigrees, datamatrix, loci)
res1 = paramlinkConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                           truePed = 1, available = "GS", ref=2, seed=17)
# Simulate father instead
res1 = paramlinkConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                           truePed = 1, available = "FAT", ref=2, seed=17)
LR.H1 = table(res1$LR.All.Markers[,1])/Nsim 

res2 = paramlinkConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                           truePed = 2, available = "FAT", ref=2, seed=17)
LR.H2 = table(res2$LR.All.Markers[,1])/Nsim

### Example 
data(Demo3Markers)
pedigrees = Demo3Markers$pedigrees
datamatrix = Demo3Markers$datamatrix
loci = Demo3Markers$loci
persons = rownames(datamatrix)
Nsim = 5
res1 = paramlinkConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                           truePed = 1, available = "Mother", ref=2, seed=177)
res1[[1]][,1] #Always the same LR! Why?:
x = Familias2linkdat(pedigrees, datamatrix,loci)
P.mother.H1 = oneMarkerDistribution(x[[1]],3,3, verbose=FALSE)
P.mother.H2 = oneMarkerDistribution(x[[2]][[1]],3,3, verbose=FALSE) 
round(P.mother.H1/P.mother.H2,10) == 1
# The probability distribution of the mother is the same for both hypotheses
# and therefore we always get the same LR.

# Example 
data(symmetric)
pedigrees = symmetric$pedigrees
datamatrix = symmetric$datamatrix
loci = symmetric$loci
persons = rownames(datamatrix)
res1 = paramlinkConditional(Nsim = 2, datamatrix, loci, pedigrees, 
                           truePed = 1, available = NULL, ref=2, seed=17)
# Without mutation, all LRs 1. With mutation:
res1$LR.All.Markers

# Example  Mariana's F21 example
\dontrun{ #Takes a few minutes; compares 'Familias' and 'paramlink'.
  #Results and paramlink execution times
  data(F21)
  pedigrees = F21$pedigrees
  datamatrix = F21$datamatrix
  loci = F21$loci
  persons = rownames(datamatrix)
  Nsim = 1000
  start.time <- Sys.time()
  res1.paramlink = paramlinkConditional(Nsim = Nsim, datamatrix, loci, pedigrees, 
                             truePed = 1, available = "Missing Person", ref=2, seed=17)
  res2.paramlink = paramlinkConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                             truePed = 2, available = "Missing Person", ref=2, seed=17)
  end.time <- Sys.time()
  paramlink.time <- end.time - start.time
  start.time <- Sys.time()
  res1.familias = FamiliasConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                              truePed = 1, available = "Missing Person", ref=2, seed=17)
  res2.familias = FamiliasConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                              truePed = 2, available = "Missing Person", ref=2, seed=17)
  end.time <- Sys.time()
  familias.time <- end.time - start.time
  familias.time - paramlink.time #around -0.7 mins
  LR1.familias = res1.familias[[1]][,1]
  LR1.paramlink = res1.paramlink[[1]][,1]
  aa=cbind(LR1.familias, LR1.paramlink)
  foo=apply(aa,2, quantile)
  foo = LR1.familias - LR1.paramlink
  max(foo/(0.5*(LR1.familias+LR1.paramlink)))
  LR1.familias = res2.familias[[1]][,1]
  LR1.paramlink = res2.paramlink[[1]][,1]
  aa=cbind(LR1.familias, LR1.paramlink)
  foo=apply(aa,2, quantile)
  foo = LR1.familias - LR1.paramlink
  max(foo)
 }
}


