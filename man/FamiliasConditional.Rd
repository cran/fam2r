\name{FamiliasConditional}
\alias{FamiliasConditional}
\alias{FamiliasConditionalOne}
\title{
Conditional simulation of marker data on pedigrees and 'Familias' likelihood
}
\description{
Marker data is simulated for a specified person and several markers using 
\code{FamiliasConditionalOne} for each marker.
The only difference between \code{FamiliasConditional} and 
\code{\link{paramlinkConditional}} is that the former uses
the C implementation of \code{\link[Familias]{FamiliasPosterior}} for likelihood 
calculation while the latter is based on the \code{R} 
implementation in \code{\link[paramlink]{paramlink}}.
}
\usage{
FamiliasConditional(Nsim = 5, datamatrix, loci, pedigrees, truePed = 1, available = NULL,
                  ref = 2, prior = NULL, seed = NULL)
FamiliasConditionalOne(Nsim = 5, mark = 1, ref = 2, datamatrix, 
          persons, loci, pedigrees,truePed = 1, available = NULL, 
          prior = NULL, seed = NULL)
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
An list with elements of type 'FamiliasPedigree'.
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
Integer. Index of \code{\[link[paramlink]{marker}} used in \code{FamiliasConditionalOne}
}
\item{persons}{
Character vector. Names of persons.
}
}
\details{
Marker data is simulated for a specified person and one specified marker using \code{FamiliasConditional} which calls \code{\link[paramlink]{markerSim}} of \code{\link{paramlink}}.
The marker data is then loaded into a datamatrix and likelihoods calculated using Familias.
}

\value{
\item{LR.All.Markers }{One LR for each simulation for each pedigree.}
\item{lik.All.Markers }{One likehood for each simulation for each pedigree.}
\item{LR.Per.Marker}{One LR for each simulation and marker for each pedigree.}
\item{lik.Per.Marker}{Onelikelihood for each simulation and marker for each pedigree.}
\item{first.Sim}{Marker data for first marker.}
}

\references{
Kling et al. (2017)
}
\author{
Thore Egeland <Thore.Egeland@gmail.com> and Magnus Dehli Vigeland
}

\seealso{
See also \code{\link[paramlink]{markerSim}}
}
\examples{
# Example. Mutation.
# See documentation of conditionalLR for examples without mutation
data(grandmother)
pedigrees = grandmother$pedigrees
datamatrix = grandmother$datamatrix
persons = rownames(datamatrix)
Nsim = 100 
p = c(0.5,0.5)
loci = list(FamiliasLocus(p, 1:2, "L1", MutationModel = "Proportional",
                          MutationRate=0.005))
res1 = FamiliasConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                           truePed = 1, available = "GS", ref=2, seed=17)
LR1 = table(res1$LR.All.Markers[,1])/Nsim 
# Next, we calculate genotype probabilities for GS exactly
x = Familias2linkdat(pedigrees[[1]], datamatrix, loci)
m = marker(x, 1:2)
x = addMarker(x,m)
p.GS = oneMarkerDistribution(x, 3, partialmarker=1, verbose = FALSE)
# The difference between simulated and exact probabilities are
LR1 - p.GS[c(2,3,1)]

### Example 
data(Demo3Markers)
pedigrees = Demo3Markers$pedigrees
datamatrix = Demo3Markers$datamatrix
loci = Demo3Markers$loci
x = Familias2linkdat(pedigrees, datamatrix, loci)  
P.mother.H1 = oneMarkerDistribution(x[[1]], ids=3, partialmarker=3,
  ignore.affection=TRUE)
P.mother.H2 = oneMarkerDistribution(x[[2]][[1]], ids=3,
  partialmarker=3, ignore.affection=TRUE)
all(round(P.mother.H1/P.mother.H2,10) == 1)
# The probability distribution of the mother is the same for both hypotheses
# and therefore we always get the same LR.

# Example 
data(symmetric)
pedigrees = symmetric$pedigrees
datamatrix = symmetric$datamatrix
loci = symmetric$loci
persons = rownames(datamatrix)
\dontrun{
# Plot with newdev=TRUE, resize plot window and then plot with newdev=FALSE
x = Familias2linkdat(pedigrees, datamatrix, loci)
plotPedList(x, newdev=FALSE, frametitles=c("H1", "H2", "H3"), 
  available="shaded", marker=1:2, skip.empty.genotypes = TRUE)
}
res1 = FamiliasConditional(Nsim = 2, datamatrix,  loci, pedigrees, 
                truePed = 1, available = "MO", ref = 2, seed = 17)
# Without mutation, all LRs 1. With mutation:
res1$LR.All.Markers
}


