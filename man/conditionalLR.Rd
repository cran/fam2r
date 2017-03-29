\name{conditionalLR}
\alias{conditionalLR}
\alias{LRwrap}
\alias{label2num}
\alias{removePersons}
\title{
Simulates marker data on pedigrees conditionally on typed individuals and calculates likelihoods
}
\description{
Marker data are simulated on pedigrees,
conditional on existing genotypes and likelihoods are calculated
using \code{\link{FamiliasConditional}} or 
\code{\link{paramlinkConditional}}.
\code{transferMarkerdata2} transfers markerdata from a \code{\link[paramlink]{linkdat}}  object to
a \code{\link[paramlink]{linkdat}}  object or a list thereby 
generalising \code{\link[paramlink]{transferMarkerdata}}.
}
\usage{
conditionalLR(Nsim=5, datamatrix, loci, pedigrees, file = NULL , program = "Familias", 
	prior=NULL, available=NULL, seed=NULL, ref=NULL, truePeds = NULL,
	verbose = TRUE, simplify = FALSE)

# Remove persons, e.g., singletons for all hypotheses
removePersons(pedigrees, datamatrix, ids=NULL)

label2num(label, familiasped)
}


\arguments{
  \item{Nsim}{
Integer. Number of simulations.
}
  \item{datamatrix}{
A data frame. The row names must be the names of the persons you have data for. 
The columns contain the alleles, two columns for each marker, in the same order used in the loci list.
}
  \item{loci}{
A FamiliasLocus object or a list of such objects.
}
  \item{pedigrees}{
An object of type 'FamiliasPedigree' or 'pedigree', or a list of such objects.
}
  \item{familiasped}{
An object of type 'FamiliasPedigree' or 'pedigree', or a list of such objects.
}
\item{file}{
Character. First part of name for output file. If NULL nothing is written to files.
}
\item{program}{
Character. Specifies program used for likelihood calculation.
Either 'Familias' (default) or 'paramlink'.
}
  \item{prior}{
Double vector. Not currently relevant as only LRs are reported.
}
  \item{available}{
Character or integer identifying person to be simulated.
}
  \item{seed}{
Integer.
}
  \item{ref}{
Integer. Index of pedigree in numerator of LR. Set to last if \code{NULL}.
}

\item{truePeds}{
Indices of pedigrees to be simulated from. If \code{NULL}, all.
}

\item{simplify}{
Logical. \code{simplify=TRUE} is used for the standard cases,
those with two hypotheses.
}
\item{verbose}{
logical. If \code{TRUE} output is explained.
}
\item{label}{
A character vector, 'Familias' ID label.
}

\item{ids}{
A character or integer vector.
}
}
\value{
The number of hypotheses corresponds to the number of pedigrees. 
In many cases, there will only be two hypotheses. The output is then simplified
if one specifies
\code{simplify = TRUE}. By default the reference hypothesis is number 2. i.e., \code{ref=2})
and the likelihood ratio is LR = \code{Pr(data|H1)/Pr(data|H2)}. 
The output will then be a matrix with columns. The first column
consists of the simulated LR-s when H1 is assumed true, the second one when H2 is true.
When there are more than two hypothes, an array is returned.
\code{LR[,,i]} are the LR values.
when simulations are conditioned on pedigree \code{i}.
There is one row for each simulation and one column for each pedigree.
In other words \code{LR[k,l,i]} is the \code{k}-th simulated value 
of \code{Pr(data|Hl)/Pr(dat|Href)} when hypothesis \code{Hi} is true, the one simulated from.
File(s) with simulated values are written (if variable \code{file} is not \code{NULL}).
If \code{simplify=TRUE} one file is written, otherwise there will be one
one for each pedigree \code{i}. 

}
\details{If \code{truePeds} is a subset of all pedigrees, only the files corresponding to
\code{truePeds} are written. In this case, \code{LR[,,i]} contains missing values if \code{i}
is not in \code{truePeds} }
\references{
Kling et al. (2017).
}
\author{
Thore Egeland <Thore.Egeland@gmail.com>
}

\examples{
# Example
# Computational details for (currently) Example 1 of Kling et al (2017) are provided.
# The purpose is to explain the code and also check results against exact formulae and compare
# to a previous implementation, i.e., library(famr). # There is one marker with alleles 1, 2 
# and 3 having frequencies p1, p2, and p3. # We diregard complicating factors like mutation 
# (see next example), # and silent alleles. # One person, the grandmother GM, is genotyped as 1/1. 
# The grandson GS is to be simulated.

# The dataset 'grandmother' loaded below is a list with three components explained
# in the documentation of FamiliasPosterior, i.e, pedigrees, datamatrix and loci.
# The paramlink function 'Familias2linkdat' converts from 'Familias' format to a format
# suitable for plotting and conditional simulation, likelihood calculations etc, 
# using 'paramlink'; below 'plotPedList' is used to plot.

data(grandmother)
pedigrees = grandmother$pedigrees
datamatrix = grandmother$datamatrix
loci = grandmother$loci
persons = rownames(datamatrix)
\dontrun{
# Plot with newdev=TRUE, resize plot window and then plot with newed=FALSE
x = Familias2linkdat(pedigrees, datamatrix, loci)
plotPedList(x, newdev=TRUE, frametitles=c("H1", "H2"), 
  available ="shaded", marker=1, skip.empty.genotypes = TRUE)
}
# In this case there are two (=length(pedigrees)) 
# hypotheses H1 and H2  or equivalently two pedigrees.
# We will be interested in the likelihood ratio (LR) defined by Pr(data|H1)/Pr(data|H2) 
# or rather the random variables LR(H1) and LR(H2) where  H1 and H2 indicate the true hypotheses.
# In other words, we estimate the distribution of LR(H1) by simulating assuming H1 to be true
# and similarly for H2. Obviously, we also condition on genotyped individuals, GM in this case.
# Assume first H1 to be true. Then, as explained in the Kling et al. (2017), 
# there are three possible values for the likelihood ratios, 
# namely y1=1/2, y2=0.5+1/(4*p1), and y3=0.5+0.5/p1
# occuring with the probabilities py1, py2 and py3 calculated below

p = loci[[1]]$alleles
p1 = p[1]; p2=p[2]; p3=p[3]
py2 = (1-p1)*(p1+0.5)
py3 = 0.5*p1*(1+p1)
py1 =1-py2-py3
y1 = 0.5; y2 = 0.5+1/(4*p1); y3 =0.5+0.5/p1; LRs = c(y1, y2, y3)
LR.H1.exact = c(py1,py2,py3)
names(LR.H1.exact) = paste(LRs)

# The above probability distribution, LR.H1.exact, can be approximated by simulation,
# using 'markerSim' followed by likelihood calculation in 'Familias' or 'paramlink'.

Nsim = 5; seed = 17; avail = "GS"
res1 = conditionalLR(Nsim=Nsim, datamatrix, loci, pedigrees, available=avail,
                     seed=seed, program = "Familias", simplify=TRUE)
res2 = conditionalLR(Nsim=Nsim, datamatrix, loci, pedigrees, available=avail,
                     seed=seed, program = "paramlink",  simplify=TRUE) #Change

LR.H1.Familias = table(res1[,1])/Nsim 
LR.H1.paramlink = table(res2[,1])/Nsim
stopifnot(round(LR.H1.Familias-LR.H1.paramlink,12)==0)
# We see that the two implementations give the same result and agree well
# with the theoretical result. The previous implementation also agrees well.
\dontrun{
  #Try old code
  install.packages("http::/familias.name/famr_1.0-zip")
  library(famr)
  res3 <- conditionalSimulationWrite(nsim = Nsim, datamatrix, persons, 
    loci, pedigrees, available = 3, seed = 1482659436, ref = 2, file = NULL)
  LR.H1.old = table(res3[,1])/Nsim
  LR.H2.old = table(res3[,2])/Nsim
  # The above code is limited in some respects: It assumes 
  # that there are two pedigrees and that there
  # are no mutations. With ref=2, H2 is the denominator of the LR. Column 1 of the output,
  # res3[,1] above is simulated assumed H1 to be true, res[,2] assuming H2 to be true.
}

# Consider next simulation under H2 calculated for H2
LR.H2.Familias = table(res1[,2])/Nsim 
LR.H2.paramlink = table(res2[,2])/Nsim 
stopifnot(round(LR.H2.Familias-LR.H2.paramlink,12)==0)

# Obviously, the possible values for LR are the same for H1 and H2.
# Note that the largest value of LR occurs with probability py3=p1^2= 1e-04
# This value may not be reached in the simulations. Again results agree well and
# also with the previous implementation if Nsim=1000. The probability distribution of LR 
# conditionally on H2 is calculated as

py3 = p1^2
py1 = (1-p1)^2
py2 = 2*p1*(1-p1)
LR.H2.exact = c(py1,py2,py3)
names(LR.H2.exact) = paste(LRs)

# The previous example continues, but we will now 
# model mutations and for simplicity assume a SNP marker. 
\dontrun{
# Let the mutation rate be 0.05
# The mutation rate is chosen (too) high to see some impact.
# in a 'proportional mutation model, i.e.,
p = c(0.2, 0.8); R=0.05
loci = list(FamiliasLocus(p, 1:2, "L1", MutationModel = "Proportional",
                          MutationRate=R))
x = Familias2linkdat(pedigrees[[1]], datamatrix, loci)
m = marker(x,  1, c(1,1), alleles=1:2, afreq=p)
x = addMarker(x,m)
p.GS = oneMarkerDistribution(x, 3, partialmarker=1, verbose = FALSE)
# For instance
p22.one = p.GS["2/2"]
# is the probability of the grandson being 2/2 when
# the grandmother is 1/1; without mutation this would be
p22.ind = 0.5*p[2]^2
# We next check the exact result by an exact formula 
# Egeland, Pinto and Amorim (2017, submitted) and also
# using simulation. Let
H = 1-sum(p^2); k = R/H
p22 = 0.5*p[2]^2*(2-(1-k)^2) # LR = 0.5+0.5*(1-(1-k)^2) with  
stopifnot(round(p22.one-p22, 10)==0)
# probability p22. Below we simulate to check
Nsim = 1000; seed=177
res = conditionalLR(Nsim=Nsim, datamatrix, loci, pedigrees, available="GS",
  seed=seed, program = "Familias", verbose = FALSE, simplify=TRUE)
LR.H1.mut = table(res[,1])/Nsim 
(LR.H1.mut[1]-p22)/p22 #relative difference
}

# Example  Missing grandchild example
\dontrun{ #Takes 3-4 minutes
  data(F21)
  pedigrees = F21$pedigrees
  datamatrix = F21$datamatrix
  loci = F21$loci
  persons = rownames(datamatrix)
  x = Familias2linkdat(pedigrees, datamatrix, loci) 
  Nsim = 1000
  res1 = conditionalLR(Nsim = Nsim, datamatrix,  loci, pedigrees, program ="Familias",
                              available = "Missing Person", seed=17, verbose = FALSE, simplify=TRUE)
  LR = data.frame(LR.H1=res1[,1], LR.H2=res1[,2])
  length(LR[,1][LR[,1]>10^5])/Nsim
  length(LR[,2][LR[,2]==0])/Nsim #PE estimate
  res1 = FamiliasConditional(Nsim = Nsim, datamatrix,  loci, pedigrees,
                              available = "Missing Person", seed=17)
                              
  res2 = conditionalLR(Nsim = Nsim, datamatrix,  loci, pedigrees, program ="Familias",
                       available = "Missing Person", seed=17, verbose = FALSE)
  res = cbind(res1[[1]][,1], res2[,,1][,1])
  
  boxplot(log(res)); title ("log LR(H1), Familias and paramlink (right)")
}

\dontrun{
# Example
data(Demo3Markers)
pedigrees = Demo3Markers$pedigrees
datamatrix = Demo3Markers$datamatrix
loci = Demo3Markers$loci
persons = rownames(datamatrix)
Nsim = 5
res.Familias = conditionalLR(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                              file = NULL, program = "Familias", truePed = NULL, 
                              available = "Mother", ref=NULL, seed=177, simplify = TRUE)
res.paramlink = conditionalLR(Nsim = Nsim, datamatrix,  loci, pedigrees, 
                              file = NULL, program = "paramlink", truePed = NULL, 
                              available = "Mother", ref=NULL, seed=177, simplify=TRUE)

# Always the same LR as the probabilities of the genotypes of the person simulated ("Mother") 
# is the same for both hypotheses:
stopifnot(round(res.Familias[,1]-res.paramlink[1,1],12)==0)
stopifnot(round(res.Familias[,2]-res.paramlink[1,2],12)==0)

# Example 
data(symmetric)
pedigrees = symmetric$pedigrees
datamatrix = symmetric$datamatrix
loci = symmetric$loci
persons = rownames(datamatrix)
truePeds = 1:3
res.Familias = conditionalLR(Nsim = 5, datamatrix,  loci, pedigrees, 
                     file = NULL, program = "Familias", truePeds = truePeds, 
                     available = NULL, ref=2, seed=177, verbose = FALSE)
res.paramlink = conditionalLR(Nsim = 5, datamatrix,  loci, pedigrees, 
                   file = NULL, program = "paramlink", truePeds = truePeds, #Change
                   available = NULL, ref=2, seed=177, verbose = FALSE)
stopifnot(round(res.Familias[,,truePeds]-res.paramlink[,,truePeds],12)==0)

# Example 
ped = list(singleton(5,1), nuclearPed(2))
x = markerSim(ped, N=5, alleles=1:5, verbose=FALSE, available=5)
y = nuclearPed(3)
y2 = transferMarkerdata(x, y)
y2
}
}


