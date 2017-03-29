paramlinkConditional = function(Nsim = 5, datamatrix,  loci, pedigrees, 
	truePed = 1, ref=NULL, available = NULL,  prior=NULL, seed = NULL){
	Npeds = length(pedigrees)
	Nmark = length(loci)
	persons = rownames(datamatrix)
	lik = LR = array(1,dim=c(Nmark, Nsim, Npeds))
	set.seed(seed)
	simData = matrix(ncol=2, nrow=length(persons))
	for (s in 1:Nmark){
		foo = paramlinkConditionalOne(Nsim = Nsim, mark=s, datamatrix, loci, pedigrees, 
			truePed = truePed, available = available, seed = NULL, ref=ref)
		lik[s,,] = foo$lik*lik[s,,]
		LR[s,,] = foo$LR*LR[s,,]
		simData = cbind(simData, foo$simDatamatrix[,1:2])
	}
    LRall = likall = matrix(1, nrow=Nsim, ncol = Npeds)
	for (i in 1:Nmark){
		LRall  = LRall*LR[i,,]
		likall = likall*lik[i,,]
	}
    list(LR.All.Markers = LRall, lik.All.Markers = likall, 
		LR.Per.Marker = LR, lik.Per.Marker =lik, first.Sim = simData[,-c(1:2)])
}

paramlinkConditionalOne = function(Nsim = 5, mark=1,  ref=2, datamatrix,
	loci, pedigrees, truePed = 1, available = NULL,	prior=NULL, seed = NULL){
  
	# Simulates markerdata for truePed and marker mark
	x = Familias2linkdat(pedigrees, datamatrix, loci) #fjernet truePed
	target.num = if(is.character(available)) label2num(available,
				pedigrees) else available
	sim1 = markerSim(x[[truePed]], Nsim, available=target.num, partialmarker=mark, verbose=F, seed=seed) # la til
	sim.df = NULL 
	Npeds = length(pedigrees)
	lik = matrix(NA, ncol=Npeds, nrow=Nsim)
	for(i in 1:Npeds){
	  ped = if(inherits(x[[i]], 'linkdat')) removeMarkers(x[[i]]) else
	    lapply(x[[i]], removeMarkers)
	  ped = transferMarkerdata(sim1, ped)
	  lik[,i] = sapply(1:Nsim, function(j) likelihood(ped, j))
	}
  LR = t(apply(lik,1, function(x, ref) x/x[ref], ref))
  list(LR=LR, lik=lik, simData=sim1, 	simDatamatrix=sim.df)
} 