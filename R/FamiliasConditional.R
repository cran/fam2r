FamiliasConditional= function(Nsim = 5, datamatrix, loci, pedigrees, 
	truePed = 1, available = NULL, ref=2, prior=NULL, seed = NULL){
	Npeds = length(pedigrees)
	Nmark = length(loci)
	persons = rownames(datamatrix)
	lik = LR = array(1,dim=c(Nmark, Nsim, Npeds))
	set.seed(seed)
	simData = matrix(ncol=2, nrow=length(persons))
	for (s in 1:Nmark){
		foo = FamiliasConditionalOne(Nsim = Nsim, mark=s, datamatrix, persons, loci, pedigrees, 
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

FamiliasConditionalOne = function(Nsim = 5, mark=1,  ref=2, datamatrix, persons,  
	loci, pedigrees, truePed = 1, available = NULL, prior=NULL, seed = NULL){
	# Simulates markerdata for truePed and marker mark
  x = Familias2linkdat(pedigrees, datamatrix, loci) 
  target.num = if(is.character(available)) label2num(available,pedigrees) else available
  sim1 = markerSim(x[[truePed]], Nsim, available=target.num, 
    partialmarker=mark, verbose=FALSE, seed=seed)
  
  # Fills sim.df, datamatrix, for Familias
  sim.df = if(inherits(sim1, 'linkdat')) as.data.frame(sim1) else
    do.call(rbind, lapply(sim1, as.data.frame))

  persons.ordered = if(inherits(sim1, 'linkdat')) sim1$plot.labels else
    unlist(lapply(sim1,function(z) z$plot.labels))
  rownames(sim.df) = persons.ordered
  sim.df = sim.df[intersect(persons.ordered, persons),-c(1:5)]
  sim.df[sim.df==0]=NA
      
  #One locus made for each simulation
	reploci = rep(loci[mark], Nsim)   
	repnames = paste(loci[[mark]]$locusname, 1:Nsim, sep="_")
	names(reploci) = repnames
	for(i in seq_along(reploci)) reploci[[i]]$locusname = repnames[i]  
    if(is.null(prior)) {
		lp = length(pedigrees)
		prior = rep(1/lp,lp)
	}
	res = FamiliasPosterior(pedigrees, reploci, sim.df, ref=ref, prior=prior)
    list(LR=res$LRperMarker, lik=res$likelihoodsPerSystem, simData=sim1, 
		simDatamatrix=sim.df)
} 