conditionalLR = function(Nsim=5, datamatrix, loci, pedigrees,
    file = NULL , program = "Familias", prior=NULL, available=NULL, seed=NULL, ref=NULL, 
	  truePeds = NULL, verbose = TRUE, simplify = FALSE){
  Npeds = length(pedigrees)
  persons = rownames(datamatrix)
  if(Nsim<1) stop("No of simulations should be an integer 1 or greater")
	if(!is.null(seed)) set.seed(seed)
	if(is.null(ref)) ref = Npeds
	if(is.null(truePeds)) truePeds = 1:Npeds
	LR = array(dim=c( Nsim, Npeds, Npeds))
	if(program!="Familias" & program != "paramlink") 
		stop("variable program must be Familias or paramlink")
	if (program == "Familias")
	  for (i in truePeds)
	   	LR[,,i] = FamiliasConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
  			truePed = i, available = available, ref=ref, prior=prior, seed = seed)[[1]]
	
	if (program == "paramlink")
	  for (i in truePeds)
	    LR[,,i] = paramlinkConditional(Nsim = Nsim, datamatrix,  loci, pedigrees, 
	     truePed = i, available = available, ref=ref, prior=prior, seed = seed)[[1]]
	
	if(simplify) {
	  if(Npeds>2) stop("simplify=TRUE not possible with 3 or more hypotheses")
	  LR = LRwrap(LR, ref=ref)
	  truePeds = 1
	}
	if(!is.null(file) &!simplify) {
	  for (i in truePeds)
	  write.table(LR[,,i], file=paste(file,i,".txt", sep=""), quote=FALSE, row.names=FALSE, 
	              col.names=paste("LR.1",1:Npeds, sep=""))
	}
	if(!is.null(file) & simplify) 
	    write.table(LR, file=file, quote=FALSE, row.names=FALSE)

	if(verbose & !simplify){ 
		cat("LR[,,i] is the likelihood ratio conditioned on pedigree i","\n")
		cat("LR[,,i] is a matrix with one row for each simulation and one column","\n")
		cat("for each pedigree. The denominator of the LR is pedigree no ",paste(ref),"\n")
	}
    LR
}


LRwrap = function(arrayLR, file= NULL, ref=1){
  dims = dim(arrayLR)
  get = (ref==2)*1+(ref==1)*2
  LR = matrix(nrow =dims[1], ncol=dims[2])
  for (s in 1:dims[3])
    LR[,s] = arrayLR[,,s][,get]
  dimnames(LR) <- list(NULL, paste("LR.H", 1:dims[3], sep = ""))
  if (!is.null(file)) 
    write.table(LR, file = file, quote = FALSE, row.names = FALSE, 
                col.names = TRUE)
  LR
}

label2num = function (label, familiasped){
  if (is.list(familiasped) && class(familiasped[[1]]) == "FamiliasPedigree") 
    familiasped = familiasped[[1]]
  match(label, familiasped$id)
}



removePersons = function(pedigrees, datamatrix, ids=NULL){
  if(!is.null(ids)){
    n =length(pedigrees)
    target.num = if(is.character(ids)) label2num(ids,pedigrees) else ids
    for (i in 1:n){
      pedigrees[[i]]$id = pedigrees[[i]]$id[-target.num]
      pedigrees[[i]]$findex = pedigrees[[i]]$findex[-target.num]
      pedigrees[[i]]$findex[pedigrees[[i]]$findex >= target.num] =pedigrees[[i]]$findex[pedigrees[[i]]$findex >= target.num]-1
      pedigrees[[i]]$mindex = pedigrees[[i]]$mindex[-target.num]
      pedigrees[[i]]$mindex[pedigrees[[i]]$mindex >= target.num] =pedigrees[[i]]$mindex[pedigrees[[i]]$mindex >= target.num]-1
      pedigrees[[i]]$sex = pedigrees[[i]]$sex[-target.num]
   }
      datamatrix = datamatrix[-target.num,]
  }
  list(pedigrees=pedigrees, datamatrix=datamatrix)
}
