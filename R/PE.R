PE = function(pedigrees, datamatrix, loci, claim=1, true=2, 
              available = NULL, file= NULL, ignore = FALSE){
  if(ignore){
    for (i in 1:length(loci)){
      alleles = names(loci[[i]]$alleles)
      n=length(alleles)
      M = diag(length(loci[[1]]$alleles))
      loci[[i]]$maleMutationMatrix = M
      loci[[i]]$femaleMutationMatrix = M
  }
}
  x = Familias2linkdat(pedigrees, datamatrix, loci) 
  target.num = if(is.character(available)) label2num(available,pedigrees) else available
  PEall = sapply(1:length(loci), function(i) exclusionPower(ped_claim=x[[claim]],
                  ped_true=x[[true]], ids=target.num, markerindex=i, plot=FALSE)) 
  
  PEcombined = 1- prod(1-PEall)
  PES = c(PEall, PEcombined )
  nn = unlist(c(lapply(loci, function(x) x$locusname), "Combined"))
  PES = data.frame(marker=nn, PE=PES) 
  if(!is.null(file)) {
    write.table(PES, file=file, quote=FALSE, row.names=FALSE)
  }
  PES
}
