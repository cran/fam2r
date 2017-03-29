LRparamlink = function(x, ref, markers) {
  st = proc.time()
  # get marker names
  y = if(is.linkdat(x[[1]])) x[[1]] else x[[1]][[1]]
  if(missing(markers)) markers = seq_len(y$nMark)
  markernames = sapply(y$markerdata[markers], attr, "name")
  # break all loops (NB: would use rapply, but doesnt work since is.list(linkdat) = TRUE
breaklp = function(a) if(class(a)[[1]] == "linkdat") 
  breakLoops(a,verbose=F) else a # avoid problems with singletons
x_loopfree = lapply(x, function(xx) if(is.linkdat.list(xx))
  lapply(xx, breaklp) else breaklp(xx))

# compute likelihoods
liks = lapply(x_loopfree, function(xx) vapply(markers, function(i)
  likelihood(xx, locus1=i), FUN.VALUE=1))
likelihoodsPerSystem = do.call(cbind, liks)

# LR per marker and total
LRperMarker = do.call(cbind, lapply(1:length(x), function(j)
  liks[[j]]/liks[[ref]]))
LR = apply(LRperMarker,2,prod)

# output
time = proc.time()-st
rownames(likelihoodsPerSystem) = rownames(LRperMarker) = markernames
list(LR=LR, LRperMarker=LRperMarker,
     likelihoodsPerSystem=likelihoodsPerSystem, time=time)
}

