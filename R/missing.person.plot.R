missing.person.plot = function(ped_related, missing, id.labels=NULL,
	available="shaded",marker=NULL, width=c(4,4,1), newdev=TRUE,
	frametitles=c("H1: POI related", "H2:POI unrelated"), ...) {
    if(!is.linkdat(ped_related)) stop("Expecting a connected pedigree as H1")

    miss_internal = internalID(ped_related, missing)

    # Labels
    if(is.null(id.labels)) {
        if(!is.null(plotlabs <- ped_related$plot.labels)) {
            stopifnot(length(plotlabs)==ped_related$nInd)
            id.labels = plotlabs
        }
        else
            id.labels = ped_related$orig.ids
    }
    else if(identical(id.labels, "num")) {
        id.labels = ped_related$orig.ids
        id.labels[miss_internal] = "MP"
    }
    else
        id.labels = rep(id.labels, length=ped_related$nInd)

    # Pedigree 1: Related
    labels1 = id.labels
    misslab = labels1[miss_internal]
    labels1[miss_internal] = ifelse(misslab=="", "POI", paste(misslab, "= POI"))

    # Color POI red
    col1 = ifelse(1:ped_related$nInd == miss_internal, 2, 1)

    plot1 = list(ped_related, id.labels=labels1, col=col1)

    # Pedigree 2: Unrelated
    labels2 = id.labels
    #labels2[miss_internal] = "MP"
    plot2 = list(ped_related, id.labels=labels2)

    s = singleton(id=missing, sex = getSex(ped_related, missing))
    if(!is.null(marker))
        s = transferMarkerdata(from=ped_related, to=s)
    plot3 = list(s, id.labels="POI", col=2)

    plotPedList(list(plot1, plot2, plot3), frames=list(1,2:3), available=available,
        marker=marker, skip.empty.genotypes=TRUE,
frametitles=frametitles, newdev=newdev, ...)
}


internalID = function(x, orig.ids){
    internal_ids = match(orig.ids, x$orig.ids)
    if (any(is.na(internal_ids))) 
        stop(paste("Indicated ID(s) not among original ID labels:", 
            paste(orig.ids[is.na(internal_ids)], collapse = ",")))
    internal_ids
}

getSex = function (x, orig.ids) 
as.vector(x$pedigree[internalID(x, orig.ids), "SEX"])