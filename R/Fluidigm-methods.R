setAs('SingleCellAssay', 'FluidigmAssay', function(from)  new("FluidigmAssay",env=from@env,mapping=from@mapping,id=from@id, cellKey=from@cellKey, featureData=from@featureData, phenoData=from@phenoData, cellData=from@cellData, description=from@description))

expavg <- function(x) mean(2^x-1)
logmean <- function(x) log2(mean(x)+1)

#try to throw an error if groups isn't in cellData
#groups can be character vector or symbol (quote(Group1:Group2), used for lattice)
checkGroups <- function(sc, groups){
  if(!is.character(groups) || is.factor(groups))
    stop("'groups' must be character or factor")
  sd <- setdiff(groups, names(cData(sc))) 
  if(length(sd)>0) stop(sprintf('%s not found in %s object', paste(sd, collapse=', '), class(sc)))
  invisible(TRUE)
}

freq <- function(sc){
 stopifnot(inherits(sc, 'SingleCellAssay'))
 apply(exprs(sc)>0, 2, mean)
}

condmean <- function(sc){
stopifnot(inherits(sc, 'SingleCellAssay'))
exprsNA <- exprs(sc)
exprsNA[exprsNA==0] <- NA
apply(exprsNA, 2, mean, na.rm=TRUE)
}

getConcordance <- function(singleCellRef, singleCellcomp, groups=NULL, fun.natural=expavg, fun.cycle=logmean){
  ## vector of groups over which we should aggregate
  ## stopifnot(inherits(singleCellRef, 'FluidigmAssay') && inherits(singleCellcomp, 'FluidigmAssay'))
  mapL <- list(getMapping(singleCellRef), getMapping(singleCellcomp))
  scL <- list(singleCellRef, singleCellcomp)
  castL <- list()
  
  for(i in seq_along(scL)){
    checkGroups(scL[[i]], groups)
    terms1 <- union(groups, mapL[[i]]$ncells)
    lhs1 <- paste(c(terms1, mapL[[i]]$primerid), collapse="+")
    firstForm <- formula(sprintf("%s ~.", lhs1))
    ##should look like Patient.ID + ... + n.cells ~ PrimerID
    tmp <- cast(melt(scL[[i]]), firstForm, fun.aggregate=fun.natural, value=mapL[[i]]$measurement)
    ##exponential average per gene, scaled by number of cells
    tmp["(all)"] <- tmp["(all)"]/tmp[mapL[[i]]$ncells]
    rhs2 <- union(groups, mapL[[i]]$primerid)
    terms2 <- sprintf("%s ~ .", paste(rhs2, collapse="+"))
    secondForm <- formula(terms2)
    nexp = cast(melt(scL[[i]]), secondForm, fun.aggregate=function(x){sum(x>0)}, value=mapL[[i]]$measurement)
    castL[[i]] <- cast(melt(tmp), secondForm, fun.aggregate=fun.cycle)
    renamestr <- c('primerid', 'et')
    names(renamestr) <- c(mapL[[i]]$primerid, '(all)')
    castL[[i]] <- rename(castL[[i]], renamestr)
    castL[[i]] <- cbind(castL[[i]], nexp=nexp$`(all)`)
  }
  concord <- merge(castL[[1]], castL[[2]], by=c(groups, 'primerid'), suffixes=c(".ref", ".comp"), all=T)
  concord[is.na(concord)] <- 0
  concord
}

getwss <- function(concord, nexp){
  mean((concord$et.ref - concord$et.comp)^2*nexp, na.rm=TRUE)
}


getss <- function(concord){
  mean((concord$et.ref - concord$et.comp)^2, na.rm=TRUE)
}

##' Concordance correlation coefficient lin 1989
##' 
##' Return's Lin's intraclass correlation coefficient for concordance
##' @title getrc
##' @param concord data.frame returned by getConcordance
##' @return numeric
getrc <- function(concord){
with(concord, 2*cov(et.ref, et.comp)/(var(et.ref)+var(et.comp)+(mean(et.ref)-mean(et.comp))^2))
}

##' @import lattice
panel.shifts <- function(x, y, groups, subscripts, comparison, ...){
  #vargs = list(...)
  #comparison = vargs$comparison
  #print(str(subscripts))
 panel.xyplot(x, y, subscripts=subscripts, ...)
 panel.segments(x, y, comparison[subscripts,1], comparison[subscripts,2], ..., alpha=.6)
 panel.abline(a=0, b=1, col='gray')
}


concordPlot <- function(concord0, concord1){
  #preconditions concord0, concord1 have row-correspondence
  #concord0 is plotted reference
  p<-xyplot(et.ref ~ et.comp, concord0)
  
}


#TODO: remove multiple cells
filter <- function(sc, groups=NULL, filt_control=NULL, apply_filter=TRUE){

    if (is.null(filt_control)){
    filt_control <- list(filter=T, nOutlier=2, sigmaContinuous=7, sigmaProportion=7, sigmaSum=NULL, K=1.48)
  }

  if (!is.null(groups)) {
    checkGroups(sc, groups)
      scL <- split(sc, groups)
      lapp <- lapply(scL, filter, groups=NULL, filt_control, apply_filter)
      ## Do various things with lapp:
      if(apply_filter){
        ## list of SingleCellAssays
        out <- do.call(combine, lapp)
      } else if(filt_control$filter){
        out <- do.call(rbind, lapp)     #Fix order, argh
        scKey <- split(getcellKey(sc), cData(sc)[,groups])
        scKeybound <- do.call(c, scKey)
        out <- out[order(scKeybound),] #test this
      } else{                           #I'd reckon it's an unapplied filterset, we'll just keep it as a list
        out <- lapp
      }
      return(out)
  }
    
  exprs <- exprs(sc)
  filtered <- do.call(.internalfilter, c(list(exprs), filt_control))
  if(apply_filter){
    anyfilter <- apply(filtered, 1, any)
    scout <- sc[[!anyfilter]]
  } else{
   scout <- filtered 
  }
  scout
}

#rows are cells, columns are genes
#exprs is a matrix of cells x genes for the strata in which the filtering should occur
#filter optionally return a logical vector or the statistics
.internalfilter <- function(exprs, filter=T, nOutlier=2, sigmaContinuous=7, sigmaProportion=7, sigmaSum=NULL,K=1.48){
  exprs[exprs==0]<-NA
  medianC = apply(exprs, 2, median, na.rm=T)
  medianDev = apply(t(abs(t(exprs) - medianC)), 2, median,na.rm=T)
  z.exprs = t((t(exprs) - medianC)/(medianDev*K+.001))
  
  ztot =apply(abs(z.exprs)>sigmaContinuous, 1, sum, na.rm=T)
  
  outlier = ztot >= nOutlier

  if(!is.null(sigmaSum) && !is.na(sigmaSum)){
    zsum <- apply(abs(z.exprs), 1, sum, na.rm=T)
    outlier <- zsum > sigmaSum
  }
  
  fet0 <- apply(!is.na(exprs), 1, mean)
  null<- fet0==0
  fet0 <- asin(sqrt(fet0))
  medianFet0 <- median(fet0[!null], na.rm=T)
  medianDev.fet0 <- median(abs(fet0[!null] - medianFet0), na.rm=T)
  z.fet0 <- (fet0 - medianFet0)/(medianDev.fet0*K + .001)
  
  unlikely = abs(z.fet0) > sigmaProportion


  if(filter){
	  #intout - outlier based on z-score
	  #pctout - outlier based on frequency of expression
    return(data.frame(intout=outlier,null=null, pctout=unlikely))
  } else{
    return(list(z.exprs=z.exprs, fet0=fet0, z.fet0=z.fet0, et.med=medianC, et.mad=medianDev))
  }
}

