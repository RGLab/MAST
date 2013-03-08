setAs('SingleCellAssay', 'FluidigmAssay', function(from)  new("FluidigmAssay",env=from@env,mapping=addMapping(from@mapping,list(ncells=NULL)),id=from@id, wellKey=from@wellKey, featureData=from@featureData, phenoData=from@phenoData, cellData=from@cellData, description=from@description))

## ###===========Generics===============
## setGeneric('filter')


##' Exponential average
##'
##' Puts log transformed values onto natural scale and takes mean of vector.
##' Calculates mean(2^x - 1)
##' @param x \code{numeric}
##' @return \code{numeric}
##' @export
expavg <- function(x) mean(2^x-1)

##' Log mean
##'
##' Takes mean of natural scaled values and then logrithm
##' Approximately the inverse operation of \code{\link{expavg}}
##' Calculates log2(mean(x) + 1)
##' @param x \code{numeric}
##' @return \code{numeric}
##' @export
logmean <- function(x) log2(mean(x)+1)

logshift <- function(x) log2(x+1)

#try to throw an error if groups isn't in cellData
#groups can be character vector or symbol (quote(Group1:Group2), used for lattice)
checkGroups <- function(sc, groups){
if(!missing(groups) && !is.null(groups)){
  if(!is.character(groups) || is.factor(groups))
    stop("'groups' must be character or factor")
  sd <- setdiff(groups, names(cData(sc)))
  if(length(sd)>0) stop(sprintf('%s not found in %s object', paste(sd, collapse=', '), class(sc)))
 }
  invisible(TRUE)
}

##' Report the proportion of expression for each gene
##'
##' NAs can be optionally removed
##' @param sc SingleCellAssay
##' @param na.rm should NAs be removed, or carried through?
##' @return vector of proportions
freq <- function(sc, na.rm=TRUE){
 stopifnot(inherits(sc, 'SingleCellAssay'))
 apply(exprs(sc)>0, 2, mean, na.rm=na.rm)
}
##' Report the mean et value for each gene
##'
##' NAs are always removed
##' @param sc SingleCellAssay
##' @return vector of means
##' @export
condmean <- function(sc){
stopifnot(inherits(sc, 'SingleCellAssay'))
exprsNA <- exprs(sc)
exprsNA[exprsNA==0] <- NA
apply(exprsNA, 2, mean, na.rm=TRUE)
}

##' Report number of expressing cells per gene
##'
##' NAs are removed
##' @param sc \code{SingleCellAssay}
##' @return \code{numeric} vector
##' @export
numexp <- function(sc){
stopifnot(inherits(sc, 'SingleCellAssay'))
apply(exprs(sc)>0, 2, sum, na.rm=TRUE)
}

##' Get the concordance between two
##'
##' Return the concordance between two assays (i.e. single cell and hundred cell)
##' The "average" of \code{singleCellRef} (after adjusting for the number of cells) and
##' \code{singleCellComp} are taken per gene, per \code{groups}.
##' A \code{data.frame} with one row per gene-\code{groups} is returned with some additional columns.
##' @title getConcordance
##' @param singleCellRef "reference" SingleCellAssay
##' @param singleCellcomp "comparison" SingleCellAssay
##' @param groups character vector giving variable(s) on which the comparison is conditioned
##' @param fun.natural function to transform the SingleCellAssays to a mRNA proportional level
##' @param fun.cycle inverse function of fun.natural
##' @return concordance between two assays
##' @author Andrew McDavid
##' @export getConcordance
##' @importFrom plyr is.formula
##' @importFrom reshape cast
getConcordance <- function(singleCellRef, singleCellcomp, groups=NULL, fun.natural=expavg, fun.cycle=logmean){
  ## vector of groups over which we should aggregate
  if(!(inherits(singleCellRef, 'FluidigmAssay') && inherits(singleCellcomp, 'FluidigmAssay'))){
    stop("singleCellRef and singleCellComp should be SingleCellAssay objects")
  }
  scL <- list(singleCellRef, singleCellcomp)
  castL <- list()

  for(i in seq_along(scL)){
    checkGroups(scL[[i]], groups)
    terms1 <- union(groups, getMapping(scL[[i]],"ncells")[[1]])
    lhs1 <- paste(c(terms1, getMapping(scL[[i]],"primerid")[[1]]), collapse="+")
    firstForm <- formula(sprintf("%s ~.", lhs1))
    ##should look like Patient.ID + ... + n.cells ~ PrimerID
    tmp <- cast(melt(scL[[i]]), firstForm, fun.aggregate=fun.natural, value=getMapping(scL[[i]],"measurement")[[1]])
    ##exponential average per gene, scaled by number of cells
    if(class(melt(scL[[i]])[[getMapping(scL[[i]],"ncells")]]) == 'factor'){
      warning("ncells is a factor rather than numeric.\n I'll continue, but this may cause problems down the line")
    }
    tmp["(all)"] <- tmp["(all)"]/as.numeric(as.character(tmp[[getMapping(scL[[i]],"ncells")[[1]]]]))
    rhs2 <- union(groups, getMapping(scL[[i]],"primerid")[[1]])
    terms2 <- sprintf("%s ~ .", paste(rhs2, collapse="+"))
    secondForm <- formula(terms2)
    nexp = cast(melt(scL[[i]]), secondForm, fun.aggregate=function(x){sum(x>0)}, value=getMapping(scL[[i]],"measurement")[[1]])
    #put back on Et scale. fun.cycle adds 1 so -Inf becomes 0 on natural scale
    #Not sure this is what we want? Okay.. this will be fine
    castL[[i]] <- cast(melt(tmp), secondForm, fun.aggregate=fun.cycle)
    renamestr <- c('primerid', 'et')
    names(renamestr) <- c(getMapping(scL[[i]],"primerid")[[1]], '(all)')
    castL[[i]] <- rename(castL[[i]], renamestr)
    castL[[i]] <- cbind(castL[[i]], nexp=nexp$`(all)`)
  }
  concord <- merge(castL[[1]], castL[[2]], by=c(groups, 'primerid'), suffixes=c(".ref", ".comp"), all=T)
  #This should not be done. Missing values should remain missing, not be artificially set to zero. We'll see what it does to the wss etc..
 # concord[is.na(concord)] <- 0
  concord
}

##' @title getwss
##' @param concord output from \code{getConcordance}
##' @param nexp number of expressed cells per row in \code{concord}
##' @return weighted sum of squares
##' @export getwss
getwss <- function(concord, nexp){
  mean((concord$et.ref - concord$et.comp)^2*nexp, na.rm=TRUE)
#    log2(mean(((2^concord$et.ref-1)-(2^concord$et.comp-1))^2*nexp,na.rm=TRUE))
    ##another WSS function that takes into account the disagreement between wells.
    #mean(with(concord,{a<-(et.ref>0);b<-(et.comp>0);a*b*(et.ref-et.comp)^2+a*(1-b)+b*(1-a)}),na.rm=TRUE)
}

##' @title getss
##' @param concord output from \code{getConcordance}
##' @return sum of squares
##' @export getss
##' @note Ditto.. compute on the natural scale
getss <- function(concord){
  mean((concord$et.ref - concord$et.comp)^2, na.rm=TRUE)
  #  log2(mean(((2^concord$et.ref-1)-(2^concord$et.comp-1))^2,na.rm=TRUE))
#  mean(with(concord,{a<-(et.ref>0);b<-(et.comp>0);a*b*(et.ref-et.comp)^2+a*(1-b)+b*(1-a)}),na.rm=TRUE)

}

##' Concordance correlation coefficient lin 1989
##'
##' Return's Lin's intraclass correlation coefficient for concordance
##' @title getrc
##' @param concord data.frame returned by getConcordance
##' @return numeric
##' @export getrc
getrc <- function(concord){
with(concord, {foo<-na.omit(cbind(et.ref=et.ref,et.comp=et.comp));2*cov(foo[,"et.ref"], foo[,"et.comp"])/(var(foo[,"et.ref"])+var(foo[,"et.comp"])+(mean(foo[,"et.ref"])-mean(foo[,"et.comp"]))^2)})
    #concordance on log scale but treating NAs
#    with(concord,{foo<-na.omit(2^cbind(et.ref,et.comp)-1);2*cov(foo[,"et.ref"],foo[,"et.comp"])/(var(foo[,"et.ref"])+var(foo[,"et.comp"])+(mean(foo[,"et.ref"])-mean(foo[,"et.comp"]))^2)})
}

## This might have been made obsolete by plotSCAConcordance
## This is a panel function that draws lines between (x, y) pairs and comparison
## comparison is a data.frame output from getConcordance
panel.shifts <- function(x, y, groups, subscripts, comparison, ...){
  #vargs = list(...)
  #comparison = vargs$comparison
  #print(str(subscripts))
 panel.xyplot(x, y, subscripts=subscripts, ...)
 panel.segments(x, y, comparison[subscripts,1], comparison[subscripts,2], ..., alpha=.5)
 panel.abline(a=0, b=1, col='gray')
}



concordPlot <- function(concord0, concord1){
  #preconditions concord0, concord1 have row-correspondence
  #concord0 is plotted reference
  p<-xyplot(et.ref ~ et.comp, concord0)
}


###TODO: remove multiple cells,
###make this S3 generic so we don't clobber filter in R
##' Function that filters a single cell assay
##'
##' The function filters wells that don't pass filtering criteria described in filt_control. filt_control is a list with named elements nOutliers (minimum nmber of outlier cells for a cell to be filtered. sigmaContinuous (the z-score outlier threshold for the continuous part of the signal), and sigmaProportion (the z-score outlier threshold for the discrete part of the signal).
##' @title Filter a SingleCellAssay or Fluidigm Assay
##' @param sc The \code{SingleCellAssay} object
##' @param groups The \code{character} naming the grouping variable (optional)
##' @param filt_control The \code{list} with configuration parameters for the filter.
##' @param apply_filter \code{logical} should the filter be applied?
##' @return A filtered result
##' @author Andrew McDavid
##' @export filter
filter <- function(sc, groups=NULL, filt_control=NULL, apply_filter=TRUE){
  default_filt <- list(filter=T, nOutlier=2, sigmaContinuous=7, sigmaProportion=7, sigmaSum=NULL, K=1.48)
    if (is.null(filt_control)){
    filt_control <- list()
  }
  missingControl <- setdiff(names(default_filt), names(filt_control))
  filt_control[missingControl] <- default_filt[missingControl]

  if (!is.null(groups)) {
    checkGroups(sc, groups)
      scL <- split(sc, groups)
      lapp <- lapply(scL, filter, groups=NULL, filt_control, apply_filter)
      ## Do various things with lapp:
      if(apply_filter && filt_control$filter){
        ## list of SingleCellAssays
        out <- .SingleCellAssayCombine(lapp)
      } else if(filt_control$filter){
        out <- do.call(rbind, lapp)     #Fix order, argh
        scKey <- split(getwellKey(sc), cData(sc)[,groups])
        scKeybound <- do.call(c, scKey)
        out <- out[match(getwellKey(sc), scKeybound),] #test this
      } else{                           #I'd reckon it's an unapplied filterset, we'll just keep it as a list
        out <- lapp
      }
      return(out)
  }

  exprs <- exprs(sc)
  filtered <- do.call(SingleCellAssay:::.internalfilter, c(list(exprs), filt_control))
  if(apply_filter && filt_control$filter){
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

##' Estimate Gene Frequency from experiments with pools with varying number-of-cells
##'
##' We use maximum likelihood estimation of a censored binomial distribution
##' We return the standard error, as determined from the observed Fisher Information
##' @title Pooled Gene Frequency Estimation
##' @param fd FluidigmAssay object with ncells set appropriately
##' @return Matrix with point estimates of frequency and standard deviation of estimate
##' @author andrew
freqFromPools <- function(fd){
  stopifnot(inherits(fd, 'FluidigmAssay'))
ncellID <- getMapping(fd, 'ncells')
geneID <- getMapping(fd, 'primerid')
genes <- fData(fd)[,geneID]
ee <- exprs(fd)
ncells <- cData(fd)[,ncellID]
est.and.se <- apply(ee, 2, function(col){
    est <- optimize(pooledModel, c(0, 1), maximum=TRUE, observed=col, ncells=ncells)$max
    se <- 1/sqrt(pooledModel(est, col, ncells, info=TRUE, loglik=FALSE))
    c(estimate=est, standard.err=se)
})
est.and.se
}

##' Summarize expression parameters
##'
##' The mean of positive cells, mu, proportion of gene expression pi,
##' and number of expressing cells per \code{groups} per gene is returned as a list
##' @param fd object inheriting from SingleCellAssay
##' @param groups character vector of grouping variables
##' @return list of mu, pi and num.
setMethod('summary', 'SingleCellAssay', function(object, groups){
  if(!missing(groups) && !is.null(groups)){
checkGroups(object, groups)
sp <- split(object, groups)
mu <- lapply(sp, condmean)
pi <- lapply(sp, freq)
num <- lapply(sp, numexp)
mu <- t(do.call(rbind, mu))
pi <- t(do.call(rbind, pi))
num <- t(do.call(rbind, num))
} else{
  mu <- data.frame(condmean(object))
pi <- data.frame(freq(object))
num <- data.frame(numexp(object))
}
list(mu=mu, pi=pi, num=num)
})


primerAverage <- function(fd, geneGroups, fun.natural=expavg, fun.cycle=logshift){
  fVars <- fData(fd)[, geneGroups, drop=FALSE]
  geneset <- split(1:nrow(fVars), fVars)
  fdOrder <- order(unlist(lapply(geneset, min)))
  geneset <- geneset[fdOrder]           #maintain order of genes from original
  genefirst <- unlist(lapply(geneset, function(x) x[1]))
  skeleton <- fd[[TRUE,genefirst]]       #make smaller version of fd, then we'll replace exprs with the mean per geneset

  exprs.orig <- exprs(fd)
  exprs.new <- lapply(geneset, function(colidx){
    cols <- exprs.orig[,colidx, drop=FALSE]
    if(ncol(cols)>1){
    return(fun.cycle(apply(cols, 1,fun.natural)))
  }
    return(cols)
  })
  exprs.new <- do.call(cbind, exprs.new)
  exprs(skeleton) <- exprs.new
  m <- addMapping(getMapping(fd), list(primerid=geneGroups), replace=TRUE)
  skeleton@mapping <- m
  skeleton
  ## TODO: fix primerids
}
