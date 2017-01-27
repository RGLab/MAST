##' Exponential average
##'
##' Puts log transformed values onto natural scale and takes mean of vector.
##' Calculates mean(2^x - 1)
##' @param x \code{numeric}
##' @return \code{numeric}
##' @examples
##' x <- 1:10
##' logmean(expavg(x))
##' @export
expavg <- function(x) mean(2^x-1)

##' Log mean
##'
##' Takes mean of natural scaled values and then logrithm
##' Approximately the inverse operation of \code{\link{expavg}}
##' Calculates log2(mean(x) + 1)
##' @param x \code{numeric}
##' @return \code{numeric}
##' @examples
##' x <- 1:10
##' expavg(logmean(x))
##' @export
logmean <- function(x) log2(mean(x)+1)

logshift <- function(x) log2(x+1)

                                        #try to throw an error if groups isn't in cellData
                                        #groups can be character vector or symbol (quote(Group1:Group2), used for lattice)
checkGroups <- function(sc, groups){
    if(!missing(groups) && !is.null(groups)){
        if(!is.character(groups) || is.factor(groups))
            stop("'groups' must be character or factor")
        sd <- setdiff(groups, names(colData(sc)))
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
##' @examples
##' data(vbetaFA)
##' freq(vbetaFA)
##' @export
freq <- function(sc, na.rm=TRUE){
    stopifnot(is(sc, 'SingleCellAssay'))
    apply(exprs(sc)>0, 2, mean, na.rm=na.rm)
}

##' Report the mean et value for each gene
##'
##' NAs are always removed
##' @param sc SingleCellAssay
##' @return vector of means
##' @examples
##' data(vbetaFA)
##' condmean(vbetaFA)
##' @export
condmean <- function(sc){
    stopifnot(is(sc, 'SingleCellAssay'))
    exprsNA <- exprs(sc)
    exprsNA[exprsNA==0] <- NA
    apply(exprsNA, 2, mean, na.rm=TRUE)
}

##' Report standard deviation of et, for positive et for each gene
##'
##' NAs are always removed
##' @param sc SingleCellAssay
##' @return vector of standard deviations
condSd <- function(sc){
    stopifnot(is(sc, 'SingleCellAssay'))
    exprsNA <- exprs(sc)
    exprsNA[exprsNA==0] <- NA
    sqrt(apply(exprsNA, 2, var, na.rm=TRUE))
}

##' Report number of expressing cells per gene
##'
##' NAs are removed
##' @param sc \code{SingleCellAssay}
##' @return \code{numeric} vector
numexp <- function(sc){
    stopifnot(is(sc, 'SingleCellAssay'))
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
##' @seealso plotSCAConcordance
##' @examples
##' data(vbetaFA)
##' sca1 <- subset(vbetaFA, ncells==1)
##' sca100 <- subset(vbetaFA, ncells==100)
##' concord <- getConcordance(sca1, sca100)
##' getss(concord)
##' getrc(concord)
##' @export 
getConcordance <- function(singleCellRef, singleCellcomp, groups=NULL, fun.natural=expavg, fun.cycle=logmean){
    ## vector of groups over which we should aggregate
    if(!(is(singleCellRef, 'FluidigmAssay') && is(singleCellcomp, 'FluidigmAssay'))){
        stop("singleCellRef and singleCellComp should be SingleCellAssay objects")
    }
    scL <- list(singleCellRef, singleCellcomp)
    castL <- list()

    for(i in seq_along(scL)){
        checkGroups(scL[[i]], groups)
        terms1 <- union(groups, "ncells")
        lhs1 <- paste(c(terms1, "primerid"), collapse="+")
        firstForm <- formula(sprintf("%s ~.", lhs1))
        ##should look like Patient.ID + ... + n.cells ~ PrimerID
        m <- as(scL[[i]], 'data.table')
        tmp <- reshape2::dcast(m, firstForm, fun.aggregate=fun.natural)
        ##exponential average per gene, scaled by number of cells
        if(class(m[['ncells']]) == 'factor'){
            warning("ncells is a factor rather than numeric.\n I'll continue, but this may cause problems down the line")
            tmp[['ncells']] <- as.numeric(as.character(tmp[["ncells"]]))
        }
        tmp[["."]] <- tmp[["."]]/tmp[['ncells']]
        rhs2 <- union(groups, "primerid")
        terms2 <- sprintf("%s ~ .", paste(rhs2, collapse="+"))
        secondForm <- formula(terms2)
        nexp = reshape2::dcast(m, secondForm, fun.aggregate=function(x){sum(x>0)}, value.var="value")
                                        #put back on Et scale. fun.cycle adds 1 so -Inf becomes 0 on natural scale
        castL[[i]] <- reshape2::dcast(melt(tmp), secondForm, fun.aggregate=fun.cycle)
        renamestr <- c('.'='et')
        castL[[i]] <- plyr::rename(castL[[i]], renamestr)
        castL[[i]] <- cbind(castL[[i]], nexp=nexp[['.']])
    }
    concord <- merge(castL[[1]], castL[[2]], by=c(groups, 'primerid'), suffixes=c(".ref", ".comp"), all=TRUE)
    concord
}

if(getRversion() >= "2.15.1") globalVariables(c('et.ref', 'et.comp'))

##' @describeIn getConcordance getrc the sum of squares, weighted by nexp
##' @param nexp number of expressed cells per row in \code{concord}
##' @export
getwss <- function(concord, nexp){
    mean((concord$et.ref - concord$et.comp)^2*concord$nexp.ref, na.rm=TRUE)
}

##' @describeIn getConcordance return the sum of squares
##' @export
getss <- function(concord){
    mean((concord$et.ref - concord$et.comp)^2, na.rm=TRUE)
                                        #  log2(mean(((2^concord$et.ref-1)-(2^concord$et.comp-1))^2,na.rm=TRUE))
                                        #  mean(with(concord,{a<-(et.ref>0);b<-(et.comp>0);a*b*(et.ref-et.comp)^2+a*(1-b)+b*(1-a)}),na.rm=TRUE)

}

##' @describeIn getConcordance Return Lin's (1989) concordance correlation coefficient
##' @param concord data.frame returned by getConcordance
##' @export
getrc <- function(concord){
    with(concord, {foo<-na.omit(cbind(et.ref=et.ref,et.comp=et.comp));2*cov(foo[,"et.ref"], foo[,"et.comp"])/(var(foo[,"et.ref"])+var(foo[,"et.comp"])+(mean(foo[,"et.ref"])-mean(foo[,"et.comp"]))^2)})
                                        #concordance on log scale but treating NAs
                                        #    with(concord,{foo<-na.omit(2^cbind(et.ref,et.comp)-1);2*cov(foo[,"et.ref"],foo[,"et.comp"])/(var(foo[,"et.ref"])+var(foo[,"et.comp"])+(mean(foo[,"et.ref"])-mean(foo[,"et.comp"]))^2)})
}



##' Filter a SingleCellAssay
##'
##' Remove, or flag wells that are outliers in discrete or continuous space.
##'
##' The function filters wells that don't pass filtering criteria described in filt_control.
##' filt_control is a list with named elements \code{nOutlier}
##' (minimum nmber of outlier cells for a cell to be filtered [default = 2]
##' \code{sigmaContinuous} (the z-score outlier threshold for the continuous part of the signal) [default = 7]
##' and \code{sigmaProportion} (the z-score outlier threshold for the discrete part of the signal) [default = 7].
##'
##' If \code{groups} is provided, the filtering is calculated within each level of the group, then combined again as output.
##' @param sc The \code{SingleCellAssay} object
##' @param groups An optional \code{character} naming the grouping variable
##' @param filt_control The \code{list} with configuration parameters for the filter.
##' @param apply_filter \code{logical} should the filter be applied, or should a matrix of booleans giving if a well would be subject to a filtering criteria be returned?
##' @return A filtered result
##' @author Andrew McDavid
##' @seealso burdenOfFiltering
##' @examples
##' data(vbetaFA)
##' ## Split by 'ncells', apply to each component, then recombine
##' vbeta.filtered <- filter(vbetaFA, groups='ncells')
##' ## Returned as boolean matrix
##' was.filtered <- filter(vbetaFA, apply_filter=FALSE)
##' ## Wells filtered for being discrete outliers
##' head(subset(was.filtered, pctout))
##' @export
filter <- function(sc, groups=NULL, filt_control=NULL, apply_filter=TRUE){
    default_filt <- list(filter=TRUE, nOutlier=2, sigmaContinuous=7, sigmaProportion=7, sigmaSum=NULL, K=1.48)
    if (is.null(filt_control)){
        filt_control <- list()
    }
    missingControl <- setdiff(names(default_filt), names(filt_control))
    filt_control[missingControl] <- default_filt[missingControl]

    if (!is.null(groups)) {
        checkGroups(sc, groups)
        scL <- split(sc, groups)
        lapp <- lapply(scL, filter, groups=NULL, filt_control=filt_control, apply_filter=apply_filter)
        ## Do various things with lapp:
        if(apply_filter && filt_control$filter){
            ## list of SingleCellAssays
            out <- do.call(cbind, lapp)
        } else if(filt_control$filter){
            out <- do.call(rbind, lapp)     #Fix order, argh
            out <- out[match(getwellKey(sc), row.names(out)),] #test this
        } else{                           #I'd reckon it's an unapplied filterset, we'll just keep it as a list
            out <- lapp
        }
        return(out)
    }

    exprs <- exprs(sc)
    internalfilter<-get(".internalfilter",environment(SingleCellAssay))
    filtered <- do.call(internalfilter, c(list(exprs), filt_control))
    if(apply_filter && filt_control$filter){
        anyfilter <- apply(filtered, 1, any)
        scout <- sc[,!anyfilter]
    } else{
        scout <- filtered
    }

                                        #if filt_control$filter==TRUE & apply_filter==FALSE
                                        #then rename the retured rows using the wellKey mapping
                                        #if(filt_control$filter && !apply_filter){
                                        #  wk<-as.list(names(getwellKey(sc)))
                                        #  names(wk)<-getwellKey(sc)
                                        #  rownames(scout)<-data.frame(idvars=do.call(rbind,wk[rownames(scout)]))[,1]
                                        #}
    scout
}

                                        #rows are cells, columns are genes
                                        #exprs is a matrix of cells x genes for the strata in which the filtering should occur
                                        #filter optionally return a logical vector or the statistics
.internalfilter <- function(exprs, filter=TRUE, nOutlier=2, sigmaContinuous=7, sigmaProportion=7, sigmaSum=NULL,K=1.48){
    exprs[exprs==0]<-NA
    medianC = apply(exprs, 2, median, na.rm=TRUE)
    medianDev = apply(t(abs(t(exprs) - medianC)), 2, median,na.rm=TRUE)
    z.exprs = t((t(exprs) - medianC)/(medianDev*K+.001))

    ztot =apply(abs(z.exprs)>sigmaContinuous, 1, sum, na.rm=TRUE)

    outlier = ztot >= nOutlier

    if(!is.null(sigmaSum) && !is.na(sigmaSum)){
        zsum <- apply(abs(z.exprs), 1, sum, na.rm=TRUE)
        outlier <- zsum > sigmaSum
    }

    fet0 <- apply(!is.na(exprs), 1, mean)
    null<- fet0==0
    fet0 <- asin(sqrt(fet0))
    medianFet0 <- median(fet0[!null], na.rm=TRUE)
    medianDev.fet0 <- median(abs(fet0[!null] - medianFet0), na.rm=TRUE)
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



#' Average within duplicated genes/primers
#'
#' .
#' @param fd \code{SingleCellAssay} or subclass 
#' @param geneGroups \code{character} naming a column in the \code{featureData} that keys the duplicates
#' @param fun.natural transformation to be used to collapse the duplicate expression values
#' @param fun.cycle transformation to be used after collapsing
#' @return collapsed version of \code{fd}.
#' @export
primerAverage <- function(fd, geneGroups, fun.natural=expavg, fun.cycle=logshift){
    fVars <- mcols(fd)[, geneGroups, drop=FALSE]
    geneset <- split(1:nrow(fVars), fVars)
    fdOrder <- order(unlist(lapply(geneset, min)))
    geneset <- geneset[fdOrder]           #maintain order of genes from original
    genefirst <- unlist(lapply(geneset, function(x) x[1]))
    skeleton <- fd[genefirst,]       #make smaller version of fd, then we'll replace exprs with the mean per geneset

    exprs.orig <- assay(fd)
    exprs.new <- lapply(geneset, function(rowidx){
        rows <- exprs.orig[rowidx,, drop=FALSE]
        if(nrow(rows)>1){
            return(fun.cycle(apply(rows, 2,fun.natural)))
        }
        return(rows)
    })
    exprs.new <- do.call(rbind, exprs.new)
    dimnames(skeleton) = dimnames(exprs.new)
    assay(skeleton) <- exprs.new
    rowData(skeleton)$primerid <- mcols(skeleton)[,geneGroups]
    # skeleton <- sort(skeleton,rows=FALSE)
    return(skeleton)
}
