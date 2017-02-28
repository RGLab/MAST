##' Construct a SingleCellAssay from a matrix or array of expression
##'
##' If the gene expression measurements are already in a rectangular form,
##' then this function allows an easy way to construct a SingleCellAssay object while
##' still doing some sanity checking of inputs.
##' @param class desired subclass of object.  Default \code{SingleCellAssay}.
##' @param exprsArray matrix or array, columns are cells, rows are genes
##' @param cData cellData an object that can be coerced to a DataFrame, ie, data.frame, AnnotatedDataFrame.  Must have as many rows as \code{ncol(exprsArray)}
##' @param fData featureData an object that can be coerced to a DataFrame, ie, data.frame, AnnotatedDataFrame.  Must have as many rows as \code{nrow(exprsArray)}.
##' @return an object of class \code{class}
##' @export
##' @examples
##' ncells <- 10
##' ngenes <- 5
##' fData <- data.frame(primerid=LETTERS[1:ngenes])
##' cData <- data.frame(wellKey=seq_len(ncells))
##' mat <- matrix(rnorm(ncells*ngenes), nrow=ngenes)
##' sca <- FromMatrix(mat, cData, fData)
##' stopifnot(inherits(sca, 'SingleCellAssay'))
##' stopifnot(inherits(sca, 'SummarizedExperiment'))
##' ##If there are mandatory keywords expected by a class, you'll have to manually set them yourself
##' cData$ncells <- 1
##' fd <- FromMatrix(mat, cData, fData)
##' stopifnot(inherits(fd, 'SingleCellAssay'))
FromMatrix <- function(exprsArray, cData, fData, class='SingleCellAssay'){
    can <- checkArrayNames(exprsArray, cData, fData)
    nslice <- dim(can$exprsArray)[3]
    assays <- vector('list', length=nslice)
    for(i in seq_len(nslice)){
        assays[[i]] <- can$exprsArray[,,i,drop=FALSE]
        dim(assays[[i]]) <- dim(assays[[i]])[-3] # only drop last index
        dimnames(assays[[i]]) <- dimnames(can$exprsArray)[-3]
    }
    names(assays) <- dimnames(can$exprsArray)[3]
    obj <- SummarizedExperiment(assays=assays, colData=as(can$cData, 'DataFrame'))
    mcols(obj) <- as(can$fData, 'DataFrame')
    as(obj, class)
}

as3dArray <- function(matOrArray){
    dn <- dimnames(matOrArray)
    dm <- dim(matOrArray)
    if(length(dm)>3 || length(dm)<2 ) stop('`exprsArray` must be matrix or 3-d array')
    if(length(dm)==3) return(matOrArray)
    ## length(dm)==2
    dim(matOrArray) <- c(dim(matOrArray), 1)
    dimnames(matOrArray) <- c(dn, list(NULL))
    return(matOrArray)
}

checkArrayNames <- function(exprsArray, cData, fData){
    if(!is.numeric(exprsArray)) stop('`exprsArray` must be numeric')
    exprsArray <- as3dArray(exprsArray)
    dn <- dimnames(exprsArray)
    noDimnames <- is.null(dn) || is.null(dn[[1]]) || is.null(dn[[2]])

    pidDefault <- if(is.null(dn[[1]])) sprintf('p%0*d', ceiling(log10(nrow(exprsArray)+1)), seq_len(nrow(exprsArray))) else dn[[1]]
    wkDefault <- if(is.null(dn[[2]])) sprintf('wk%0*d', ceiling(log10(ncol(exprsArray)+1)), seq_len(ncol(exprsArray))) else dn[[2]]
    
    if(missing(fData)) fData <- S4Vectors::DataFrame(primerid=pidDefault)
    if(missing(cData)) cData <- S4Vectors::DataFrame(wellKey=wkDefault)
    
    
    if(ncol(exprsArray) != nrow(cData)) stop('`cData` must contain as many rows as `exprsArray`')
    if(nrow(exprsArray) != nrow(fData)) stop('`fData` must contain as many columns as `exprsArray`')

    if(!('primerid' %in% names(fData))){
        message("`fData` has no primerid.  I'll make something up.")
        fData[['primerid']] <- pidDefault
    } else{
        fData[['primerid']] <- as.character(fData$primerid)
    }
    row.names(fData) <- fData$primerid

    if(!('wellKey' %in% names(cData))){
        message("`cData` has no wellKey.  I'll make something up.")
        cData[['wellKey']] <- wkDefault
    } else{
        cData[['wellKey']] <- as.character(cData$wellKey)
    }
    row.names(cData) <- cData$wellKey
    

    if(noDimnames){
        message('No dimnames in `exprsArray`, assuming `fData` and `cData` are sorted according to `exprsArray`')
        dn <- list(primerid=row.names(fData), wellkey=row.names(cData), 'et')  
    }
    
    if(!isTRUE(all.equal(dn[[2]], cData$wellKey))) stop('Order of `exprsArray` and `cData` doesn\'t match')
    if(!isTRUE(all.equal(dn[[1]], fData$primerid))) stop('Order of `exprsArray` and `fData` doesn\'t match')
    dimnames(exprsArray) <- dn
    list(exprsArray=exprsArray, cData=cData, fData=fData)
}

setMethod('fData', 'SingleCellAssay', function(object){
    .Deprecated('mcols')
    mcols(object)
})



##' "Melt" a \code{SingleCellAssay} matrix
##'
##' Return a molten (flat) representation, taking the
##' cross-product of the expression values, the \code{colData} (column meta data),
##' and the feature data (\code{mcols}).
##' @param data \code{SingleCellAssay}
##' @param ... ignored
##' @param na.rm ignored
##' @param value.name name of 'values' column in returned value
##' @return A \code{data.table}, with the cartesian product of the
##' row and column attributes and the expression values
##' @examples
##' data(vbetaFA)
##' melt.SingleCellAssay(vbetaFA[1:10,])
##' as(vbetaFA[1:10,], 'data.table')
melt.SingleCellAssay<-function(data,...,na.rm=FALSE, value.name='value'){
    featdata <- as.data.table(mcols(data))
    celldata <- as.data.table(colData(data))

    exprs <- do.call(cbind, lapply(assays(data), as.vector))
    if( (ncol(exprs)==1 && !is.null(value.name)) || is.null(colnames(exprs))|| colnames(exprs)=='NULL'){ #??
        colnames(exprs) <- make.unique(rep(value.name, ncol(exprs)))
    }
    
    m <- cbind( featdata[rep(seq_len(nrow(featdata)), nrow(celldata)),,drop=FALSE],
               celldata[rep(seq_len(nrow(celldata)), each=nrow(featdata)),,drop=FALSE], exprs)

    
    return(m)
}


mkunique<-function(x,G){
    cbind(x,primerid.unk=make.unique(as.character(get(G,x))))
}


check.vars <- function(cellvars, featurevars, phenovars, dataframe, nc, nr){
    if( any(cellvars %in% featurevars))
        stop("'cellvars', 'idvars' must be disjoint from 'featurevars', 'primerid', 'geneid'")
    cvars.in <- cellvars %in% names(dataframe)
    fvars.in <- featurevars %in% names(dataframe)
    if( !all(cvars.in)) stop(cellvars[!cvars.in][1], ' not found')
    if( !all(fvars.in)) stop(featurevars[!fvars.in][1], ' not found')
    nuniquef<-nrow(uniqueModNA(dataframe[,featurevars,with=FALSE], 'primerid'))
    nuniquec<-nrow(uniqueModNA(dataframe[,cellvars,with=FALSE], 'wellKey'))  
    if(nuniquef != nc)
        stop("'featurevars' must be keyed by 'primerid'")
    if(nuniquec != nr)
        stop("'cellvars' must be keyed by 'idvars'")
}

if(getRversion() >= "2.15.1") globalVariables(c(
                                  'primerid.orig',
                                  'wellKey', 'id')) #fixdf

## might have bad complexity, but could construct one at time, then glue cheaply
## Not too bad except for deduplication.. will use data.table
## Output is sorted by primerid then wellKey
##' @importFrom data.table melt  :=  setkey  setkeyv  set %like%  dcast  data.table  rbindlist  setDT  CJ  .SD  melt  like  setorder  setnames  .N  setDF key setorderv dcast.data.table melt.data.table setattr as.data.table
fixdf <- function(df, idvars, primerid, measurement, cmap, fmap){
    df<-data.table(df)
    cn_df<-colnames(df)
    if(!is(df,"data.frame")){
        stop("Argument `dataframe` should be a data.frame.")
    }
    if(!all(idvars %in% cn_df)){
        stop("Invalid idvars column name. Not in data.frame")
    }
    if(!all(primerid %in% cn_df)){
        stop("Invalid primerid column name. Not in data.frame")
    }
    if(!all(measurement %in% cn_df)){
        stop("Invalid measurement column name. Not in data.frame")
    }
    ## FIXME: should check if cmap and fmap are in df and throw intelligible error

    bothMap <- c(cmap, fmap)
    for(nm in names(bothMap)){
        if(nm %in% cn_df && bothMap[[nm]] != nm){
            warning("renaming column ", nm)
            setnames(df,cn_df,make.unique(c(nm, colnames(df)))[-1])
        }
        if(!all(bothMap[[nm]] %in% names(df))){
            stop('could not find column named ', bothMap[[nm]], ' in dataframe')
        }
        set(df, j=nm, value=df[,bothMap[[nm]], with=FALSE])
                                        #setnames(df, bothMap[[nm]],nm)
    }

    wk <- do.call(paste, df[,idvars,with=FALSE])
    pid <- do.call(paste, df[,primerid, with=FALSE])
    if(any(is.na(wk))) warning('Dropping NAs from wellKey')
    if(any(is.na(pid))) warning('Dropping NAs from primerid')
    
                                        #df[,'wellKey'] <- wk
    df[,wellKey:=wk]
                                        #df[,'primerid'] <- pid
    df[,primerid:=pid]
    dupPrimers <- table(df$wellKey, df$primerid) #cross tab of primerid x wellKey
    duped.primers <- apply(dupPrimers>1, 2, which)
    duped.primers <- duped.primers[sapply(duped.primers, length)>0]
    incomplete <- any(dupPrimers==0)
    if(length(duped.primers)>0){
        warning("Primerid ", names(duped.primers)[1], " appears be duplicated.\n I will attempt to make it unique, but this may fail if the order of the primers is inconsistent in the dataframe.")
                                        #dt$primer.orig <- dt$primerid
        df[,primerid.orig:=primerid]
        df[,primerid:=make.unique(.SD$primerid.orig),by='wellKey']
                                        #    df <- ddply(df,'wellKey',mkunique,G='primerid')
                                        #    df[,'primerid.orig'] <- df[,'primerid']
                                        #    df[,'primerid'] <- df[,'primerid.unk']
                                        #    df[,'primerid.unk'] <- NULL
        fmap['primerid.orig'] <- 'primerid.orig'
    }
    
    if(incomplete){
        message("dataframe appears incomplete, attempting to complete it with NAs")
        skeleton <- data.table((expand.grid(unique(df[,primerid]), unique(df[, wellKey]),stringsAsFactors=FALSE)))
        setnames(skeleton,c("primerid","wellKey"))
        setkey(skeleton,primerid,wellKey);
        setkey(df,primerid,wellKey);
        df<-df[skeleton]
    }

                                        #ord <- do.call(order, df[, c("primerid", "wellKey")])
                                        #df <- df[ord,]
    keynames <- c('primerid', 'wellKey')
    keynames <- union(keynames, colnames(df))
    setkeyv(df,keynames)
    list(df=(df), rn=unique(df$wellKey), cn=unique(df$primerid), fmap=fmap, cmap=cmap)
}

##' Construct a SingleCellAssay (or derived subclass) from a `flat` (melted) data.frame/data.table
##'
##' SingleCellAssay are a generic container for such data and are simple wrappers around SummarizedExperiment objects.
##' Subclasses exist that embue the container with additional attributes, eg \link{FluidigmAssay}.
##' @param dataframe A 'flattened' \code{data.frame} or \code{data.table} containing columns giving cell and feature identifiers and  a measurement column
##' @param idvars character vector naming columns that uniquely identify a cell
##' @param primerid character vector of length 1 that names the column that identifies what feature (i.e. gene) was measured
##' @param measurement character vector of length 1 that names the column containing the measurement 
##' @param id An identifier (eg, experiment name) for the resulting object
##' @param cellvars Character vector naming columns containing additional cellular metadata
##' @param featurevars Character vector naming columns containing additional feature metadata
##' @param phenovars Character vector naming columns containing additional phenotype metadata
##' @param class character providing desired subclass to construct.
##' @param ... additional arguments are ignored
##' @export
##' @aliases SingleCellAssay
##' @aliases FluidigmAssay
##' @examples
##' data(vbeta)
##' colnames(vbeta)
##' vbeta <- computeEtFromCt(vbeta)
##' vbeta.fa <- FromFlatDF(vbeta, idvars=c("Subject.ID", "Chip.Number", "Well"),
##' primerid='Gene', measurement='Et', ncells='Number.of.Cells',
##' geneid="Gene",cellvars=c('Number.of.Cells', 'Population'),
##' phenovars=c('Stim.Condition','Time'), id='vbeta all', class='FluidigmAssay')
##' show(vbeta.fa)
##' nrow(vbeta.fa)
##' ncol(vbeta.fa)
##' head(mcols(vbeta.fa)$primerid)
##' table(colData(vbeta.fa)$Subject.ID)
##' vbeta.sub <- subset(vbeta.fa, Subject.ID=='Sub01')
##' show(vbeta.sub)
##' @return SingleCellAssay, or derived, object
FromFlatDF<-function(dataframe,idvars,primerid,measurement,id=numeric(0), cellvars=NULL, featurevars=NULL, phenovars=NULL, class='SingleCellAssay', ...){
    if(missing(dataframe) || missing(idvars) || missing(primerid) || missing(measurement)){
        stop("Must supply all of 'idvars', 'primerid', 'measurement'  and 'dataframe'")
    }
    if(!is(dataframe,"data.table")){
        dataframe<-data.table(dataframe)
    }
    setkeyv(dataframe, colnames(dataframe))
    ## fixdf: make primerid unique, generate idvar column, rename columns according to cmap and fmap, complete df
    ## keeping support for "mappings" in case we change our mind

    defaultcls <- new(class)
    fmap <- defaultcls@fmap
    cmap <- defaultcls@cmap
    dots <- list(...)
    for(a in names(dots)){
        if(a %in% names(fmap)){
            fmap[a] <- dots[[a]]
        }
        if(a %in% names(cmap)){
            cmap[a] <- dots[[a]]
        }
    }
    
    fixed <- fixdf(dataframe, idvars, primerid, measurement,
                   cmap=cmap, fmap=fmap)
    dl <- array(as.matrix(fixed$df[,measurement, with=FALSE]),
                dim=c(length(fixed$rn), length(fixed$cn), length(measurement)),
                dimnames=list(wellKey=fixed$rn, primerid=fixed$cn, layer=measurement))
    dl <- aperm(dl, c(2, 1, 3))
    cellvars <- unique(c(cellvars, names(fixed$cmap), idvars, phenovars, 'wellKey'))
    featurevars <- unique(c(names(fixed$fmap), primerid, featurevars, 'primerid'))
    check.vars(cellvars, featurevars, phenovars, fixed$df, length(fixed$cn), length(fixed$rn))
    cell.adf<-DataFrame(uniqueModNA(fixed$df[,cellvars,with=FALSE], 'wellKey') )
    row.names(cell.adf) <- unique(fixed$df$wellKey)
    ##pheno.adf <- new('AnnotatedDataFrame')
    ##need a phenokey into the melted data frame for this to make sense
    f.adf <- DataFrame(uniqueModNA(fixed$df[,featurevars, with=FALSE], 'primerid'))
    row.names(f.adf) <- unique(fixed$df$primerid)
    FromMatrix(dl, cell.adf, f.adf, class)
}

#'@export
#'@rdname FromFlatDF
FluidigmAssay <- SingleCellAssay <- function(...){
    .Deprecated("FromFlatDF")
    FromFlatDF(...)
}

uniqueModNA.old <- function(df, exclude){
                                        #browser()
    df <- as.data.frame(df)
    w.include <- names(df)
    if(ncol(df)>1){
        w.include <- setdiff(w.include, exclude)
    }
    u <- unique(df)
    allNa <- apply(is.na(u)[,w.include, drop=FALSE], 1, all)
    u[!allNa,,drop=FALSE]
}

## Get unique rows in data.frame, only counting NAs as distinct for
## columns named in `include`
## Precondition: keyed data.table
## Result will be sorted by exclude
uniqueModNA <- function(df, include){
    if(!is(df, 'data.table')){
        stop('df should be data.table')
    }
    k <- key(df)
    setkey(df,NULL)                     # so that unique operates on all columns
    setorderv(df, k)
    w.exclude <- names(df)
    if(ncol(df)>1){
        w.exclude <- setdiff(w.exclude, include)
    }
    ##u <- unique(df)
    ##anyNa <- apply(is.na(u)[,w.include, drop=FALSE], 1, all)
    ##u[!anyNa,,drop=FALSE]
    u<-unique(df)
    allNa <- apply(is.na(u)[,w.exclude, drop=FALSE], 1, all)
    u<-u[!allNa,] #either we spend time above or here..
    setorderv(u, include)
    u
}


setMethod('getwellKey', 'SingleCellAssay', function(sc) {colData(sc)$wellKey})

##' @rdname cData
##' @export
setMethod('cData', 'SingleCellAssay', function(sc){
    .Deprecated('use colData')
    colData(sc)
})

##' Subset a \code{SingleCellAssay} by cells (columns)
##'
##' Evaluates the expression in \code{...} in the context of \code{colData(x)} and returns a subsetted version of \code{x}
##' @param x \code{SingleCellAssay}
##' @param ... expression
##' @return \code{SingleCellAssay}
##' @examples
##' data(vbetaFA)
##' subset(vbetaFA, ncells==1)
##' @export
setMethod('subset', 'SingleCellAssay', function(x, ...){
    e <- substitute(...)
    asBool <- try(eval(e, colData(x), parent.frame(n=1)), silent=TRUE)
    if(is(asBool, 'try-error')) stop(paste('Variable in subset not found:', strsplit(asBool, ':')[[1]][2]))
                                        #this is a special case of "subset", not of the "[[" method, so..
    if(isTRUE(asBool)){
        x 
    }else{
        x[,asBool]
    }
})



setReplaceMethod('assayNames', c('SingleCellAssay', 'character'), function(x, i, ..., value){
    an <- assayNames(x)
    if(is.null(an)) an <- rep(NA_character_, length(assays(x)))
    if(missing(i)) i <- seq_along(assays(x))
    if(any(i != floor(i))) stop("`i` must be an integer")
    if(any(i<1 | i>length(assays(x)))) stop("`i` out of bounds")
    an[i] <- value
    names(assays(x, withDimnames=FALSE)) <- an
    x
})


##' @export
##' @rdname cData
setReplaceMethod("cData", "SingleCellAssay", function(sc, value) {
    .Deprecated('colData<-')
    colData(sc) <- value
    sc
})

##' Replace \code{colData}
##'
##' Replace \code{colData} with a \code{DataFrame}.
##' Checks to make sure that \code{row.names(value)} match \code{colnames{x}}, in contrast to the parent method
##'  Checks for a wellKey column, as well.
##' @param x \code{SingleCellAssay}
##' @param value \code{DataFrame}
##' @return modified \code{SingleCellAssay}
##' @export
setReplaceMethod("colData", c("SingleCellAssay", 'DataFrame'), function(x, value) {
    ## Only reason we over-ride
    ## Parent doesn't make this test
    if( (nrow(value) != ncol(x)) || any(row.names(value) != colnames(x))) stop('`row.names` in replacement value mismatch colnames in `x`')
    ##I guess?? Do we want/need a wellKey column anymore?
    if(!('wellKey' %in% names(value))) stop('Replacement value must contain a "wellKey" column')
    
    callNextMethod()
})



##' Split into \code{list}
##'
##' Splits a \code{SingleCellAssay} into a \code{list} by a factor (or something coercible into a factor) or a character giving a column of \code{colData(x)}
##' @param x SingleCellAssay
##' @param f length-1 character, or atomic of length ncol(x)
##' @param drop drop unused factor levels
##' @param ... ignored
##' @return List
##' @examples
##' data(vbetaFA)
##' split(vbetaFA, 'ncells')
##' fa <- as.factor(colData(vbetaFA)$ncells)
##' split(vbetaFA, fa)
##' 
##' @aliases split,SingleCellAssay,factor-method split,SingleCellAssay,list-method
##' @export
setMethod('split', signature(x='SingleCellAssay', f='character'), 
          function(x, f, drop=FALSE, ...){
              ## Split a SingleCellAssay by criteria
###f must be a character naming a cData variable
              if(length(f) != ncol(x)){
                  if(any(notin <- !(f %in% names(colData(x))))) stop('Variable ', f[notin], ' not found in `colData(x)`')
                  f <- lapply(colData(x)[,f, drop=FALSE], as.factor)
              } else{
                  f <- as.factor(f)
              }
              split(x, f, drop=drop)
          })

setMethod('split', signature(x='SingleCellAssay', f='factor'), function(x, f, drop=FALSE, ...){
    split(x, list(f))
})


setMethod('split', signature(x='SingleCellAssay', f='list'), function(x, f, drop=FALSE, ...){
    fi <- do.call(interaction, f)[,drop=drop]
    fidx <- split(seq_len(ncol(x)), f)
    ## TODO: could use a CompressedList here, might be more efficient.
    setNames(lapply(seq_along(fidx), function(i) x[,fidx[[i]]]), levels(fi))
})


##' @export
##' @rdname cData
##' @param x \code{SingleCellAssay}
##' @param y \code{SingleCellAssay}
##' @param ... \code{SingleCellAssay}
##' @details  \code{combine(x, y, ...)}: Concatenate two experiments along rows/columns
setMethod('combine', signature(x='SingleCellAssay', y='SingleCellAssay'), function(x, y,  ...){
    ## Not using .Deprecated because I want to maintain test coverage for this function.
    warning('Deprecated: use rbind/cbind')
    if(ncol(x) == ncol(y) ){
        do.call(rbind, list(x, y, ...))
    } else if(nrow(x) == nrow(y)){
        do.call(cbind, list(x, y, ...))
    } else{
        stop("Neither row nor column dimensions match")
    }
})

##' @export
##' @rdname cData
setMethod('combine', signature(x='SingleCellAssay', y='ANY'), function(x, y,  ...){
    if(ncol(x) == nrow(y) || ncol(x) == length(y)){
        cd <- cbind(colData(x), y)
        colData(x) <- cd
    } else if(nrow(x) == nrow(y)){
        mc <- cbind(mcols(x), y)
        mcols(x) <- mc
    } else{
        stop("Neither row nor column dimensions match")
    }
    mc
})


## obsolete
getMapping <- function(x, map){
    stop('Obsolete')
}





if(getRversion() >= "2.15.1") globalVariables(c(
                                  'wellKey',
                                  'primerid', 
                                  'variable')) #setAs('SingleCellAssay', 'data.table')

setAs('SingleCellAssay', 'data.table', function(from){
    melt.SingleCellAssay(from)
})

setMethod('exprs', 'SingleCellAssay', function(object) t(assay(object)))

setReplaceMethod('exprs', 'SingleCellAssay', function(object, value){
    assay(object, 1) <- t(value)
    object
})
