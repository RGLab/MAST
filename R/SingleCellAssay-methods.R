##' Methods for analysing single cell assay data
##'
##' This packages provides data structures and functions for statistical analysis of single-cell assay data such as Fluidigm single cell gene expression assays.
##'
##' @name MAST-package
##' @aliases MAST-package
##' @docType package
##' @title "Tools for Single Cell Assay Analysis"
##' @keywords package
##' @rdname MAST-package
NULL


##' Construct a SummarizedExperiment from a matrix or array of expression
##'
##' If the gene expression measurements are already in a rectangular form,
##' then this function allows an easy way to construct a SummarizedExperiment object while
##' still doing some sanity checking of inputs.
##' @param exprsArray matrix or array, columns are cells, rows are genes
##' @param cData cellData data.frame, AnnotatedDataFrame or DataFrame
##' @param fData featureData data.frame, AnnotatedDataFrame or DataFrame
##' @return an object of class \code{class}
##' @export
##' @examples
##' ncells <- 10
##' ngenes <- 5
##' fData <- data.frame(primerid=LETTERS[1:ngenes])
##' cData <- data.frame(wellKey=seq_len(ncells))
##' mat <- matrix(rnorm(ncells*ngenes), nrow=ngenes)
##' sca <- FromMatrix(mat, cData, fData)
##' stopifnot(inherits(sca, 'SummarizedExperiment0'))
##' ##If there are mandatory keywords expected by a class, you'll have to manually set them yourself
##' cData$ncells <- 1
##' fd <- FromMatrix('FluidigmAssay', mat, cData, fData)
##' stopifnot(inherits(fd, 'FluidigmAssay'))
FromMatrix <- function(exprsArray, cData, fData){
    can <- checkArrayNames(exprsArray, cData, fData)
    nslice <- dim(can$exprsArray)[3]
    assays <- vector('list', length=nslice)
    for(i in seq_len(nslice)){
        assays[[i]] <- can$exprsArray[,,i,drop=FALSE]
        dim(assays[[i]]) <- dim(assays[[i]])[-3] # only drop last index
        dimnames(assays[[i]]) <- dimnames(can$exprsArray)[-3]
    }
    obj <- SummarizedExperiment(assays=assays, colData=can$cData)
    mcols(obj) <- can$fData
    obj
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

setMethod('fData', 'SummarizedExperiment0', function(object){
    warning('Deprecated: use mcols')
    mcols(object)
})



##' Melt a rectangular array
##'
##' Return a 'molten' (flat) representation of a rectangular array
##'
##' @param data A rectangular array, with attributes attached to its rows and
##' columns
##' @param ... ignored
##' @return A \code{data.frame} typically, with the cartesian product of the
##' row and column attributes and the values from the rectangular array
##' 
##' @export
melt.SummarizedExperiment0<-function(data,...,na.rm=FALSE, value.name='value'){
    featdata <- as.data.table(mcols(data))
    celldata <- as.data.table(colData(data))
    m <- cbind( featdata[rep(seq_len(nrow(featdata)), nrow(celldata)),,drop=FALSE],
          celldata[rep(seq_len(nrow(celldata)), each=nrow(featdata)),,drop=FALSE],
          value=as.vector(assay(data)))      
    
  setnames(m, 'value', value.name) 
  return(m)
}

##'@export
melt.SingleCellAssay <- melt.SummarizedExperiment0

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
                 'wellKey')) #fixdf

## might have bad complexity, but could construct one at time, then glue cheaply
## Not too bad except for deduplication.. will use data.table
## Output is sorted by primerid then wellKey
##' @import data.table
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
      setnames(df,cn_df,make.unique(c(nm, colnames(df)))[-1])
    warning("renaming column ", nm)
    }
      if(!all(bothMap[[nm]] %in% names(df))){
        stop('could not find column named ', bothMap[[nm]], ' in dataframe')
      }
    rn <- nm
    names(rn) <- as.character(bothMap[[nm]])
    if(keep.names){
        df[,eval(rn):=get(names(rn))]
    } else{ setnames(df, names(rn),nm)}
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

##' Construct a SummarizedExperiment from a `flat` (melted) data.frame/data.table
##'
##' @param dataframe A 'flattened' \code{data.frame} or \code{data.table} containing columns giving cell and feature identifiers and  a measurement column
##' @param idvars character vector naming columns that uniquely identify a cell
##' @param primerid character vector of length 1 that names the column that identifies what feature (i.e. gene) was measured
##' @param measurement character vector of length 1 that names the column containing the measurement 
##' @param id An identifier (eg, experiment name) for the resulting object
##' @param cellvars Character vector naming columns containing additional cellular metadata
##' @param featurevars Character vector naming columns containing additional feature metadata
##' @param phenovars Character vector naming columns containing additional phenotype metadata
##' @param ... additional arguments are ignored
##' @export
##' @aliases SingleCellAssay
##' @examples
##' ## See FluidigmAssay for examples
##' ##' @examples
##' data(vbeta)
##' colnames(vbeta)
##' vbeta <- computeEtFromCt(vbeta)
##' vbeta.fa <- FluidigmAssay(vbeta, idvars=c("Subject.ID", "Chip.Number", "Well"),
##' primerid='Gene', measurement='Et', ncells='Number.of.Cells',
##' geneid="Gene",cellvars=c('Number.of.Cells', 'Population'),
##' phenovars=c('Stim.Condition','Time'), id='vbeta all')
##' show(vbeta.fa)
##' nrow(vbeta.fa)
##' ncol(vbeta.fa)
##' head(fData(vbeta.fa)$primerid)
##' table(cData(vbeta.fa)$Subject.ID)
##' vbeta.sub <- subset(vbeta.fa, Subject.ID=='Sub01')
##' show(vbeta.sub)
##' @return SummarizedExperiment object
FromFlatDF<-function(dataframe,idvars,primerid,measurement,id=numeric(0), cellvars=NULL, featurevars=NULL, phenovars=NULL, ...){
    if(missing(dataframe) || missing(idvars) || missing(primerid) || missing(measurement)){
        stop("Must supply all of 'idvars', 'primerid', 'measurement'  and 'dataframe'")
    }
    if(!is(dataframe,"data.table")){
        dataframe<-data.table(dataframe)
    }
    setkeyv(dataframe, colnames(dataframe))
    ## fixdf: make primerid unique, generate idvar column, rename columns according to cmap and fmap, complete df
    ## keeping support for "mappings" in case we change our mind
    fmap <- 'primerid'
    cmap <- 'wellKey'
    fixed <- fixdf(dataframe, idvars, primerid, measurement,
                   cmap=cmap, fmap=fmap)
    dl <- array(as.matrix(fixed$df[,measurement, with=FALSE]),
                dim=c(length(fixed$rn), length(fixed$cn), length(measurement)),
                dimnames=list(wellKey=fixed$rn, primerid=fixed$cn, layer=measurement))
    dl <- aperm(dl, c(2, 1, 3))
    cellvars <- unique(c(cellvars, cmap, idvars, phenovars))
    featurevars <- unique(c(fmap, primerid, featurevars))
    check.vars(cellvars, featurevars, phenovars, fixed$df, length(fixed$cn), length(fixed$rn))
    cell.adf<-DataFrame(uniqueModNA(fixed$df[,cellvars,with=FALSE], 'wellKey') )
    row.names(cell.adf) <- unique(fixed$df$wellKey)
    ##pheno.adf <- new('AnnotatedDataFrame')
    ##need a phenokey into the melted data frame for this to make sense
    f.adf <- DataFrame(uniqueModNA(fixed$df[,featurevars, with=FALSE], 'primerid'))
    row.names(f.adf) <- unique(fixed$df$primerid)
    FromMatrix(dl, cell.adf, f.adf)
}

FluidigmAssay <- SingleCellAssay <- function(...){
    warning("Deprecated: Use FromFlatDF")
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


setMethod('getwellKey', 'SingleCellAssay', function(sc) {cData(sc)$wellKey})



##' @describeIn cData
##' @export
setMethod('cData', 'SummarizedExperiment0', function(sc){
    warning('Deprecated: use colData')
    colData(sc)
})

setMethod('subset', 'SummarizedExperiment0', function(x, ...){
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


##' @describeIn cData
##' @export
setReplaceMethod("cData", "SingleCellAssay", function(sc, value) {
    warning('Deprecated: use colData<-')
    colData(sc) <- value
    sc
})


##' Split into SimpleList
##'
##' Splits a \code{SummarizedExperiment0} into a \code{SimpleList} by a factor (or something coercible into a factor) or a character giving a column of \code{colData(x)}
##' @param x SingleCellAssay
##' @param f length-1 character or factor of length nrow(x)
##' @return SimpleList
##' @examples
##' data(vbetaFA)
##' split(vbetaFA, 'ncells')
##' fa <- as.factor(colData(vbetaFA)$ncells)
##' split(vbetaFA, fa)
##' @export
setMethod('split', signature(x='SummarizedExperiment0', f='character'), 
          function(x, f, drop=FALSE, ...){
              ## Split a SingleCellAssay by criteria
###f must be a character naming a cData variable
              if(length(f) != ncol(x)){
                  f <- lapply(colData(x)[,f, drop=FALSE], as.factor)
              } else{
                  f <- as.factor(f)
              }
              callNextMethod(x, f, drop, ...)
})


## obsolete
getMapping <- function(x, map){
  stop('Obsolete')
}





if(getRversion() >= "2.15.1") globalVariables(c(
                  'wellKey',
                 'primerid', 
                  'variable')) #setAs('SingleCellAssay', 'data.table')

setAs('SummarizedExperiment0', 'data.table', function(from){
    melt.SummarizedExperiment0(from)
})

setMethod('exprs', 'SummarizedExperiment0', function(object) t(assay(object)))

setReplaceMethod('exprs', 'SummarizedExperiment0', function(object, value){
    assay(object, 1) <- t(value)
    object
})
