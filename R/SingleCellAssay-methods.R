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
    obj <- SummarizedExperiment(assays=can$exprsArray, colData=can$cData)
    mcols(obj) <- can$fData
    obj
}

checkArrayNames <- function(exprsArray, cData, fData){
    if(!is.numeric(exprsArray)) stop('`exprsArray` must be numeric')
    if(length(dim(exprsArray))!=2) stop('`exprsArray` must be matrix or 2-d array')
    dn <- dimnames(exprsArray)
    noDimnames <- is.null(dn) || is.null(dn[[1]]) || is.null(dn[[2]])

    pidDefault <- if(is.null(dn[[1]])) sprintf('p%0*d', ceiling(log10(nrow(exprsArray)+1)), seq_len(nrow(exprsArray))) else dn[[1]]
    wkDefault <- if(is.null(dn[[2]])) sprintf('wk%0*d', ceiling(log10(ncol(exprsArray)+1)), seq_len(ncol(exprsArray))) else dn[[2]]
    
    if(missing(fData)) fData <- S4Vectors::DataFrame(primerid=pidDefault)
    if(missing(cData)) cData <- S4Vectors::DataFrame(wellKey=wkDefault)
    
    
    if(ncol(exprsArray) != nrow(cData)) stop('`cData` must contain as many rows as `exprsArray`')
    if(nrow(exprsArray) != nrow(fData)) stop('`fData` must contain as many columns as `exprsArray`')

    if(!('primerid' %in% names(fData))){
        warning("`fData` has no primerid.  I'll make something up.")
        fData[['primerid']] <- pidDefault
    } else{
        fData[['primerid']] <- as.character(fData$primerid)
    }
    row.names(fData) <- fData$primerid

    if(!('wellKey' %in% names(cData))){
        warning("`cData` has no wellKey.  I'll make something up.")
        cData[['wellKey']] <- wkDefault
    } else{
        cData[['wellKey']] <- as.character(cData$wellKey)
    }
    row.names(cData) <- cData$wellKey
    

    if(noDimnames){
        message('No dimnames in `exprsArray`, assuming `fData` and `cData` are sorted according to `exprsArray`')
        dn <- list(primerid=row.names(fData), wellkey=row.names(cData))  
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
##' @import reshape
##' @export
melt.SingleCellAssay<-function(data,...){
  m <- cbind(
    cData(data)[rep(seq_len(nrow(data)), ncol(data)),,drop=FALSE],
    fData(data)[rep(seq_len(ncol(data)), each=nrow(data)),,drop=FALSE],
    value=as.vector(exprs(data)))
  rn <-  c('value'=dimnames(data)[[3]][layer(data)])
  if(data@keep.names) return(rename(m,rn))
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
                 'wellKey')) #fixdf

## might have bad complexity, but could construct one at time, then glue cheaply
## Not too bad except for deduplication.. will use data.table
##' @import data.table
fixdf <- function(df, idvars, primerid, measurement, cmap, fmap, keep.names){
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

##' @importFrom plyr ddply
## unnamed arguments get passed along to callNextMethod
## which eventually just sets the slots
setMethod('initialize', 'SingleCellAssay',
          function(.Object, dataframe, idvars, primerid, measurement, exprsMatrix, cellvars=NULL, featurevars=NULL, phenovars=NULL, sort=TRUE, ...){
            ##message(class(.Object), ' calling SingleCellAssay Initialize')  #DEBUG
            .Object <- callNextMethod(.Object, ...)
            if(sort) .Object <- sort(.Object)
            if(!missing(dataframe)){              #called using melted dataframe
              ##message('...with dataframe') #DEBUG
              if(missing(idvars) || missing(primerid) || missing(measurement)){
                stop("Must supply all of 'idvars', 'primerid' and 'measurement' if 'dataframe' is passed")
              }
              if(nrow(.Object) > 0 || ncol(.Object)>0) warning('slots will be overwritten when dataframe is provided')
              
              ## fixdf: make primerid unique, generate idvar column, rename columns according to cmap and fmap, complete df
              if(!is(dataframe,"data.table")){
                dataframe<-data.table(dataframe)
              } else{
                  ## worry about things being passed by reference
                  #dataframe <- data.table::copy(dataframe)
              }
              setkeyv(dataframe, colnames(dataframe))
              
              fixed <- fixdf(dataframe, idvars, primerid, measurement, .Object@cmap, .Object@fmap, .Object@keep.names)
              .Object@fmap <- fixed$fmap        #needed if we deduplicated primerid
              .Object@cmap <- fixed$cmap
              dl <- array(as.matrix(fixed$df[,measurement, with=FALSE]),
                          dim=c(length(fixed$rn), length(fixed$cn), length(measurement)),
                          dimnames=list(wellKey=fixed$rn, primerid=fixed$cn, layer=measurement))
              .Object@.Data <- dl
              cellvars <- union(cellvars, c('wellKey', idvars, phenovars, names(.Object@cmap))) #fixme when we support phenovars
              featurevars <- union(c('primerid', primerid, featurevars), names(.Object@fmap))
              if(.Object@keep.names){
                featurevars <- union(featurevars, unlist(.Object@fmap))
                cellvars <- union(cellvars, unlist(.Object@cmap))
              }
              check.vars(cellvars, featurevars, phenovars, fixed$df, length(fixed$cn), length(fixed$rn))
              cell.adf  <- new("AnnotatedDataFrame")
              ## fixed$df should be keyed by cellvars, primerid, ...
              pData(cell.adf)<-uniqueModNA(fixed$df[,cellvars,with=FALSE], 'wellKey') 
              sampleNames(cell.adf) <- unique(fixed$df$wellKey)
              ##pheno.adf <- new('AnnotatedDataFrame')
              ##need a phenokey into the melted data frame for this to make sense
              f.adf <- new('AnnotatedDataFrame')
              pData(f.adf) <- uniqueModNA(fixed$df[,featurevars, with=FALSE], 'primerid')
              sampleNames(f.adf) <- unique(fixed$df$primerid)

              ## currently take care of this in fixdf
              ## if(.Object@keep.names){
              ##   pData(cell.adf)[,unlist(.Object@cmap)] <- pData(cell.adf)[,names(.Object@cmap)]
              ##   pData(f.adf)[,unlist(.Object@fmap)] <- pData(f.adf)[,names(.Object@fmap)]
              ##   pData(f.adf)[,primerid] <- pData(f.adf)[,'primerid']
              ## }
  
              .Object@cellData<-cell.adf
              .Object@featureData <- f.adf
            }
            
            .Object
          })


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

## Get unique rows in data.frame, not counting NAs as distinct for
## columns named in exclude
## Precondition: keyed data.table
uniqueModNA <- function(df, exclude){
    if(!is(df, 'data.table')){
        stop('df should be data.table')
    }
    setkey(df,NULL)                     # so that unique operates on all columns
    w.include <- names(df)
    if(ncol(df)>1){
        w.include <- setdiff(w.include, exclude)
}
    ##u <- unique(df)
    ##anyNa <- apply(is.na(u)[,w.include, drop=FALSE], 1, all)
    ##u[!anyNa,,drop=FALSE]
    u<-unique(df)
    allNa <- apply(is.na(u)[,w.include, drop=FALSE], 1, all)
    u<-u[!allNa,] #either we spend time above or here..
    u
}


setMethod('getwellKey', 'SingleCellAssay', function(sc) {cData(sc)$wellKey})



##' @describeIn cData
##' @export
setMethod('cData', 'SummarizedExperiment0', function(sc){
    warning('Deprecated: use colData')
    colData(sc)
})



##' @describeIn cData
##' @export
setReplaceMethod("cData", "SingleCellAssay", function(sc, value) {
    warning('Deprecated: use colData<-')
  colData(sc) <- value
})


##' Split into SCASet
##'
##' Splits a \code{SingleCellAssay} into a \code{SCASet} by a factor (or something coercible into a factor) or a character giving a column of the melted SingleCellAssay
##' @param x SingleCellAssay
##' @param f length-1 character or factor of length nrow(x)
##' @return SCASet
##' @aliases split,SingleCellAssay,ANY-method
##' @aliases split,DataLayer,ANY-method
##' @examples
##' data(vbetaFA)
##' split(vbetaFA, 'ncells')
##' fa <- as.factor(cData(vbetaFA)$ncells)
##' split(vbetaFA, fa)
##' @export
setMethod('split', signature(x='SingleCellAssay'), 
          function(x, f, drop=FALSE, ...){
  ## Split a SingleCellAssay by criteria
  ###f must be a character naming a cData variable
  if(is(f, 'character')){
    if(length(f) != nrow(x)){
      f <- lapply(cData(x)[,f, drop=FALSE], as.factor)
    } else{
      f <- list(as.factor(f))
    }
  } else if(is(f, 'factor')){
    f <- list(f)
  }
     all.factor <- all(sapply(f, is.factor))
     if(!all.factor) stop('f must be character vector naming columns of x; or factor; or list of factors')
     all.length <- all(sapply(f, length)==nrow(x))
     if(!all.length) stop('each element in f must be length nrow(x)')
  out <- callNextMethod()
  cD <- split.data.frame(x@cellData, f)
  for(i in seq_along(out)){
    out[[i]]@cellData <- cD[[i]]
    out[[i]]@id <- names(out)[i]
  }
  new('SCASet', set=out)
})


## obsolete
getMapping <- function(x, map){
  warning('Obsolete')
  return(list(map))
}

##' @describeIn show
setMethod('show', 'SingleCellAssay', function(object){
callNextMethod()
cat(' id: ', object@id, '\n')
})




if(getRversion() >= "2.15.1") globalVariables(c(
                  'wellKey',
                 'primerid', 
                  'variable')) #setAs('SingleCellAssay', 'data.table')

setAs('SummarizedExperiment0', 'data.table', function(from){
    ex <- data.table(wellKey=getwellKey(from), exprs(from))
    fd <- setkey(data.table(fData(from)), primerid)
    cd <- setkey(data.table(cData(from)), wellKey)
    M <- melt.data.table(ex, 'wellKey')[,primerid:=variable][,variable:=NULL]
    setkey(M, primerid, wellKey)
    M <- merge(merge(M, cd), fd, by='primerid')
    stopifnot(nrow(M) == ncol(from)*nrow(from))
    M
})
