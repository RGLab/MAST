##' Methods for analysing single cell assay data
##'
##' This packages provides data structures and functions for statistical analysis of single-cell assay data such as Fluidigm single cell gene expression assays.
##'
##' @name SingleCellAssay-package
##' @aliases SingleCellAssay-package
##' @docType package
##' @title "Tools for Single Cell Assay Analysis"
##' @keywords package
##' @rdname SingleCellAssay-package
##' @seealso \code{\link{SCASet}}
NULL

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
##' 
##' @rdname melt
##' @title melt
##' @aliases melt
##' @keywords transformation
##' @importFrom reshape melt
##' @importFrom reshape melt.default
##' @importFrom reshape melt.array
##' @importFrom reshape melt.cast_df
##' @importFrom reshape melt.list
##' @importFrom reshape melt.matrix
##' @importFrom reshape melt.cast_matrix
##' @importFrom reshape melt.data.frame
##' @importFrom reshape melt.table
##' @export
## melt.SingleCellAssay<-function(data,...){
##   m <- melt.data.frame(cbind(cData(data), exprs(data)), id.vars=names(cData(data)), variable_name='primerid')
##   m <- merge(m, fData(data), by='primerid')
##   if(data@keep.names) return(rename(m, c('value'=dimnames(data)[[3]][layer(data)])))
##   return(m)
## }
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
##' @importFrom reshape expand.grid.df
##' @importFrom reshape rename
## unnamed arguments get passed along to callNextMethod
## which eventually just sets the slots
setMethod('initialize', 'SingleCellAssay',
          function(.Object, dataframe, idvars, primerid, measurement, exprsMatrix, cellvars=NULL, featurevars=NULL, phenovars=NULL, sort=TRUE, ...){
            ##message(class(.Object), ' calling SingleCellAssay Initialize')  #DEBUG
            .Object <- callNextMethod()
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

sort.SingleCellAssay <- function(x, decreasing=FALSE, ...){
  if(nrow(cData(x))>0){
    rowOrd <- order(cData(x)$wellKey)
    x@.Data <- as(x, 'DataLayer')[rowOrd,]
    x@cellData <- x@cellData[rowOrd,]
  }
  if(nrow(fData(x))>0){
     colOrd <- order(fData(x)$primerid)
    x@.Data <- as(x, 'DataLayer')[,colOrd]
    x@featureData <- x@featureData[colOrd,]
  }
  stopifnot(validObject(x))
  x
}

setMethod('sort', signature=c(x='SingleCellAssay'), sort.SingleCellAssay)

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

setGeneric("melt",function(data,...){
standardGeneric("melt")
#  UseMethod(generic="melt",data)
 },useAsDefault=reshape::melt)



setMethod('getwellKey', 'SingleCellAssay', function(sc) {cData(sc)$wellKey})



##' @describeIn cData
##' @export
setMethod('cData', 'SingleCellAssay', function(sc)  pData(sc@cellData))



##' @describeIn cData
##' @export
setReplaceMethod("cData", "SingleCellAssay", function(sc, value) {
  if (is.data.frame(value)) {
    value <- as(value, "AnnotatedDataFrame")
  }
  if (!is(value, "AnnotatedDataFrame")) {
    stop("'value' must be either a data.frame or an AnnotatedDataFrame")
  }
  mandatory <- c('wellKey', names(sc@cmap)) #Must contain the expected fields for the class
  missing <- setdiff(mandatory, varLabels(value)) 
  if(length(missing)>0) stop('cellData is missing mandatory field ', paste(missing, collapse=','))
  inputOrder <- match(getwellKey(sc), value$wellKey)
  if(any(is.na(inputOrder))) stop('cellData is missing some wellkeys')
  if(any(inputOrder != seq_along(value$wellKey))) warning("sorting cellData by wellKey")
  value <- value[inputOrder,] #sort
  sc@cellData <- value
  return(sc)
})


##' @export cellData
setMethod('cellData', 'SingleCellAssay', function(sc) sc@cellData)


##' @rdname fData-methods
##' @aliases fData,SingleCellAssay-method
##' @export fData
setMethod('fData', 'SingleCellAssay', function(object) pData(object@featureData))

##' @rdname featureData-methods
##' @aliases featureData,SingleCellAssay-method
##' @exportMethod featureData
setMethod('featureData', 'SingleCellAssay', function(object)  object@featureData)

##' @rdname melt
##' @details \code{signature(data="SingleCellAssay")}: return a \code{data.frame}, which contains a melted version of \code{data}.
##' @aliases melt,SingleCellAssay-method
##' @exportMethod melt
setMethod("melt","SingleCellAssay",melt.SingleCellAssay )

.scaSubset <- function(x, i, j, ..., drop=FALSE){
  if(missing(i)){
    i<-1:nrow(x)
  }
  if(any(is.na(i))) stop("NAs not permitted in 'i' index")
  if(is.factor(i)) stop("Factors not permitted in 'i' index")
  
  if(is(i,"character")){
    wk<-getwellKey(x)
    if(length(setdiff(i, wk))>0) stop('wellKeys \n', paste(setdiff(i, wk), sep=','), '\n not found!')
    i <- match(i, wk)
  }
  
  if(!missing(j)){
    if(any(is.na(j))) stop("NAs not permitted in 'j' index")
    if(is.factor(j)) stop("Factors not permitted in 'j' index")
    pk<-fData(x)$primerid
    
    if(is(j,"character")){
      J <- match(j, pk)
      if(!(all(j%in%pk))){
        stop("feature names \n",paste(j[!j%in%as.matrix(pk)],collapse=" "), "\n not found!");
      }
      j<-J
    }
    newfdf <- featureData(x)[j,] 
  }else {                               #j missing
    j <- TRUE
  }

  
  newcdf <- cellData(x)[i,]
  if(!exists("newfdf")){
    newfdf<-x@featureData
  }
  ## Not callNextMethod, because we want to dispatch [ not [[
  ## And we don't want to wrap this method 
  .Data <- selectMethod('[', signature='DataLayer')(x, i, j, ..., drop=drop)
  x@featureData <- newfdf
  x@cellData <- newcdf
  x@.Data <- .Data
  x
}


setMethod('[[', signature(x="SingleCellAssay"), .scaSubset)

##' Subset a SingleCellAssay
##' @details \code{signature(x="SingleCellAssay", i="ANY")}: \code{x[i]}, where \code{i} is a logical, integer, or character vector, recycled as necessary to match \code{nrow(x)}. Optional \code{x[[i,j]]} where j is a logical, integer or character vector selecting the features based on ``primerid'' which is unique, while ``geneid'' or gene name is not necessarily unique.
##'
##' @note
##' x[[i,j]] functions similarly, but this behavior may change in future releases.
##' @param x SingleCellAssay
##' @param i logical, integer or character (naming wellKeys)
##' @param j logical, integer or character (naming primerid)
##' @param ... ignored
##' @param drop ignored
##' @return SingleCellAssay, suitably subsetted
##' @aliases [,SingleCellAssay,ANY-method
##' @aliases [,SingleCellAssay-method
##' @aliases [[,SingleCellAssay-method
##' @exportMethod [
setMethod("[", signature(x="SingleCellAssay"), .scaSubset)
  


##' Subset a SingleCellAssay by data in cellData
##'
##' @param x SingleCellAssay
##' @param thesubset expression, which when evaluated in cellData environment which returns a logical
##' @rdname subset
##' @export
##' @examples
##' data(vbetaFA)
##' subset(vbetaFA, ncells==1)
#\code{signature(x='SingleCellAssay', thesubset='ANY')}: Return a new SingleCellAssay consisting of cells in which thesubset is TRUE
setMethod('subset', 'SingleCellAssay', function(x, thesubset, ...){
  e <- substitute(thesubset)
  asBool <- try(eval(e, cData(x), parent.frame(n=2)), silent=TRUE)
  if(is(asBool, 'try-error')) stop(paste('Variable in subset not found:', strsplit(asBool, ':')[[1]][2]))
  #this is a special case of "subset", not of the "[[" method, so..
  if(length(asBool)==1){
    if(asBool==TRUE){
      x 
    }else{
      x[[asBool]]
    }
  }else{
    x[[asBool]]
  }
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

##'Combines a SCASet into a unified SingleCellAssay
##'
##' No error checking is currently done to insure that objects in the SCASet conform with each other,
##' so mysterious errors may result if they do not.
##' @importMethodsFrom BiocGenerics combine
##' @export
##' @aliases combine,SCASet,missing-method
setMethod('combine', signature=c(x='SCASet', y='missing'), function(x, y, ...){
    
    skeleton <- x[[1]]
    if(length(x) == 1) return(skeleton)
    for(i in seq(from=2, to=length(x))){
        skeleton <- combine(skeleton, x[[i]])
    }
    return(skeleton)
})

##'Combine two SingleCellAssay or derived classes
##'
##' Combines two Single Cell-like objects provided they have the same number of Features and Layers.
##' The union of columns from featureData will be taken
##' The union (padded if necessary with NA) will be taken from cellData.
##' @importMethodsFrom BiocGenerics combine
##' @import abind
##' @export
##' @aliases combine,DataLayer,Datalayer-method
##' @note
##' We might also wish to combine features along each cell/row but a use case for this
##' hasn't arrived yet.
##' @aliases combine,SingleCellAssay,SingleCellAssay-method
setMethod('combine', signature(x='SingleCellAssay', y='SingleCellAssay'), function(x, y, ...) {
  proto <- callNextMethod()
  cellData <- combine(cellData(x), cellData(y))
  ## Combine wants to do a rbind-like operation with ADFs
  ## But we want a Cbind
  nx <- names(fData(x))
  ny <- names(fData(y))
  featureData <- featureData(x)
  pData(featureData) <- cbind(fData(x), fData(y)[,setdiff(ny, nx)])
  proto@cellData <- cellData
  proto@featureData <- featureData
  proto
})

## obsolete
getMapping <- function(x, map){
  return(list(map))
}

##' @describeIn show
setMethod('show', 'SingleCellAssay', function(object){
callNextMethod()
cat(' id: ', object@id, '\n')
})

##' Combine a SingleCellAssay and a vector, data.frame or AnnotatedDataFrame
##'
##' When \code{x} is a SingleCellAssay and \code{y} is a  vector, data.frame or AnnotatedDataFrame,
##' an attempt is made bind it to the columns of the cellData or featureData.
##' The behavior depends on the number of rows/length of \code{y}.
##'If \code{nrow(y) == nrow(x)}, then the cellData is used.
##' If \code{nrow(y) == ncol(x)}, then the cellData is used.
##' It is an error if neither mathces.
##' @param x SingleCellAssay
##' @param y vector, data.frame or AnnotatedDataFrame
##' @param ... ignored
##' @aliases combine,SingleCellAssay,ANY-method
##' @aliases combine,SingleCellAssay,data.frame-method
##' @aliases combine,SingleCellAssay,AnnotatedDataFrame-method
##' @return SingleCellAssay
setMethod('combine', signature=c(x='SingleCellAssay', y='ANY'), function(x, y, ...){
  adf <- new('AnnotatedDataFrame')
  df <- data.frame(y)
  names(df) <- deparse(substitute(y))
  pData(adf) <- df
  selectMethod('combine', c(x=class(x), y='AnnotatedDataFrame'))(x, adf)
})

setMethod('combine', signature=c(x='SingleCellAssay', y='data.frame'), function(x, y, ...){
  adf <- new('AnnotatedDataFrame')
  pData(adf) <- y
  selectMethod('combine', c(x=class(x), y='AnnotatedDataFrame'))(x, adf)
})

setMethod('combine', signature=c(x='SingleCellAssay', y='AnnotatedDataFrame'), function(x, y, ...){
    if(nrow(x) == ncol(x)) stop("x has same number of rows and columns, must explicitly specify 'along'")
    if(nrow(y) == nrow(x)) along <- 'cellData'
    else if(nrow(y) == ncol(x)) along <- 'featureData'
    else stop('Dimension mismatch between y and x')
    
  if(length(intersect(along , c('cellData', 'featureData')))!=1) stop("If specified, along must be either 'cellData', or 'featureData'")
  newdata <- slot(x, along)
  pData(newdata) <- cbind(pData(newdata), pData(y))
  slot(x, along) <- newdata
  x
})

setAs('ExpressionSet', 'SingleCellAssay', function(from){
    ## just a transposed version
    ex <- t(exprs(from))
    dn <- dimnames(ex)
    names(dn) <- c('wellKey', 'primerid')
    dim(ex) <- c(dim(ex), 1)
    pd <- phenoData(from)
    pData(pd)[,'wellKey'] <- sampleNames(pd)
    fd <- featureData(from)
    fd$primerid <- sampleNames(fd)
    dimnames(ex) <- c(dn, layer='ExpressionSet')
    DL <- new('DataLayer', .Data=ex)
    new('SingleCellAssay', .Data=DL, featureData=fd, cellData=pd, sort=FALSE)
})

melt.data.table <- function(dt, id.var){
    ## Well...that's unfortunate
    dt[, list(variable = names(.SD), value = unlist(.SD, use.names = FALSE)), keyby =eval(bquote(.(id.var)))]
}

if(getRversion() >= "2.15.1") globalVariables(c(
                  'wellKey',
                 'primerid', 
                  'variable')) #setAs('SingleCellAssay', 'data.table')

setAs('SingleCellAssay', 'data.table', function(from){
    ex <- data.table(wellKey=getwellKey(from), exprs(from))
    fd <- setkey(data.table(fData(from)), primerid)
    cd <- setkey(data.table(cData(from)), wellKey)
    M <- melt.data.table(ex, 'wellKey')[,primerid:=variable][,variable:=NULL]
    setkey(M, primerid, wellKey)
    M <- merge(merge(M, cd), fd, by='primerid')
    stopifnot(nrow(M) == ncol(from)*nrow(from))
    M
})
