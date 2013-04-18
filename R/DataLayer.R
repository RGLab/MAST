setGeneric('conform', function(dl, other) standardGeneric('conform'))
setGeneric('nlayer', function(x) standardGeneric('nlayer'))
setGeneric('layer', function(x) standardGeneric('layer'))
setGeneric('layer<-', function(x, value) standardGeneric('layer<-'))
setGeneric('addlayer', function(x, name) standardGeneric('addlayer'))
#setGeneric('dellayer', function(x, i) standardGeneric('dellayer'))



setMethod('addlayer', signature(x='DataLayer', name='character'), function(x, name){
  newLayer <- array(NA, dim=c(nrow(x), ncol(x), 1), dimnames=c(dimnames(x)[-3], name))
  x@.Data <- abind(x@.Data, newLayer)
  x
})

##' Get or set a matrix of measurement values in a \code{SingleCellAssay}
##'
##' Return or set a matrix of the measurement: cells by primerids
##' @title exprs
##' @name exprs
##' @param object DataLayer
##' @return numeric matrix
##' @docType methods
##' @rdname exprs-methods
##' @aliases exprs,DataLayer-method
##' @importMethodsFrom Biobase exprs
##' @return a \code{matrix} of measurement values with wells on the rows and features on the columns of the default layer
##' @export
setMethod("exprs",signature(object="DataLayer"),function(object){
  o <- object@.Data[,,layer(object), drop=FALSE]
  dn <- dimnames(o)
  dim(o) <- dim(o)[-3]
  dimnames(o) <- dn[-3]
  o
})


setMethod('initialize', 'DataLayer',
          function(.Object, ...){
            ##message('init DataLayer') #DEBUG
            .Object <- getMethod('initialize', 'ANY')(.Object, ...)
            .Object
          })

##' @importMethodsFrom Biobase "pData<-"
##' @importMethodsFrom Biobase pData
##' @importMethodsFrom Biobase "exprs<-"
##' @rdname exprs-methods
##' @name exprs
##' @exportMethod "exprs<-"
##' @docType methods
##' @aliases exprs<-,DataLayer,ANY-method

setReplaceMethod('exprs', c('DataLayer', 'ANY'),
                 function(object, value){
                   #if(!conform(object, value)) stop('Replacement must be same dimension as target')
                   object[[,]] <- value
                   object@valid <- FALSE
                   object
                 })


##' Do two objects conform in dimension and type?
##'
##' Returns false if:
##' other does not inherit from Matrix
##' Otherwise returns a bitmask, which is the sum of:
##' 1: objects have same number of rows
##' 2: objects have same number of columns
##' 4: objects have same number of layers (when compareLayers is TRUE)
##' @title conform
##' @param dl DataLayer
##' @param other Another object
##' @param compareLayers Should the number of layers be tested?
##' @return bitmask containing number of dimensions that agree
##' @author andrew
setMethod('conform', c('DataLayer', 'ANY'),
          function(dl, other){
            if(!(inherits(other, 'matrix') || inherits(other, 'array'))) return(FALSE)
            (nrow(dl) == nrow(other))*1 + (ncol(dl) == ncol(other))*2 +
              ifelse(inherits(other, 'DataLayer'),(nlayer(dl)==nlayer(other))*4, 0)
          })

setMethod('ncol', 'DataLayer',
          function(x){
            #if(length(x)==0) return(0)
            ncol(x@.Data[,,x@layer,drop=FALSE])
          })

setMethod('nrow', 'DataLayer',
          function(x){
            #if(length(x)==0) return(0)
            nrow(x@.Data[,,x@layer,drop=FALSE])
          })


setMethod('nlayer', 'DataLayer',
          function(x){
            dim(x@.Data)[3]
          })

setMethod('[', 'DataLayer', function(x, i, j, ..., drop=FALSE){
  if(!missing(i) && is.matrix(i)) stop('Only rectangular selections permitted')
 out <- .subsetHelper(x, i, j, ..., drop=drop)
  new('DataLayer', out, valid=x@valid, layer=x@layer)
})

.subsetHelper <- function(x, i, j, ..., drop=FALSE){
 vargs <- list(...)
  if(length(vargs)>0 || (!missing(i) && is.matrix(i) && ncol(i)>2)) stop('incorrect number of dimensions')
  #if(length(x)==0) return(numeric(0))
  out <- x@.Data[i,j,,drop=FALSE]
}


##' @name [[
##' @title subset methods
##' @details Returns the matrix representation (of the current layer) of the \code{DataLayer}
##' @return \code{matrix}
##' @exportMethod '[['
##' @aliases [[,DataLayer,ANY-method
setMethod('[[', 'DataLayer', function(x, i, j, ..., drop=FALSE){
    if(!missing(i) && is.matrix(i) && ncol(i)==2){                        #matrix indexing
    ## i <- cbind(i[rep(1:nrow(i), times=nlayer(x)),], #make nlayer copies of i, appending 1...nlayer onto it
    ##            rep(seq_len(nlayer(x)), each=nrow(i)))
    i <- cbind(i, layer(x))
    out <- x@.Data[i,drop=FALSE]
  } else{
  out <- .subsetHelper(x, i, j, ..., drop=drop)
  out <- out[,,layer(x), drop=drop]
}
  if(!drop) dim(out) <- dim(out)[-length(dim(out))] #kill layer dimension
  out
})



setMethod('[[<-', 'DataLayer', function(x, i, j, ..., value){
   vargs <- list(...)
  if(length(vargs)>0 || (!missing(i) && is.matrix(i) && ncol(i)>2)) stop('incorrect number of dimensions')
  if(!missing(i) && is.matrix(i) && ncol(i)==2){                        #matrix indexing
    i <- cbind(i, layer(x))
    x@.Data[i] <- value
  } else{
   x@.Data[i,j,layer(x)] <- value
}
  x@valid <- FALSE
  x
})





##' show methods
##' @exportMethod show
##' @aliases show,DataLayer-method
##' @rdname show-methods
##'
setMethod("show","DataLayer",function(object){
  cat(class(object), '\n', nlayer(object), " Layers; ", nrow(object), " wells; ", ncol(object), " features\n")
  invisible(NULL)
})

setMethod('get', c('DataLayer', 'ANY'), function(x, pos){
  #if(length(x)==0) return(numeric(0))
  x[,,pos]
})

setMethod('layer', c('DataLayer'), function(x){
  x@layer
})

setReplaceMethod('layer', c('DataLayer', 'numeric'), function(x, value){
  if(round(value)!=value) stop('Index must be integer')
  if(value < 1 || value > nlayer(x)) stop('Index out of range')
  x@layer <- value
  x
})

setReplaceMethod('layer', c('DataLayer', 'character'), function(x, value){
  if(length(intersect(value, dimnames(x)[[3]]))!=1) stop('Bad index ', value)
  x@layer <- match(value, dimnames(x)[[3]])
  x
})


##'Combine two SingleCellAssay or derived classes
##'
##' Combines two Single Cell-like objects provided they have the same number of Features and Layers.
##' The union of columns from featureData will be taken
##' The union (padded if necessary with NA) will be taken from cellData.
##' @importMethodsFrom BiocGenerics combine
##' @importFrom abind abind
##' @exportMethod combine
##' @aliases combine,SingleCellAssay,SingleCellAssay-method
##' @docType methods
##' @rdname combine-methods
setMethod('combine', signature(x='DataLayer', y='DataLayer'), function(x, y, ...) {
   if(!conform(x, y)>=6){
     stop('Objects must have same number of features and layers; x has dim ', paste(dim(x), ' '), '; y has dim ', paste(dim(y),' '))
   }
   .Data <- abind(x@.Data, y@.Data, along=1)
   proto <- new(class(x), .Data=.Data, valid=x@valid&y@valid)
     for( sl in setdiff(names(getSlots(class(x))), c('valid', '.Data'))){
        slot(proto, sl) <- slot(x, sl)
    }
   proto
 })
