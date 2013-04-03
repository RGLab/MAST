setGeneric('conform', function(dl, other) standardGeneric('conform'))
setGeneric('nlayer', function(x) standardGeneric('nlayer'))
setGeneric('layer', function(x) standardGeneric('layer'))
setGeneric('layer<-', function(x, value) standardGeneric('layer<-'))


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
##' @export exprs
setMethod("exprs",signature(object="DataLayer"),function(object){
  o <- object@.Data[,,layer(object), drop=FALSE]
  dn <- dimnames(o)
  dim(o) <- dim(o)[-3]
  dimnames(o) <- dn[-3]
  o
})


setMethod('initialize', 'DataLayer',
          function(.Object, ...){
            .Object <- getMethod('initialize', 'ANY')(.Object, ...)
            if(length(.Object@.Data)==1)
              dim(.Object@.Data) <- c(1, 1, 1)
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
                   if(!conform(object, value)) stop('Replacement must be same dimension as target')
                   object[,] <- value
                   object@valid <- FALSE
                   object
                 })

setMethod('conform', c('DataLayer', 'ANY'),
          function(dl, other){
            inherits(other, 'matrix') && nrow(dl) == nrow(other) && ncol(dl) == ncol(other)
          })

setMethod('ncol', 'DataLayer',
          function(x){
            if(length(x)==0) return(0)
            ncol(x@.Data[,,x@layer,drop=FALSE])
          })

setMethod('nrow', 'DataLayer',
          function(x){
            if(length(x)==0) return(0)
            nrow(x@.Data[,,x@layer,drop=FALSE])
          })


setMethod('nlayer', 'DataLayer',
          function(x){
            dim(x@.Data)[3]
          })

setMethod('[', 'DataLayer', function(x, i, j, ..., drop=TRUE){
  vargs <- list(...)
  if(length(vargs)>0 || (is.matrix(i) && ncol(i)>2)) stop('incorrect number of dimensions')
  if(length(x)==0) return(numeric(0))
  
  if(is.matrix(i) && ncol(i)==2){                        #matrix indexing
    i <- cbind(i, layer(x))
    out <- x@.Data[i,drop=drop]
  } else{
  out <- x@.Data[i,j,layer(x),drop=drop]
}
  if(!drop) dim(out) <- dim(out)[-length(dim(out))] #kill layer dimension
  out
})


setMethod('[<-', 'DataLayer', function(x, i, j, ..., value){
   vargs <- list(...)
  if(length(vargs)>0 || (is.matrix(i) && ncol(i)>2)) stop('incorrect number of dimensions')
  if(is.matrix(i) && ncol(i)==2){                        #matrix indexing
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
  if(length(x)==0) return(numeric(0))
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
