setAs("SCASet","list",function(from)from@set)

##' @importMethodsFrom BiocGenerics lapply
##' @export lapply
##' @rdname lapply-methods
##' @docType methods
##' @aliases lapply,SCASet-method
##' @title lapply
##' @name lapply
setMethod("lapply",c("SCASet"),function(X,FUN,...){
  FUN<-match.fun(FUN)
  X<-as(X,"list")
  .Internal(lapply(X,FUN))
})



##' Index into an SCASet object
##'
##' Return a \code{SingleCellAssay} at index/indices i
##' @export
##' @docType methods
##' @name [[
##' @param x object to be subscripted
##' @param i index
##' @param ... Ignored
##' @return subscripted SingleCellAssay or derived class
##' @details \code{signature(x="SCASet", i="ANY")}: \code{x[[i]]}, where \code{i} is length-1 integer or character matching sampleNames.  Returns the SingleCellAssay at position or with sampleName \code{i}.
##' @rdname doubleAngleBracket-methods
##' @aliases [[,SCASet,ANY-method
##' @keywords transform
##' @export
try({
setMethod("[[",signature("SCASet","ANY"),function(x,i,...){
  if(length(i)!=1)
    stop("subscript out of bounds (index must have length 1)")
  sca<-x@set[[if(is.numeric(i))
              sampleNames(x)[[i]]
  else i]]
  return(sca)
})})

##' Subset an SCASet to a smaller SCASet
##' @export
##' @name [
##' @docType methods
##' @aliases [,SCASet,ANY-method
##' @keywords transform
##' @rdname singleAngleBracket-methods
##' @param x \code{SCASet} to be subscripted
##' @param i index
##' @return An \code{SCASet}
##' @details \code{signature(x="SCASet",i="ANY")}: \code{x[i]}, where \code{i} is the vector of integers or characters matching sampleNames. Returns an SCASet of \code{length(i)} of elements at positions in \code{i}.
try({
  setMethod("[",signature("SCASet","ANY"),function(x,i,...drop=FALSE){
    ret<-x
    ret@set <- x@set[if(is.numeric(i))
                       sampleNames(x)[i]] 
    return(ret)
  })
})

##'get the sample names in an SCASet
##' @name sampleNames
##' @export
##' @aliases sampleNames
##' @aliases sampleNames,SCASet-method
##' @docType methods
##' @keywords transform
##' @param object The \code{SCASet} object
##' @importMethodsFrom Biobase sampleNames
##' @importMethodsFrom Biobase "sampleNames<-"
##' @return a \code{character} vector of the names of the samples in the \code{SCASet}
setMethod("sampleNames","SCASet",function(object){
  unlist(lapply(object@set,function(x)x@id),use.names=FALSE)
})


##' @rdname show-methods
##' @exportMethod show
##' @aliases show,SCASet-method
setMethod("show","SCASet",function(object){
  cat("SCASet of size ",length(object@set),"\n")
  cat("Samples ",paste(unlist(lapply(object,function(x)x@id),use.names=T),collapse=", "),"\n")
  invisible(NULL)
})

##' Get the length of an SCASet
##' @aliases length,SCASet-method
##' @exportMethod length
##' @rdname length-methods
##' @name length
##' @title length of an SCASet
try({setMethod('length', 'SCASet', function(x){
  length(x@set)
})}, silent=TRUE)
