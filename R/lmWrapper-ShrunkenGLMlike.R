##' @include AllClasses.R
##' @include AllGenerics.R
setMethod('initialize', 'ShrunkenGLMlike', function(.Object, ...){
    callNextMethod()
})


setMethod('fit', signature=c(object='ShrunkenGLMlike', response='missing'), function(object, response, silent=TRUE, ...){
    object <- callNextMethod()
    if(object@fitted['C']){
    thisdispersion <- object@fitC$dispersionMLE
    df.total <- object@fitC$df.null+1
    df.residual <- object@fitC$df.residual
    object@fitC$dispersionMLEnoshrink <- thisdispersion
    object@fitC$dispersionMLE <- (thisdispersion*df.total + object@priorVar*object@priorDOF)/(df.total+object@priorDOF)
    object@fitC$dispersion <- (object@fitC$dispersion*df.residual+object@priorVar*object@priorDOF)/(df.residual+object@priorDOF)
}
    object
    
})
