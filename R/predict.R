diag_safeContrastQF = function(contrast, vc){
    ## diag(t(a) %*% X %*% a) = rowSums(t(a) %*% X * a)
    ## we don't want to calculate the off-diagonal entries!
    dp = safeContrastDP(contrast, vc)
    colSums(t(contrast) * dp)
}


#' Return predictions from a ZlmFit object.
#'
#' @param object A \code{ZlmFit}
#' @param newdata The data to predict from. Currently ignored, will use the data in the object.
#' @param modelmatrix The model matrix specifying the linear combination of coefficients.
#' @param ... ignored
#' @return Predictions (on the link scale) and standard errors.
#' @export
#' @examples
#' ##See stat_ell
#' example(stat_ell)
predict.ZlmFit <- function(object,newdata = NULL, modelmatrix=NULL, ...){
    if(is.null(modelmatrix)) modelmatrix = object@LMlike@modelMatrix
    
    C = coef(object,"C")[,colnames(modelmatrix)]
    D = coef(object,"D")[,colnames(modelmatrix)]
    
    ## fitted values
    predC = safeContrastDP(modelmatrix, C)
    predD = safeContrastDP(modelmatrix, D)
    
    ## variance of prediction at fit
    contrCovC = apply(vcov(object, "C")[colnames(modelmatrix), colnames(modelmatrix),],3,function(x) diag_safeContrastQF(modelmatrix, x))
    contrCovD = apply(vcov(object, "D")[colnames(modelmatrix), colnames(modelmatrix),],3,function(x) diag_safeContrastQF(modelmatrix, x))
    #dimnames=list(rownames(C),rownames(modelmatrix),rownames(modelmatrix)),dim=c(nrow(contrCovC),ncol(contrCovC),ncol(contrCovC)))
    #dimnames=list(rownames(D),rownames(modelmatrix),rownames(modelmatrix)),dim=c(nrow(contrCovD),ncol(contrCovD),ncol(contrCovD)))

    #predC = array(uncomplexify(predC),dimnames=list(rownames(predC),colnames(predC)),dim=c(nrow(predC),ncol(predC)))
    #predD = array(uncomplexify(predD),dimnames=list(rownames(predD),colnames(predD)),dim=c(nrow(predD),ncol(predD)))
    #predC=melt(predC)
    #contrCovC=melt(sqrt(contrCovC))
    #predD=melt(predD)
    #contrCovD=melt(sqrt(contrCovD))
    # colnames(contrCovC) = c("primerid","sample","seC")
    # colnames(contrCovD) = c("primerid","sample","seD")
    # colnames(predC) = c("primerid","sample","muC")
    # colnames(predD) = c("primerid","sample","etaD")
    m = data.table(muC = as.vector(predC), etaD = as.vector(predD), seC = as.vector(contrCovC), seD = as.vector(contrCovD), CJ(sample = rownames(modelmatrix), primerid = rownames(D)))
    m
}


if(getRversion() >= "2.15.1") globalVariables(c('muC', 'etaD', 'seC'))
#' impute missing continuous expression for plotting
#' 
#' If there are no positive observations for a contrast, it is generally not estimible.
#' However, for the purposes of testing we can replace it with the least favorable value with respect to the contrasts that are defined.
#' @param object Output of predict
#' @param groupby Variables (column names in predict) to group by for imputation (facets of the plot)
#'
#' @return data.table
#' @export
#' @examples
#' ##See stat_ell
#' example(stat_ell)
impute <- function(object,groupby){
    setDT(object)
    object[,missing:=any(is.na(muC))&!is.na(etaD),eval(groupby)]
    object[missing==TRUE,muC:=replace(muC,is.na(muC),mean(muC,na.rm=TRUE)),eval(groupby)]
    object[missing==TRUE&!is.nan(muC),seC:=replace(seC,is.na(seC),max(abs(na.omit((muC[is.na(seC)]-c(muC-seC,muC+seC)))))),eval(groupby)]
    object[,missing:=NULL]
    return(na.omit(object))
}
                                        # impute=function(object,groupby){
                                        #       setDT(object)
                                        #       object[,missing:=any(is.na(muC))&!is.na(etaD),eval(groupby)]
                                        #       object[missing==TRUE,muC:=replace(muC,is.na(muC),mean(muC,na.rm=TRUE)),eval(groupby)]
                                        #       object[missing==TRUE,seC:=replace(seC,is.na(seC),{browser();sum(seC,na.rm=TRUE)}),eval(groupby)]
                                        #       object[,missing:=NULL]
                                        #       return(na.omit(object))
                                        # }
