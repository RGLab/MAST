#' Return predictions from a ZlmFit object.
#'
#' @param object A \code{ZlmFit}
#' @param newdata The data to predict from. Currently ignored, will use the data in the object.
#' @param modelmatrix The model matrix specifying the linear combination of coefficients.
#' @param ... ignored
#' @return Predictions and standard errors.
#' @export
#' @examples
#' ##See stat_ell
#' example(stat_ell)
predict.ZlmFit <- function(object,newdata = NULL, modelmatrix=NULL, ...){
    C = coef(object,"C")[,colnames(modelmatrix)]
    D = coef(object,"D")[,colnames(modelmatrix)]
    C = complexifyNA(C)
    D = complexifyNA(D)
    predC = C%*%t(modelmatrix)
    predD = D%*%t(modelmatrix)
    contrCovC = array(apply(vcov(object, "C")[colnames(modelmatrix), colnames(modelmatrix), ],3,function(x)(((modelmatrix) %*% complexifyNA(x) %*% t(modelmatrix)))),c(ncol(modelmatrix),ncol(modelmatrix),nrow(C)))
    contrCovD = array(apply(vcov(object, "D")[colnames(modelmatrix), colnames(modelmatrix), ],3,function(x)(((modelmatrix) %*% complexifyNA(x) %*% t(modelmatrix)))),c(ncol(modelmatrix),ncol(modelmatrix),nrow(C)))
    contrCovC=aperm(contrCovC,c(3,1,2))
    contrCovD=aperm(contrCovD,c(3,1,2))
    contrCovC=array(uncomplexify(contrCovC),dimnames=list(rownames(C),rownames(modelmatrix),rownames(modelmatrix)),dim=c(nrow(contrCovC),ncol(contrCovC),ncol(contrCovC)))
    contrCovD=array(uncomplexify(contrCovD),dimnames=list(rownames(D),rownames(modelmatrix),rownames(modelmatrix)),dim=c(nrow(contrCovD),ncol(contrCovD),ncol(contrCovD)))
    predC = array(uncomplexify(predC),dimnames=list(rownames(predC),colnames(predC)),dim=c(nrow(predC),ncol(predC)))
    predD = array(uncomplexify(predD),dimnames=list(rownames(predD),colnames(predD)),dim=c(nrow(predD),ncol(predD)))
    predC=melt(predC)
    contrCovC=melt(aaply(contrCovC,1,function(x)sqrt(diag(x))))
    predD=melt(predD)
    contrCovD=melt(aaply(contrCovD,1,function(x)sqrt(diag(x))))
    colnames(contrCovC) = c("primerid","sample","seC")
    colnames(contrCovD) = c("primerid","sample","seD")
    colnames(predC) = c("primerid","sample","muC")
    colnames(predD) = c("primerid","sample","muD")
    m = merge(merge(merge(predC,contrCovC,by=c("primerid","sample")),predD,by=c("primerid","sample")),contrCovD,by=c("primerid","sample"))
    setDT(m)
    m
}


if(getRversion() >= "2.15.1") globalVariables(c('muC', 'muD', 'seC'))
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
    object[,missing:=any(is.na(muC))&!is.na(muD),eval(groupby)]
    object[missing==TRUE,muC:=replace(muC,is.na(muC),mean(muC,na.rm=TRUE)),eval(groupby)]
    object[missing==TRUE&!is.nan(muC),seC:=replace(seC,is.na(seC),max(abs(na.omit((muC[is.na(seC)]-c(muC-seC,muC+seC)))))),eval(groupby)]
    object[,missing:=NULL]
    return(na.omit(object))
}
                                        # impute=function(object,groupby){
                                        #       setDT(object)
                                        #       object[,missing:=any(is.na(muC))&!is.na(muD),eval(groupby)]
                                        #       object[missing==TRUE,muC:=replace(muC,is.na(muC),mean(muC,na.rm=TRUE)),eval(groupby)]
                                        #       object[missing==TRUE,seC:=replace(seC,is.na(seC),{browser();sum(seC,na.rm=TRUE)}),eval(groupby)]
                                        #       object[,missing:=NULL]
                                        #       return(na.omit(object))
                                        # }
