#' Return predictions from a ZlmFit object.
#'
#' @param object A \code{ZlmFit}
#' @param newdata The data to predict from. Currently ignored, will use the data in the object.
#' @param modelmatrix The model matrix specifying the linear combination of coefficients.
#' @param sepinto Separate the sample column into groups with these named columns
#' @param groupby variables to facet by for the plot, will impute missing based on facets.
#' @return Predictions and standard errors.
#' @export
#'
predict.ZlmFit <- function(object,newdata = NULL, modelmatrix=NULL){
	C = coef(object,"C")[,colnames(modelmatrix)]
	D = coef(object,"D")[,colnames(modelmatrix)]
	C = MAST:::complexifyNA(C)
	D = MAST:::complexifyNA(D)
	predC = C%*%t(modelmatrix)
	predD = D%*%t(modelmatrix)
	contrCovC = array(apply(vcov(object, "C")[colnames(modelmatrix), colnames(modelmatrix), ],3,function(x)(((modelmatrix) %*% complexifyNA(x) %*% t(modelmatrix)))),c(ncol(modelmatrix),ncol(modelmatrix),nrow(C)))
	contrCovD = array(apply(vcov(object, "D")[colnames(modelmatrix), colnames(modelmatrix), ],3,function(x)(((modelmatrix) %*% complexifyNA(x) %*% t(modelmatrix)))),c(ncol(modelmatrix),ncol(modelmatrix),nrow(C)))
	contrCovC=aperm(contrCovC,c(3,1,2))
	contrCovD=aperm(contrCovD,c(3,1,2))
	contrCovC=array(MAST:::uncomplexify(contrCovC),dimnames=list(rownames(C),rownames(modelmatrix),rownames(modelmatrix)),dim=c(nrow(contrCovC),ncol(contrCovC),ncol(contrCovC)))
	contrCovD=array(MAST:::uncomplexify(contrCovD),dimnames=list(rownames(D),rownames(modelmatrix),rownames(modelmatrix)),dim=c(nrow(contrCovD),ncol(contrCovD),ncol(contrCovD)))
	predC = array(MAST:::uncomplexify(predC),dimnames=list(rownames(predC),colnames(predC)),dim=c(nrow(predC),ncol(predC)))
	predD = array(MAST:::uncomplexify(predD),dimnames=list(rownames(predD),colnames(predD)),dim=c(nrow(predD),ncol(predD)))
	predC=melt(predC)
	contrCovC=melt(aaply(contrCovC,1,function(x)sqrt(diag(x))))
	predD=melt(predD)
	contrCovD=melt(aaply(contrCovD,1,function(x)sqrt(diag(x))))
	colnames(contrCovC) = c("primer","sample","seC")
	colnames(contrCovD) = c("primer","sample","seD")
	colnames(predC) = c("primer","sample","muC")
	colnames(predD) = c("primer","sample","muD")
   m = merge(merge(merge(predC,contrCovC,by=c("primer","sample")),predD,by=c("primer","sample")),contrCovD,by=c("primer","sample"))
   setDT(m)
   m
}


#' impute missing continuous expression for plotting
#' 
#' If there are no positive observations for a contrast, it is generally not estimible.
#' However, for the purposes of testing we can replace it with the least favorable value with respect to the contrasts that are defined.
#' @param object Output of predict
#' @param groupby Variables (column names in predict) to group by for imputation (facets of the plot)
#'
#' @return data.table with missing con
#' @export
impute <- function(object,groupby){
	setDT(object)
	object[,missing:=any(is.na(muC))&!is.na(muD),eval(groupby)]
	object[missing==TRUE,muC:=replace(muC,is.na(muC),mean(muC,na.rm=TRUE)),eval(groupby)]
	object[missing==TRUE,seC:=replace(seC,is.na(seC),sum(seC,na.rm=TRUE)),eval(groupby)]
	object[,missing:=NULL]
	return(na.omit(object))
}
