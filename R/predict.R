#' Return predictions from a ZlmFit object.
#'
#' @param object A \code{ZlmFit}
#' @param newdata The data to predict from. Currently ignored, will use the data in the object.
#' @param modelmatrix The model matrix specifying the linear combination of coefficients.
#'
#' @return Predictions and standard errors.
#' @export
#'
#' @examples
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
  merge(merge(merge(predC,contrCovC,by=c("primer","sample")),predD,by=c("primer","sample")),contrCovD,by=c("primer","sample"))
}

