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
	C = coef(fit,"C")[,colnames(modelmatrix)]
	D = coef(fit,"D")[,colnames(modelmatrix)]
	C = complexifyNA(C)
	D = complexifyNA(D)
	predC = C%*%t(modelmatrix)
	predD = D%*%t(modelmatrix)
	contrCovC=aaply(.data = vcov(object,"C")[colnames(modelmatrix),colnames(modelmatrix),],.margins = 3,.fun = function(x) {((modelmatrix)%*%complexifyNA(x)%*%t(modelmatrix))})
	contrCovD=aaply(.data = vcov(object,"D")[colnames(modelmatrix),colnames(modelmatrix),],.margins = 3,.fun = function(x) {((modelmatrix)%*%complexifyNA(x)%*%t(modelmatrix))})
	contrCovC=array(uncomplexify(contrCovC),dimnames=list(rownames(contrCovC),colnames(contrCovC),colnames(contrCovC)),dim=c(nrow(contrCovC),ncol(contrCovC),ncol(contrCovC)))
	contrCovD=array(uncomplexify(contrCovD),dimnames=list(rownames(contrCovD),colnames(contrCovD),colnames(contrCovD)),dim=c(nrow(contrCovD),ncol(contrCovD),ncol(contrCovD)))
	predC = array(uncomplexify(predC),dimnames=list(rownames(predC),colnames(predC)),dim=c(nrow(predC),ncol(predC)))
	predD = array(uncomplexify(predD),dimnames=list(rownames(predD),colnames(predD)),dim=c(nrow(predD),ncol(predD)))
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

