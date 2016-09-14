##' Bootstrap a zlmfit
##'
##' Sample cells with replacement to find bootstrapped distribution of coefficients
##' @param cl a \code{cluster} object created by \code{makeCluster}
##' @param zlmfit class \code{ZlmFit}
##' @param R number of bootstrap replicates
##' @return array of bootstrapped coefficients
##' @export
pbootVcov1<-function (cl,zlmfit, R = 99)
{
    sca <- zlmfit@sca
    N <- ncol(sca)
    LMlike <- zlmfit@LMlike
    parallel::clusterEvalQ(cl,require(MAST))
    ## clusterEvalQ(cl,require(abind))
    parallel::clusterExport(cl,"N",envir=environment())
    parallel::clusterExport(cl,"LMlike",envir=environment())
    parallel::clusterExport(cl,"sca",envir=environment())
    manyvc <- parallel::parSapply(cl,1:R, function(i,...){
        s <- sample(N, replace = TRUE)
        newsca <- sca[, s]
        LMlike <- update(LMlike, design=colData(newsca))
        zlm.SingleCellAssay(sca = newsca, LMlike = LMlike, onlyCoef=TRUE)
    })
    
    d<-dim(coef(zlmfit,"D"))
    manyvc<-aperm(array(manyvc,c(d,2,R)),c(4,1,2,3))
    dimnames(manyvc)<-c(list(NULL),dimnames(coef(zlmfit,"D")),list(c("C","D")))
    manyvc
}

##' Bootstrap a zlmfit
##'
##' Sample cells with replacement to find bootstrapped distribution of coefficients
##' @param zlmfit class \code{ZlmFit}
##' @param R number of bootstrap replicates
##' @return array of bootstrapped coefficients
##' @importFrom plyr raply
##' @export
bootVcov1 <- function(zlmfit, R=99){
    sca <- zlmfit@sca
    N <- ncol(sca)
    LMlike <- zlmfit@LMlike
    manyvc <- raply(R, {
        s <- sample(N, replace=TRUE)
        newsca <- sca[,s]
        LMlike <- update(LMlike, design=colData(newsca))
        zlm.SingleCellAssay(sca=newsca, LMlike=LMlike, onlyCoef=TRUE)
    })

    manyvc
    
}

