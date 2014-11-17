##' Bootstrap a zlmfit
##'
##' Sample cells with replacement to find bootstrapped distribution of coefficients
##' @param cl a \code{cluster} object created by \code{makeCluster}
##' @param zlmfit class \code{ZlmFit}
##' @param R number of bootstrap replicates
##' @return array of bootstrapped coefficients
##' @importFrom parallel parSapply
##' @export
pbootVcov1<-function (cl,zlmfit, R = 999)
{
    sca <- zlmfit@sca
    N <- nrow(sca)
    LMlike <- zlmfit@LMlike
    clusterEvalQ(cl,require(SingleCellAssay))
    clusterEvalQ(cl,require(abind))
    clusterExport(cl,"N",envir=environment())
    clusterExport(cl,"LMlike",envir=environment())
    clusterExport(cl,"sca",envir=environment())
    manyvc <- parSapply(cl,1:R, function(i,...){
        s <- sample(N, replace = TRUE)
        newsca <- sca[s, ]
        LMlike <- update(LMlike, design=cData(newsca))
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
bootVcov1 <- function(zlmfit, R=999){
    sca <- zlmfit@sca
    N <- nrow(sca)
    LMlike <- zlmfit@LMlike
    manyvc <- raply(R, {
        s <- sample(N, replace=TRUE)
        newsca <- sca[s,]
        LMlike <- update(LMlike, design=cData(newsca))
        zlm.SingleCellAssay(sca=newsca, LMlike=LMlike, onlyCoef=TRUE)
    })

   manyvc
    
}

## use dfbetas? Incomplete...
bootVcov2 <- function(zlmfit, R=999){
    sca <- zlmfit@sca
    N <- nrow(sca)
    z <- zlm.SingleCellAssay(sca=newsca, LMlike=zlmfit@LMlike, hook=function(lml){
        cont <- lml@fitC
        disc <- lml@fitD
        class(disc) <-class(cont) <- 'glm'
        dfc <- lm.influence(cont)$coefficients
        dfd <- lm.influence(disc)$coefficients
        list(cont=dfc, disc=dfd)
        })
    ## do stuff to this
}


