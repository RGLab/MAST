##' @param cl a \code{cluster} object created by \code{makeCluster}
##' @param zlmfit class \code{ZlmFit}
##' @param R number of bootstrap replicates
##' @return array of bootstrapped coefficients
##' @describeIn bootVcov1 parallel version of bootstrapping
##' @export
pbootVcov1<-function (cl,zlmfit, R = 99)
{
    sca <- zlmfit@sca
    N <- ncol(sca)
    LMlike <- zlmfit@LMlike
    if('exprs_values' %in% slotNames(zlmfit)) { 
        exprs_values = zlmfit@exprs_values
    } else{
        exprs_values = assay_idx(sca)$aidx
    }
    parallel::clusterEvalQ(cl,require(MAST))
    ## clusterEvalQ(cl,require(abind))
    parallel::clusterExport(cl,"N",envir=environment())
    parallel::clusterExport(cl,"LMlike",envir=environment())
    parallel::clusterExport(cl,"sca",envir=environment())
    parallel::clusterExport(cl,"exprs_values",envir=environment())
    manyvc <- parallel::parSapply(cl,1:R, function(i,...){
        s <- sample(N, replace = TRUE)
        newsca <- sca[, s]
        LMlike <- update(LMlike, design=colData(newsca))
        zlm(sca = newsca, LMlike = LMlike, onlyCoef = TRUE, exprs_values = exprs_values)
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
##' @param boot_index list of indices to resample.  Only one of R or boot_index can be offered.
##' @return array of bootstrapped coefficients
##' @importFrom plyr raply
##' @examples
##' data(vbetaFA)
##' zlmVbeta <- zlm(~ Stim.Condition, subset(vbetaFA, ncells==1)[1:5,])
##' #Only run 3 boot straps, which you wouldn't ever want to do in practice...
##' bootVcov1(zlmVbeta, R=3)
##' @export
bootVcov1 <- function(zlmfit, R=99, boot_index = NULL){
    sca <- zlmfit@sca
    N <- ncol(sca)
    LMlike <- zlmfit@LMlike
    if('exprs_values' %in% slotNames(zlmfit)) { 
       exprs_values = zlmfit@exprs_values
    } else{
        exprs_values = assay_idx(sca)$aidx
    }
    if(!is.null(boot_index)){
        if(!missing(R)) stop("Provide only one of `boot_index` or `R`.")
        if(!is.list(boot_index) || length(boot_index) == 0L ||
           any(lengths(boot_index) == 0L) ||
           any(vapply(boot_index, function(i) !is.numeric(i) || anyNA(i) ||
                      any(!is.finite(i)) || any(i != as.integer(i)) || any(i < 1 | i > N), logical(1)))){
            stop("boot_index must be a list of integer vectors between 1 and ", N)
        }
    } else if(is.numeric(R)){
        boot_index = lapply(1:R, function(i) sample(N, replace = TRUE))
    } else{
        stop("Provide one of `boot_index` or `R`.")
    }
    
    manyvc <- laply(boot_index, function(s){
        newsca <- sca[,s]
        LMlike <- update(LMlike, design=colData(newsca), keepDefaultCoef = TRUE)
        zlm(sca=newsca, LMlike = LMlike, onlyCoef = TRUE, exprs_values = exprs_values)
    })

    manyvc
    
}

