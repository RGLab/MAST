## General procedure ##
#######################
## For set s in 1,...,Nset
## Let G(s) be genes in set s
## Get summary stats on G(s)
## q[G(s)], r[G(s)]: sum of discrete and continuous t statistics
## Q[G(s)], R[G(s)]: sum-covariance matrices
## => Modules will be compared by regulation patterns

## Alternate procedures ##
##########################
## 1. Use residual/influence function to infer covariance patterns
##    This would eliminate the need to bootstrap
## 2. Standardize bootstrap coefficients by asymptotic gene covariance
##    This might make the bootstrap more efficient 

## Validation ##
################
## Look for most parsimonous grouping?
## Look for within-group significance (groups with high % individually significant genes better)?
## Look for clustering?


##' Apply a vectorized binary operation recycling over last dimension
##'
##' When x is an array of order K, and y is an array of order K-1,
##' whose dimensions otherwise agree,
##' apply FUN by recycling y as necessary over dimension K of x.
##' @param x array, order K
##' @param y array, order K-1
##' @param FUN vectorized binary operation
##' @return array, order K equal to FUN(x,y)
##' @examples
##' ##Dumb example, could be done with scale(...,scale=FALSE)
##' x0 = matrix(1:10, nrow=2)
##' y0 = rowMeans(x0)
##' dim(y0) = c(1, 2)
##' x1 = MAST:::applyFlat(x0,y0)
##' stopifnot(rowMeans(x1)==0)
applyFlat <- function(x, y, FUN="-"){
    dx <- dim(x)
    dn <- dimnames(x)
    dy <- dim(y)
    
    newrowdim <- prod(dx[-length(dx)])
    stopifnot(newrowdim==prod(dy))
    newcoldim <- dx[length(dx)]
    dim(x) <- c(newrowdim, newcoldim)
    dim(y) <- NULL
    x <- eval(call(FUN, x, y))
    dim(x) <- dx
    dimnames(x) <- dn
    x
}

##' Drop specified dimension from an array
##'
##' Like drop(x) but only dropping specified dimensions.
##' There is no testing that the specified dimensions are actually singletons.
##' @param x array of at least d dimensions
##' @param d dimension(s) to drop
##' @return array x
##' @examples
##' x = array(1:4, dim=c(1, 2, 1, 2))
##' dx = MAST:::Drop(x, 1)
##' stopifnot(all(dim(dx)==c(2,1,2)))
##' 
Drop <- function(x, d){
    dim(x) <- dim(x)[-d]
    x
}

##' Gene set analysis for hurdle model
##'
##' Modules defined in \code{sets} are tested for average differences in expression from the "average" gene.
##' By using bootstraps, the between-gene covariance of terms in the hurdle model is found, and is used to adjust for coexpression between genes.
##' We drop genes if the coefficient we are testing was not estimible in original model fit in \code{zFit} or in any of the bootstrap replicates (evidenced an \code{NA} in the bootstrap array).  This might yield overly conservative inference.
##' Since bootstrapping is a randomized procedure, the degrees of freedom of a module (and its variance parameters) might differ from run-to-run.
##' You might try setting \code{var_estimate='modelbased'} to relax this requirement by assuming independence between genes and then using the asymptotic covariance estimates, which are deterministic, but may result in overly-generous inference.
##'
##' @section \code{control}:
##' \code{control} is a list with elements:
##' \itemize{
##' \item \code{n_randomize}, giving the number of genes to sample to approximate the non-module average expression. Set to \code{Inf} to turn off the approximation (the default).
##' \item \code{var_estimate}, giving the method used to estimate the variance of the modules.  \code{bootall} uses the bootstrapped covariance matrices.  \code{bootdiag} uses only the diagonal of the bootstrapped covariance matrix (so assuming independence across genes). \code{modelbased} assumes independence across genes and uses the variance estimated from the model.}
##' @section Return Value:
##' A 4D array is returned, with  dimensions "set" (each module), "comp" ('disc'rete or 'cont'inuous), "metric" ('stat' gives the average of the coefficient, 'var' gives the variance of that average, 'dof' gives the number of genes that were actually tested in the set), "group" ('test' for the genes in test-set, "null" for all genes outside the test-set).
##' @param zFit object of class ZlmFit
##' @param boots bootstraps of zFit
##' @param sets list of indices of genes
##' @param hypothesis a \code{Hypothesis} to test. Currently only one degree \code{CoefficientHypothesis} are supported.
##' @param control parameters as provided by \code{gsea_control}.  See details.
##' @return Object of class \code{GSEATests}, containing slots \code{tests},  4D array and \code{bootR}, the number of boostrap replicates.
##' @import abind
##' @export
##' @seealso \link{calcZ}
##' @seealso summary,GSEATests-method
##' @examples
##' data(vbetaFA)
##' vb1 = subset(vbetaFA, select = ncells==1)
##' vb1 = vb1[,freq(vb1)>.1][1:15,]
##' zf = zlm(~Stim.Condition, vb1)
##' boots = bootVcov1(zf, 5)
##' sets = list(A=1:5, B=3:10, C=15, D=1:5)
##' gsea = gseaAfterBoot(zf, boots, sets, CoefficientHypothesis('Stim.ConditionUnstim'))
##' ## Use a model-based estimate of the variance/covariance.
##' gsea_mb = gseaAfterBoot(zf, boots, sets, CoefficientHypothesis('Stim.ConditionUnstim'),
##' control = gsea_control(var_estimate = 'modelbased'))
##' calcZ(gsea)
##' summary(gsea)
##' \dontshow{
##' stopifnot(all.equal(gsea@tests['A',,,],gsea@tests['D',,,]))
##' stopifnot(all.equal(gsea@tests['C','cont','stat','test'], coef(zf, 'C')[15,'Stim.ConditionUnstim']))
##' }
gseaAfterBoot <- function(zFit, boots, sets, hypothesis, control=gsea_control(n_randomize=Inf, var_estimate='bootall')){

    ## Basic idea is to average statistics (based on coefficients defined in Zfit) and find the variance of that average using the bootstraps
    ## However, don't want to naively find the covariance across all genes then keep summing up different terms in it, because that will have quadratic complexity
    ## Instead we update the sums as necessary depending on the overlap between genes in the modules
    ## An efficiency might be to accumulate the sum of the crossproducts as we interate over genes (this corresponds to the matrix product 1*boots*t(boots)*1, where 1 is the 1-vector), rather than explicitly forming boots * t(boots), then summing
    
    n_randomize <- control$n_randomize
    var_est <- control$var_estimate
    var_est <- match.arg(var_est, c('bootall', 'bootdiag', 'modelbased'))
    stopifnot(inherits(hypothesis, 'CoefficientHypothesis'))
    hypothesis <- generateHypothesis(hypothesis, colnames(zFit@coefD))
    testIdx <- hypothesis@index

    ## this will need to be abstracted if we want to support arbitrary contrasts
    CC <- coef(zFit, 'C')
    CD <- coef(zFit, 'D')
    tstat <- cbind(C=CC[,testIdx],
                   D=CD[,testIdx])
    
    if(var_est == 'modelbased'){
        ## Make a pretend bootstrap with a single replicate, for the sake of reusing the logic below
        boots = array(0, dim = c(1, # one replicate
                           nrow(CD), # genes
                           ncol(CD), #coefficients
                           2), dimnames = c(list(1), dimnames(CD), list(c('C', 'D')))) # components
    }

    ## put bootstrap replicates last
    boots <- aperm(boots, c(2,3,4,1))
    dimb <- setNames(dim(boots), c('genes', 'coef', 'comp', 'rep'))
    dnb <- setNames(dimnames(boots), names(dimb))
    if(dimb['genes']!=nrow(zFit@coefD)) stop('Bootstraps must be run on same set of genes as `zFit`')
    if(!all(is.numeric(unlist(sets)))) stop('`sets` should be indices of genes in `zFit`')
    
    ## average coefficient over replicates
    bmean <- rowMeans(boots, dims=3)
    boots <- applyFlat(boots, bmean, FUN='-')
    ## boots is now an array of centered coefficients, so we can calculate the second moment just by averaging crossproducts
    
    ## we'll take crossproducts in this set, and then update as necessary to get the mean and variance of the null set
    ## sub sampling: subsample will converge to the true sample
    ## and this will save us having to update the null set so often
    ## (will be more disjoint from test set)
    nullgenes <- if(dimb['genes']>n_randomize) sample(dimb['genes'], n_randomize, replace=FALSE) else seq_len(dimb['genes'])
    nullgenes <- sort(nullgenes)

   
    if(var_est!='modelbased'){
        ## put bootstraps of coefficient of interest into array
        bootstat <- boots[,testIdx,,]
    } else{
        ## put model-based covariances into array (with trailing dimension=1)
        dimb['rep'] <- Inf
        bootstat <- array(c(vcov(zFit, 'C')[testIdx,testIdx,],
                            vcov(zFit,'D')[testIdx,testIdx,]),
                          dim=c(dimb['genes'], 2, 1),
                          dimnames=list(genes=dnb$genes, comp=c('C', 'D'), rep='1'))
    }
    
    natstat <- is.na(tstat)
    naboots <- rowSums(is.na(bootstat), dims=2)
    ## set all bootstrap and t statistics to zero if either had an NA for a gene
    ## but keep track of where we did so and adjust the DOF accordingly
    naAny <- natstat | naboots
    tstat[naAny] <- 0
    bootstat[is.na(bootstat)] <- 0
    ## zero out NAs
    bootstat <- applyFlat(bootstat, !naAny*1, FUN='*')
    
    ## this returns the sum over idx (index set) of the coefficients, and the variance of this sum
    ## we'll average those down the line, but because we want to do online updates, it's best to maintain the sum rather than the average
    getStats <- function(idx){
        thisstat <- tstat[idx,,drop=FALSE]
        tstat <- colSums(thisstat)
        dof <- colSums(!naAny[idx,,drop=FALSE])
        ## matrix of discrete and continuous and the DOF (number genes)
        tstat <- matrix(c(tstat, dof), nrow=2,
                        dimnames=list(stat=c('stat', 'dof'), comp=c('C', 'D')), byrow=TRUE)
        ## variance of sum over idx
        vstat <- getVstat(idx, idx, returnCor=TRUE)
        vstatCor <- vstat[2,]
        vstat <- vstat[-2,]
        list(stat=abind(t=tstat, v=vstat, rev.along=0),
             cor=vstatCor)
    }

    ## methods to calculate the variance of a block of genes
    ## bootstrap based
    vstatboot <- function(idx, jdx, onlyDiag=FALSE, returnCor=FALSE){
        ## sum of (idx,jdx) block of covariance, and the DOF in that sum
        ## returnCor should always be FALSE when idx!=jdx, but we won't check that

        if(onlyDiag){
            if(!isTRUE(all.equal(idx,jdx))) idx <- jdx <- integer(0)
            dof <- colSums(!naAny[idx,,drop=FALSE])^2
        } else{
            dofi <- colSums(!naAny[idx,,drop=FALSE])
            dofj <- colSums(!naAny[jdx,,drop=FALSE])
            dof <- dofi*dofj
        }
        bsi <- bootstat[idx,,,drop=FALSE]
        bsj <- bootstat[jdx,,,drop=FALSE]
        vstat <- matrix(NA, nrow=2, ncol=2, dimnames=list(c('stat', 'avgCor'), c('C', 'D')))
        for(comp in c('C', 'D')){
            tcp <- tcrossprod(Drop(bsi[,comp,,drop=FALSE], 2),
                              Drop(bsj[,comp,,drop=FALSE], 2))
            if(onlyDiag){
                vstat['stat', comp] <- sum(diag(tcp))/dimb['rep']
            } else{
                vstat['stat', comp] <- sum(tcp)/dimb['rep']
            }

            if(returnCor && length(idx)>1){ # so we don't emit warnings or die on empty or singleton idx
                ccp <- hushWarning(cov2cor(tcp), fixed("diag(.) had 0 or NA entries"))
                vstat['avgCor', comp] <- sum(ccp[upper.tri(ccp)], na.rm=TRUE)
            }
        }   
        if(returnCor){
            return(rbind(vstat, dof=dof))
        } else{
            return(rbind(vstat[-2,,drop=FALSE], dof=dof))
        }
    }
    
    ## modelbased
    vstatmodelbased <- function(idx, jdx, returnCor=FALSE){
        if(!isTRUE(all.equal(idx,jdx))) idx <- jdx <- integer(0)
        dof <- colSums(!naAny[idx,,drop=FALSE])^2
        
        vstat <- matrix(NA, nrow=2, ncol=2, dimnames=list(c('stat', 'avgCor'), c('C', 'D')))
        for(comp in c('C', 'D')){
            vstat['stat', comp] <- sum(bootstat[idx,comp,])
        }

        if(returnCor){
            return(rbind(vstat, dof=dof))
        } else{
            return(rbind(vstat[-2,,drop=FALSE], dof=dof))
        }
    }
    

    ## define the method to calculate the variance
    if(var_est=='modelbased'){
        getVstat <- vstatmodelbased
    } else if(var_est=='bootall'){
        getVstat <- vstatboot
    } else if(var_est=='bootdiag'){
        getVstat <- function(idx, jdx, returnCor=FALSE) vstatboot(idx, jdx, onlyDiag=TRUE, returnCor)
    } else{
        stop("unreachable code")
    }

    ## subtract off genes that were in the test set from the statistics in the non-test set
    ## and divide the statistics and variance by the DOF (number of terms int he sum)
    scalenames2 <- c('stat', 'var', 'dof', 'avgCor')
    scalenames1 <- c('cont', 'disc')
    scalenames3 <- c('test', 'null')
    scaleStats <- function(test, overlap, null){
        ## adjust the null stats and DOF
        null$stat <- null$stat - overlap$stat
        nullscaled <- cbind(null$stat['stat',,]/null$stat['dof',,], dof=null$stat['dof',,'t'], cor=null$cor)
        testscaled <- cbind(test$stat['stat',,]/test$stat['dof',,], dof=test$stat['dof',,'t'], cor=test$cor)
        abind(test=testscaled, null=nullscaled, rev.along=0)
    }

    NN <- getStats(nullgenes)
    tests <- array(NA, dim=c(length(sets), length(scalenames1), length(scalenames2), length(scalenames3)), dimnames=list(set=names(sets), comp=scalenames1, metric=scalenames2, group=scalenames3))
    for(sidx in seq_along(sets)){
        Tidx <- sets[[sidx]]
        Oidx <- intersect(Tidx, nullgenes)
        OO <- getStats(Oidx)
        ## off-diagonal block for T (covariance between T and O)
        ## multiply be 2 because we want to kill TO as well as OT blocks
        OT <- 2*getVstat(Oidx, setdiff(nullgenes, Oidx))
        OO$stat[,,'v'] <- OO$stat[,,'v']+OT
        TT <- getStats(Tidx)
        tests[sidx,,,] <- scaleStats(TT, OO, NN)
    }
    new('GSEATests', tests=tests, bootR=dimb['rep'])
}

##' @export
##' @param n_randomize the number of genes to sample to approximate the non-module average expression. Set to \code{Inf} to turn off the approximation (the default).
##' @param var_estimate the method used to estimate the variance of the modules, one of \code{bootall}, \code{bootdiag}, or \code{modelbased}.
##' @describeIn gseaAfterBoot set control parameters.  See Details.
gsea_control = function(n_randomize = Inf, var_estimate = 'bootall'){
    list(n_randomize = n_randomize, var_estimate = var_estimate)
}


##match moments to get approximation of t statistic
.approxt <- function(dof){
    s <- ifelse(is.finite(dof), sqrt(dof/(dof-2)), 1) #scale of each t
    kur <- 6/(dof-4) #kurtosis of each t
    ## assume we weight the sum by inverse scale (so that we have unit variance in each var)
    skur <- sum(kur)/length(dof)^2
    ## match kurtosis of sum to dof of t approximiation
    nu <- 6/skur+4
    ## scale of sum
    ss <- sum(s^2)
    ## scale of approximation
    scale <- sqrt( (nu-2)/(nu) / #expected variance of t approx
                   ss) #variance of sum
    list(W=1/s, nu=nu, scale=scale)
}

##' Get Z or T statistics and P values after running gseaAfterBoot
##' 
##' The Z or T statistics may be reported by component (discrete/continuous) when \code{combined='no'} or combined by Fisher's or Stouffer's method (\code{combined='fisher'} or \code{combined='stouffer'}.
##' Fisher's method uses the product of the p-values, while Stouffer's method uses the sum of the Z/T scores.
##' The "Z" score returned by Fisher is the normal quantile that would yield the observed Fisher P-value, whose sign is derived from the sign of the maximum component Z score.
##' The "Z" score returned by Stouffer when \code{testType='normal'} is the sum of the Z scores, over sqrt(2).
##' When \code{testType='t'} it is a weighted combination of the Z scores, with weights correponding to the degrees of freedom in each of the t statistics.
##' A t-approximation to this sum of t-variables is derived by matching moments.  It seems to be fairly accurate in practice.
##' @param gseaObj output from \code{gseaAfterBoot}
##' @param testType either 'normal' or 't'.  The 't' test adjusts for excess kurtosis due to the finite number of bootstrap replicates used to estimate the variance of the statistics.  This will result in more conservative inference.
##' @param combined \code{character} one of 'none', 'fisher' or 'stouffer'
##' @return 3D array with dimensions set (modules) comp ('cont'inuous or 'disc'rete) and metric ('Z' stat and two sided 'P' value that P(z>|Z|)) if \code{combined='no'}, otherwise just a matrix.
##' @seealso gseaAfterBoot
##' @examples
##' ## See the examples in gseaAfterBoot
##' example(gseaAfterBoot)
##' @export
calcZ <- function(gseaObj, testType='t', combined='none'){
    if(!inherits(gseaObj, 'GSEATests')) stop('`gseaObj` must inherit from `GSEAtests`')
    tests <- gseaObj@tests
    bootR <- gseaObj@bootR
    testType <- match.arg(testType, c('t', 'normal'))
    combined <- match.arg(combined, c('none', 'fisher', 'stouffer'))
    Z <- (tests[,,'stat','test']-tests[,,'stat','null'])/sqrt(tests[,,'var','test']+tests[,,'var','null'])
    
    if(testType=='t'){
        ## satterthwaite approximation to degrees of freedom
        dof <- (bootR-1)*(tests[,,'var','test']+tests[,,'var','null'])^2/(tests[,,'var','test']^2+tests[,,'var','null']^2)
        dof[is.na(dof)] <- 0 #component we couldn't test
    } else if(testType=='normal'){
        dof <- Inf
    }
    P <- pt(abs(Z), df=dof, lower.tail=FALSE)*2
    out3d <- abind(Z=Z, P=P, rev.along=0)
    names(dimnames(out3d)) <- c('set', 'comp', 'metric')
    ntest <- rowSums(!is.na(Z))
    if(combined=='none'){
        return(out3d)
    }else if(combined =='fisher'){
        maxsign <- sign(rowSums(Z, na.rm=TRUE))
        chival <- -2*(rowSums(log(P), na.rm=TRUE))
        Pval <- pchisq(chival, df=2*ntest, lower.tail=FALSE)
        out2d <- cbind(Z=-maxsign*qnorm(Pval/2), P=Pval)
    } else if(combined =='stouffer'){
        if(testType=='normal'){
            Wmat <- matrix(1, nrow=nrow(Z), ncol=ncol(Z))
            scale <- 1/sqrt(ntest)
            dofComb <- Inf
        } else{
            tapprox <- apply(dof, 1, .approxt) #find a t approximation for sum of t stats
            Wmat <- t(sapply(tapprox, '[[', 'W'))
            scale <- sapply(tapprox, '[[', 'scale')
            dofComb <- sapply(tapprox, '[[', 'nu')
            dofComb[ntest==1] <- apply(dof[ntest==1,,drop=FALSE], 1, max)
        }
        Zcomb <- rowSums(Wmat*Z, na.rm=TRUE)*scale
        out2d <- cbind(Z=Zcomb, P=pt(abs(Zcomb), dofComb, lower.tail=FALSE)*2)
    }
    names(dimnames(out2d)) <- c('set', 'metric')
    return(out2d)
}


##' Summarize gene set enrichment tests
##'
##' Returns a \code{data.table} with one row per gene set.
##' This \code{data.table} contains columns: 
##' \describe{
#'   \item{set}{name of gene set}
#'   \item{cond_Z}{Z statistic for continuous component}
#' \item{cont_P}{wald P value}
#' \item{cont_effect}{difference in continuous regression coefficients between null and test sets (ie, the numerator of the Z-statistic.)}
#' \item{disc_Z}{Z statistic for discrete}
#' \item{disc_P}{wald P value}
#' \item{disc_effect}{difference in discrete regression coefficients between null and test sets.}
#' \item{combined_Z}{combined discrete and continuous Z statistic using Stouffer's method}
#' \item{combined_P}{combined P value}
#' \item{combined_adj}{FDR adjusted combined P value}
#' }
##' @param object A \code{GSEATests} object
##' @param ... passed to \code{calcZ}
##' @return \code{data.table}
##' @seealso gseaAfterBoot
##' @examples
##' ## See the examples in gseaAfterBoot
##' example(gseaAfterBoot)
##' @export
## The following was added only to get the package to pass BiocCheck and can be deleted hopefully in the next version
##' @importFrom reshape2 colsplit
setMethod('summary', signature=c(object='GSEATests'), function(object, ...){
    t_stat <- as.data.table(reshape2::melt(calcZ(object, combined='none', ...)))
    effect_size <- as.data.table(reshape2::melt(object@tests[,,'stat','test']-object@tests[,,'stat','null']))
    effect_size_wide <- dcast(effect_size, set ~ comp)
    setnames(effect_size_wide, c('disc', 'cont'), c('disc_effect', 'cont_effect'))
    
    t_stat_wide <- dcast(t_stat, set ~ comp + metric)
    pvalArr <- calcZ(object, combined = "stouffer", ...)
    pvals <- data.table(pvalArr)
    setnames(pvals, c('P', 'Z'), c('combined_P', 'combined_Z'))
    pvals$set <- rownames(pvalArr)
    t_stat_comb <- merge(t_stat_wide, pvals, by = "set")
    t_stat_comb <- merge(t_stat_comb, effect_size_wide, by='set')
    t_stat_comb[,combined_adj:=p.adjust(combined_P,"fdr")]
    setorder(t_stat_comb,combined_adj)
    t_stat_comb
})

if(getRversion() >= "2.15.1") globalVariables(c('combined_P', 'combined_adj', 'component'))
