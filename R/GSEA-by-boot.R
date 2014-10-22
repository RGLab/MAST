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
## 3. Combine 

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
##' x0 = matrix(1:10, ncol=2)
##' y0 = colMeans(x0)
##' dim(y0) = c(2,1)
##' x1 = applyFlat(x0,y0)
##' stopifnot(colMeans(x1)==0)
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
##' dx = Drop(x, 1)
##' stopifnot(all(dim(dx)==c(2,1,2)))
Drop <- function(x, d){
    dim(x) <- dim(x)[-d]
    x
}

##' Gene set analysis for hurdle model
##'
##' Modules defined in \code{sets} are tested for average differences in expression from the "average" gene.
##' By using bootstraps, the between-gene covariance of terms in the hurdle model is found, and is used to adjust for coexpression between genes
##' 
##' \code{control} is a list with elements \code{n_randomize}, giving the number of genes to sample to approximate the non-module average expression. Default 1000.  Set to \code{Inf} to turn off the approximation.
##'
##' A 4D array is returned, with  dimensions "set" (each module), "comp" ('disc'rete or 'cont'inuous), "metric" ('stat' gives the average of the coefficient, 'var' gives the variance of that average, 'dof' gives the number of genes in the set), "group" ('test' for the genes in test-set, "null" for all genes outside the test-set).
##' @param zFit object of class ZlmFit
##' @param boots bootstraps of zFit
##' @param sets list of indices of genes
##' @param hypothesis a \code{Hypothesis} to test. Currently only one degree \code{CoefficientHypothesis} are supported
##' @param control list of control parameters.  See details.
##' @return 4D array.  See details.
##' @import abind
##' @importFrom plyr aaply
##' @export
##' @examples
##' data(vbetaFA)
##' vb1 = subset(vbetaFA, ncells==1)
##' vb1 = vb1[,freq(vb1)>.1]
##' zf = zlm.SingleCellAssay(~Stim.Condition, vb1)
##' boots = bootVcov1(zf, 10)
##' sets=list(A=1:5, B=3:10, C=15, D=1:5)
##' gsea=gseaAfterBoot(zf, boots, sets, CoefficientHypothesis('Stim.ConditionUnstim'))
##' dimnames(gsea)
##' calcZ(gsea)
##' stopifnot(all.equal(gsea['A',,,],gsea['D',,,]))
##' stopifnot(all.equal(gsea['C','cont','stat','test'], coef(zf, 'C')[15,'Stim.ConditionUnstim']))
##' stopifnot(all.equal(gsea['C','cont','stat','test'], coef(zf, 'C')[15,'Stim.ConditionUnstim']))
gseaAfterBoot <- function(zFit, boots, sets, hypothesis, control=list(n_randomize=1000)){

    ## Basic idea is to average statistics (based on coefficients defined in Zfit) and find the variance of that average using the bootstraps
    ## However, don't want to naively find the covariance across all genes then keep summing up different terms in it, because that will have quadratic complexity
    ## Instead we update the sums as necessary depending on the overlap between genes in the modules
    ## An efficiency might be to accumulate the sum of the crossproducts as we interate over genes (this corresponds to the matrix product 1*boots*t(boots)*1, where 1 is the 1-vector), rather than explicitly forming boots * t(boots), then summing
    
    n_randomize <- control$n_randomize
    stopifnot(inherits(hypothesis, 'CoefficientHypothesis'))
    hypothesis <- SingleCellAssay:::generateHypothesis(hypothesis, colnames(zFit@coefD))
    testIdx <- hypothesis@transformed
    
    ## put bootstrap replicates last
    boots <- aperm(boots, c(2,3,4,1))
    dimb <- setNames(dim(boots), c('genes', 'coef', 'comp', 'rep'))
    dnb <- setNames(dimnames(boots), names(dimb))
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

    ## this will need to be abstracted if we want to support arbitrary contrasts
    CC <- coef(zFit, 'C')
    CD <- coef(zFit, 'D')
    tstat <- cbind(C=CC[,testIdx],
                   D=CD[,testIdx])
    bootstat <- boots[,testIdx,,]
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
        vstat <- getVstat(idx, idx)
        abind(t=tstat, v=vstat, rev.along=0)
    }

    ## sum of (idx,jdx) block of covariance, and the DOF in that sum
    getVstat <- function(idx, jdx){
        ## this ends up being quite awkward because we only want to drop one dimension
        bsi <- bootstat[idx,,,drop=FALSE]
        bsj <- bootstat[jdx,,,drop=FALSE]
        vstat <- aaply(c('C', 'D'), 1, function(comp){
            sum(tcrossprod(Drop(bsi[,comp,,drop=FALSE], 2),
                           Drop(bsj[,comp,,drop=FALSE], 2)))/dimb['rep']
        })
        dofi <- colSums(!naAny[idx,,drop=FALSE])
        dofj <- colSums(!naAny[jdx,,drop=FALSE])
        matrix(c(vstat, dofi*dofj), nrow=2, dimnames=list(stat=c('stat', 'dof'), comp=c('C', 'D')), byrow=TRUE)
    }

    ## subtract off genes that were in the test set from the statistics in the non-test set
    ## and divide the statistics and variance by the DOF (number of terms int he sum)
    scalenames2 <- c('stat', 'var', 'dof')
    scalenames1 <- c('cont', 'disc')
    scalenames3 <- c('test', 'null')
    scaleStats <- function(teststat, overlapstat, nullstat){
        ## adjust the null stats and DOF
        nullstat <- nullstat - overlapstat
        nullscaled <- cbind(nullstat['stat',,]/nullstat['dof',,], dof=nullstat['dof',,'t'])
        testscaled <- cbind(teststat['stat',,]/teststat['dof',,], dof=teststat['dof',,'t'])
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
        OT <- 2*getVstat(Oidx, nullgenes)
        OO[,,'v'] <- OO[,,'v']+OT
        TT <- getStats(Tidx)
        tests[sidx,,,] <- scaleStats(TT, OO, NN)
    }
    tests
}

##' Get Z statistics and P values after running gseaAfterBoot
##'
##' @param tests output from \code{gseaAfterBoot}
##' @return 3D array with dimensinos set (modules) comp ('cont'inuous or 'disc'rete) and metric ('Z' stat or two sided 'P' value that P(z>|Z|))
##' @export
calcZ <- function(tests){
    Z <- (tests[,,'stat','test']-tests[,,'stat','null'])/sqrt(tests[,,'var','test']+tests[,,'var','null'])
    P <- pnorm(abs(Z), lower.tail=FALSE)
    ab <- abind(Z=Z, P=P, rev.along=0)
    names(dimnames(ab)) <- c('set', 'comp', 'metric')
    ab
}
