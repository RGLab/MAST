## Approach 1: simulate
## Inputs: Pb, distribution of parameters; F, formula; C, contrasts, f(X| beta), covariate distribution
## For s in 1...Ns:
##     beta[s] ~ Pb
##     X[s] ~ f(beta[s])
##     Y := simulate(X[s], beta[s])
##     fit := zlm((X,Y), F))
##     lambda[s] := lrTest(fit, C)
##  tabulate lambda to find power as a function of (F, Pb, f(X|beta))


## These are more-or-less bespoke crafted for
## Chetan's scenario
simulateX <- function(Pclass, Nind, Ncell){
    fun <- function(betaX){
        as.matrix(cbind('(Intercept)'=1, ind=rep(seq_len(Nind), each=Ncell), hiv=rbinom(Nind*Ncell, 1, Pclass)))
    }
    fun
}

thetafun <- function(){
    betaD <- betaC <- c('(Intercept)'=10, hiv=.5, ind=.5)
    betaD['(Intercept)'] <- 0           #50% expression
    list(C=betaC, D=betaD, sigmaC=1)
}

##' Simulate zero-inflated single cell data
##'
##' @param x covariate matrix.  number of rows determines number of observations.
##' @param betaC continuous coefficient vector
##' @param betaD discrete coefficient vector
##' @param sigmaC sd of continuous component
##' @return vector of observations y
simulate1Unit <- function(x, betaC, betaD, sigmaC){
    etaC <- x %*% betaC    
    etaD <- x %*% betaD
    yC <- rnorm(length(etaC), sd=sigmaC)+etaC
    yD <- runif(length(etaD))<exp(etaD)/(1+exp(etaD))
    yC*yD
}

##' Simulate tests on a single cell gene expression data set using hurdle model
##'
##' @param nsim number of simulations to run
##' @param thetafun function returning a realization of the parameters (return a constant if parameters are fixed)
##' @param Xfun function returning a realization of the design matrix.  Must explicitly include an intercept.
##' @param Formula formula used to model the design
##' @param contrasts argument passed to `lrTest`
##' @param zlmargs arguments passed to `zlm.SingleCellAssay`
##' @return matrix containing likelihood ratio stats, p values and dof
##' @examples
##' getLambda(100, thetafun, Xfun=simulateX(.7, 3, 96), ~hiv+(1|ind), 'hiv', zlmargs=list(method='glmer'))
getLambda <- function(nsim, thetafun, Xfun, Formula, contrasts, zlmargs=NULL){
    lambdaDf <- matrix(NA, nsim, 2, dimnames=list(NULL, c('lambda', 'df', 'Pr(Chisq)')))
    for( s in seq_len(nsim)){
        theta <- thetafun()
        stopifnot(c('C', 'D', 'sigmaC') %in% names(theta) )
        X <- Xfun(theta$X)
        Y <- as.matrix(simulate1Unit(X, theta$C, theta$D, theta$sigmaC))
        Xnointercept <- as.data.frame(X[,setdiff(colnames(X), '(Intercept)')])
        sca <- suppressMessages(hushWarning(FromMatrix('SingleCellAssay', Y, Xnointercept), 'no wellKey'))
        tt <- try({
        zlmfit <- suppressMessages(do.call(zlm.SingleCellAssay, c(list(formula=Formula, data=sca), zlmargs)))
        lrt <- suppressMessages(lrTest(zlmfit, contrasts))
        lambdaDf[s,] <- lrt[,'hurdle',]
    })
        if(inherits(tt, 'try-error')) stop("Errors encountered on simulation", s)
    }
    lambdaDf
}

## Approach 2: assume normality (can solve for effect size, sample size, power, level)
##' Analytic/Asymptotic power calculations for repeated measures on single cells
##'
##' If groups A and B recieve treatments A and B, and there are `nBetween` subjects in each group, and each subject has `nWithin` cells sampled,
##' then the design is "cross-sectional"
##' If each individual receives treatments A and B, and there are `nBetween` subjects and `nWithin cells sampled *in each treatment, per individual, then
##' the design is "longitudinal"
##' In either case, the power of detecting a difference of size `logfc` at significance `sig.level` is provided
##' given within-subject standard deviation `sigWithin` (this is the average s.d. of cellular expression within in subject) and between-subject standard deviation `sigBetween` (this is the s.d. of the \emph{mean} level of expression for each subject.
##'
##' At most one of the power parameters can omitted, by setting it to \code{NULL}--
##' then it will be solved in terms of the other provided parameters.
##'
##' This analysis does not account for zero-inflation/bi-modality but could be useful as an lower bound on power when  \code{sigWithin} and \code{sigBetween} are set large enough (ie, set at zero-inflated levels).
##' @param design one of "crosssectional" or "longitudinal"
##' @param sig.level level of test
##' @param power power of test
##' @param nWithin number of cells sampled per subject (per subject per treatment for longitudinal designs)
##' @param sigWithin average s.d. of expression level within subject (within treatment, for longitudinal designs)
##' @param nBetween number of subjects (per treatment, for cross-sectional designs)
##' @param sigBetween s.d. of treatment effect, across subjects
##' @param logfc overall average log-fold change between treatment A and B
##' @param ... arguments passed to uniroot
##' @return vector of power parameters
##' @references Diggle et al 2002, chapter 2
##' @export
repeatedMeasuresAsymptoticPower <- function(design='crosssectional', sig.level=.05, power=.8, nWithin=50, sigWithin=sqrt(1.4), nBetween=16, sigBetween=sqrt(.1), logfc=NULL, ...){

    design <- match.arg(design, c('crosssectional', 'longitudinal'))

    if (sum(sapply(list(sig.level, power, nWithin, sigBetween, sigWithin, nBetween, logfc), is.null)) != 1) 
        stop("Exactly one of 'sig.level', 'power', 'nWithin', 'sigBetween', 'nBetween' and 'logfc' must be NULL.")
     
    ## Equal to zero when power is achieved
    ## Greater than zero when more sample is needed to achieve power
    ## Less than zero when less sample is needed
    cross.body <- quote({
        Za <- qnorm(1-sig.level/2)
        Zb <- qnorm(power)
        rho <- sigBetween/sqrt(sigBetween^2 + sigWithin^2)
        nBetween - 2*(Za + Zb)^2 * (1 + (nWithin-1)*rho)/ ( nWithin *(logfc/sigWithin)^2) #diggle equation 2.4.2
})
    long.body <- quote({
        Za <- qnorm(1-sig.level/2)
        Zb <- qnorm(power)
        rho <- sigBetween^2/sqrt(sigBetween^2 + sigWithin^2)
        nBetween - 2*(Za+Zb)^2*sigWithin^2*(1-rho)/(nWithin *logfc^2) / #diggle equation 2.4.1
            .25 /                           #sd of a binary variable
                2                           #because it's actually just a one-group comparison
        
    })
    
    p.body <- if(design=='crosssectional') cross.body else long.body

    ## which parameter are we solving for
    if(is.null(logfc)){ 
        logfc <- uniroot(function(logfc) eval(p.body), interval=c(1e-3, 1e3), ...)$root
    } else if(is.null(sig.level)){
        sig.level <- uniroot(function(sig.level) eval(p.body), interval=c(1e-10, 1-1e-10), ...)$root
    } else if(is.null(power)){
        power <- uniroot(function(power) eval(p.body), interval=c(1e-10, 1-1e-10), ...)$root
    } else if(is.null(nWithin)){
        nWithin <- uniroot(function(nWithin) eval(p.body), interval=c(1, 1e3), ...)$root
    }else if(is.null(sigBetween)){
        sigBetween <- uniroot(function(sigBetween) eval(p.body), interval=c(1e-3, 1e3), ...)$root
    }else if(is.null(sigWithin)){
        sigWithin <- uniroot(function(sigWithin) eval(p.body), interval=c(1e-3, 1e3), ...)$root
    }else if(is.null(nBetween)){
        nBetween <- uniroot(function(nBetween) eval(p.body), interval=c(1, 1e3), ...)$root
    } else{
        stop('Unreachable code')
    }
    
    c(sig.level=sig.level, power=power, nWithin=nWithin, sigWithin=sigWithin, nBetween=nBetween, sigBetween=sigBetween, logfc=logfc)
    }
    


## Todo:
## Approach 3: empirical asymptotics
## Inputs same as approach 1.  beta fixed.
## For n in [N]:
##     X ~ f(beta[s])
##     Y := simulate(X, beta)
##     fit := zlm((X, Y), F)
##     lambda := lrTest(fit, C)
##  tabulate lambda as a function N, suppose lambda ~ non-central-chisquare(lambda[n]) to get power
## Todo: estimate covariate and nuisance parameters from pilot data (need interface to define these).

