


## Approach 1: simulate
## Inputs: Pb, distribution of parameters; F, formula; C, contrasts, f(X| beta), covariate distribution
## For s in 1...Ns:
##     beta[s] ~ Pb
##     X[s] ~ f(beta[s])
##     Y := simulate(X[s], beta[s])
##     fit := zlm((X,Y), F))
##     lambda[s] := lrTest(fit, C)
##  tabulate lambda to find power as a function of (F, Pb, f(X|beta))


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
##' @param Xframe function returning a `design matrix`, (if the design is random), or a fixed design matrix.
##' @param XFormula formula used to translate the design into a model matrix. See details.
##' @param Formula formula used to model the design
##' @param contrasts argument passed to `lrTest`
##' @param zlmargs arguments passed to `zlm.SingleCellAssay`
##' @return matrix containing likelihood ratio stats, p values and dof
##' @export
##' @examples
##' Pclass <- .5
##' Nind <- 3
##' Ncell <- 96
##' Xframe <- data.frame(ind=factor(rep(seq_len(Nind), each=Ncell)), hiv=rbinom(Nind*Ncell, 1, Pclass))
##' betaD <- betaC <- c('(Intercept)'=10, ind2=-.5, ind3=.5, hiv=.5)
##' betaD['(Intercept)'] <- 0
##' theta <-  list(C=betaC, D=betaD, sigmaC=1)
##' getLambda(10, theta, Xframe, XFormula=~hiv+ind, ~hiv, 'hiv', zlmargs=list(method='bayesglm'))
getLambda <- function(nsim, thetafun, Xframe, XFormula, Formula, contrasts, zlmargs=NULL){
    lambdaDf <- matrix(NA, nsim, 3, dimnames=list(NULL, c('lambda', 'df', 'Pr(Chisq)')))
    for( s in seq_len(nsim)){
        theta <- if(is.function(thetafun)) thetafun() else thetafun
        stopifnot(c('C', 'D', 'sigmaC') %in% names(theta) )
        design <- if(is.function(Xframe)) Xframe(theta) else Xframe
        X <- model.matrix(XFormula, design)
        Y <- as.matrix(simulate1Unit(X, theta$C, theta$D, theta$sigmaC))
        sca <- suppressMessages(hushWarning(FromMatrix('SingleCellAssay', Y, design), 'no wellKey'))
        tt <- try({
        zlmfit <- suppressMessages(do.call(zlm.SingleCellAssay, c(list(formula=Formula, sca=sca), zlmargs)))
        lrt <- suppressMessages(lrTest(zlmfit, contrasts))
        lambdaDf[s,] <- lrt[,'hurdle',]
    })
        if(inherits(tt, 'try-error')) stop("Errors encountered on simulation", s)
    }
    lambdaDf
}




##' Estimate parameters for power studies from pilot data
##'
##' This function is intended to facilitate the use of pilot/related studies to fix some parameters in a power analysis.  It returns a list of single cell expression parameters from a combination of pilot data, in the form of a fitted zlm object \code{ZlmFit} and fixed parameters using a formula-like interface.
##'
##' @details
##' Each term in the \code{Formula} must be enclosed by a functional keyword. The functional keyword \code{est(term)}, requests estimation of continuous/discrete expression parameters specified in \code{term} from pilot data.  Hence \code{term} must be present in the formula used to fit \code{ZlmFit}.  If \emph{term} names a factor that was expanded in terms of multiple contrasts, then each contrast will be estimated.
##'
##'  The functional keyword \code{fixed(contrast, C=continuous_parameter, D=discrete_parameter)}, fixes the continuous/discrete parameters for the particular \emph{contrast} to specified values.  So if \code{contrast} specifies a term from a factor, it should be named according to the contrast applied.  See the example.
##'
##' Either \code{term} or \code{contrast} can be set to the special values `1`, which is the intercept term, and `error`, which is the error term for the continuous regression.
##'
##' Terms that were present in the ZlmFit formula (or newformula, if specified) but were omitted in the \code{Formula} are set to 0, except the continuous residual standard error is set to 1.
##'
##' Since `ZlmFit` may contains many genes, the parameter values are reduced by applying \code{reduceMethod} to the coefficient vector and its standard errors.  The default method returns a weighted average, inversely weighted by the standard errors.
##' @param Formula a \code{formula} evaluated according to a special syntax. Each term must be enclosed by \code{est} or \code{fixed}. See Details.
##' @param ZlmFit a fitted object that must at least contain terms that are to be `est`imated.
##' @param newdata an optional model frame used to construct the parameter vector.  If omitted then \code{ZlmFit} is used.  If specified, then so must \code{newformula}.
##' @param newformula an optional formula used to construct the parameter vector.   If omitted then \code{ZlmFit} is used.  If specified, then so must \code{newdata}.
##' @param reduceMethod a function, accepting M x P matrices of coefficients and standard errors, returning a px1 vector of parameters.  Considerable flexibility is possible by having reduceMethod return an unevaluated expression.  See \link{reduceSampleInverseSE} for an example.
##' @return a list of class ZlmTheta, containing components `C`, `D` and `sigmaC`.
##' @export
##' @examples
##' data(vbetaFA)
##' ZlmFit <- zlm.SingleCellAssay(~Population+Subject.ID+ncells, vbetaFA, method='bayesglm')
##' estimateZlmTheta(~est(1) + est(Population) + fixed(Subject.IDSub02, C=0, D=.5) + est(error), ZlmFit)
estimateZlmTheta <- function(Formula, ZlmFit, newdata, newformula, reduceMethod=reduceWeightedMean){
    tf <- terms.formula(Formula, specials=c('est', 'fixed'))
    tf.est <- attr(tf, 'specials')$est
    tf.fix <- attr(tf, 'specials')$fixed
    if(!missing(newformula)){
        pnames <- colnames(model.matrix(newformula, newdata))
    } else{
        pnames <- colnames(coef(ZlmFit, 'D'))
    }
    thetacomp <- setNames(rep(0, length(pnames)), pnames)
    theta <- list(C=thetacomp, D=thetacomp, sigmaC=1)
    estTerms <- terms(ZlmFit@LMlike@formula)
    estAssign <- attr(model.matrix(ZlmFit@LMlike), 'assign')
    
    paramlist <- list()

    for(t in c(tf.est, tf.fix)){
        thisest <- attr(tf, 'variables')[[t+1]] #index into the parse tree.  `list` is slot 1
        estimate <- t %in% tf.est
        if(estimate && !( length(thisest)==2))# && is.symbol(thisest[[2]])))
            stop("Should be simple call")
        if(!estimate) message('##validate call somehow')
        expr <- thisest[[2]]            #argument to `est` or variable name in `fixed`
        dexpr <- paste0(deparse(expr, width.cutoff=500), collapse='')
        if(dexpr == '1') dexpr <- '(Intercept)'
        if(dexpr == '(Intercept)'){
            thisestTerm <- 0
        } else if(dexpr=='error'){
            theta[['sigmaC']] <- if(estimate) reduceMethod(ZlmFit@dispersion[,'C', drop=FALSE], ZlmFit@dispersion[,'C', drop=FALSE]/sqrt(ZlmFit@df.resid[,'C'])) else eval(thisest[[3]])
            next
        }else{
            thisestTerm <- match(dexpr, labels(estTerms))
        }

        if(estimate){
            thisAssign <- which(estAssign==thisestTerm)
            paramlist[[dexpr]] <- sapply(c('C', 'D'), function(comp){
                reduceMethod(coef(ZlmFit, comp)[,thisAssign, drop=FALSE], se.coef(ZlmFit, comp)[, thisAssign, drop=FALSE])
                }, simplify=FALSE)
        } else{
            paramlist[[dexpr]] <- as.list(thisest[3:length(thisest)])
        }
    }

    ## recast parameter estimates into C/D components with correct order.
    for(p in seq_along(paramlist)){
        cn <- names(paramlist[[p]][['C']])
        if(is.null(cn)) cn <- names(paramlist)[p]
        theta[['C']][cn] <- paramlist[[p]]$C
        theta[['D']][cn] <- paramlist[[p]]$D
    }
    class(theta) <- 'ZlmTheta'
    theta
}

##' @param coef coefficient matrix
##' @param secoef standard error matrix
##' @describeIn estimateZlmTheta reduce a matrix of coefficients and standard errors to the weighted columns means
##' @export
reduceWeightedMean <- function(coef, secoef){
    index <- colnames(coef)
    if(is.null(index)) index <- seq_along(coef)
    sapply(colnames(coef),     function(s) weighted.mean(coef[,s], 1/secoef[,s]^2, na.rm=TRUE))
}

##' @describeIn estimateZlmTheta sample from a matrix of coefficients, inversely weighted by the standard errors.  Returns a list of unevaluated expressions, one for each column.
##' @export
reduceSampleInverseSE <- function(coef, secoef){
    index <- colnames(coef)
    if(is.null(index)) index <- seq_along(coef)
    sapply(index, function(s){
        nacoef <- is.na(coef[,s]) | is.na(secoef[,s])
        tcoef <- coef[!nacoef,s]
        tsecoef <- secoef[!nacoef,s]
        bquote(sample(.(tcoef), size=1, prob=1/.(tsecoef)^2))
    })
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

