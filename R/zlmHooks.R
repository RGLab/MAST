
residualsHook <- function(fit){
    residuals(fit, which='Marginal')
}

revealHook <- function(zlm){
    return(attr(zlm, 'hookOut'))
}

##' Residual hooks and collection methods
##'
##' After each gene is fit, a hook function can optionally be run and the output saved.
##' This allows extended computations to be done using the fitted model, without keeping it in memory.
##' Here this is used to calculate various residuals, though in some cases they can be done using only the information contained in the \code{ZlmFit}-class.
##' @importFrom plyr laply
##' @export
##' @param x \code{ZlmFit}-class
##' @param sca \code{SingleCellAssay} object to which the residuals should be added
##' @param newLayerName \code{character} name of the assay layer
##' @seealso zlm
##' @section Total residual types:
##' Each component of the model contributes several flavors of residual, which can be combined in various fashions.
##' The discrete residual can be on the response scale (thus subtracting the predicted probability of expression from the 0/1 expression value).
##' Or it can be a deviance residual, revealing something about the log-likelihood.
##' 
##' @section Partial residuals:
##' It's also possible to consider partial residuals, in which the contribution of a particular covariate is added back into the model.
##' @return copy of \code{sca} with new layer
##' @examples
##' data(vbetaFA)
##' svbeta <- subset(vbetaFA, ncells==1)
##' svbeta <- svbeta[freq(svbeta)>.4,]
##' window <- function(x1) lapply(assays(x1), function(x2) x2[1:3, 1:6])
##' #total residuals of the response
##' z1 <- zlm(~ Stim.Condition, svbeta, hook=discrete_residuals_hook)
##' window(collectResiduals(z1, svbeta))
##' z2 <- zlm(~ Stim.Condition, svbeta, hook=continuous_residuals_hook)
##' window(collectResiduals(z2, svbeta))
##' z3 <- zlm(~ Stim.Condition, svbeta, hook=combined_residuals_hook)
##' window(collectResiduals(z3, svbeta))
##' #partial residuals
##' colData(svbeta)$ngeneson <- colMeans(assay(svbeta)>0)
##' z5 <- zlm(~ Stim.Condition + ngeneson, svbeta)
##' partialScore(z5, 'Stim.Condition')
collectResiduals <- function(x, sca, newLayerName='Residuals'){
    if(any(newLayerBool <- assayNames(sca) %in% newLayerName)){
        warning('Overwriting layer', newLayerName)
        i <- which(newLayerBool)
    } else{
        i <- length(assays(sca))+1
    }
    mat <- laply(revealHook(x), function(x) x)
    assay(sca, i) <- mat
    assayNames(sca, i) <- newLayerName
    sca
}


#' @describeIn collectResiduals Hook to get the discrete residuals, ie, difference between expected probability of expression and observed
#' @export
discrete_residuals_hook<- function(x){
    if(all(x@fitted["D"])){
        class(x@fitD) <- c("bayesglm","glm","lm")
        Rd <- residuals(x@fitD,type="response")
        Rd
    }
}

#' @describeIn collectResiduals Hook to get the continuous residuals, ie, residuals for conditionally positive observations.  If an observation is zero, it's residual is defined to be zero as well.
#' @export
continuous_residuals_hook<- function(x){
    if(all(x@fitted["C"])){
        class(x@fitC) <- c("glm","lm")
        class(x@fitD) <- c("bayesglm","glm","lm")
        R <- residuals(x@fitC,type="response")
        Rd <- residuals(x@fitD,type="response")
        Rd[names(R)] <- R
                                        #keep zeros as zeros
        Rd[setdiff(names(Rd),names(R))] <- 0
        Rd
    }
}

#' @describeIn collectResiduals Hook to get the combined residuals, ie, Y-E(U)*E(V)
#' @export
combined_residuals_hook<- function(x){
    if(all(x@fitted)){
        class(x@fitC) <- c("glm","lm")
        class(x@fitD) <- c("bayesglm","glm","lm")
        fc <- x@modelMatrix%*%coef(x@fitC)
        fd <- fitted(x@fitD)
        R <- matrix((x@response-fc*fd),nrow=1)
        colnames(R) <- names(residuals(x@fitD))
        R
    }
}

.getQRlm<-function (x, ...) 
{
    if (is.null(r <- x$qr)) 
        stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
    r
}

                                        #bayesglm.influence called from influence.bayesglm
bayesglm.influence <-  function(model, do.coef = do.coef, ...)
{
    wt.res <- weighted.residuals(model)
    e <- na.omit(wt.res)
    if (model$rank == 0) {
        n <- length(wt.res)
        sigma <- sqrt(deviance(model)/df.residual(model))
        res <- list(hat = rep(0, n), coefficients = matrix(0, 
                                                           n, 0), sigma = rep(sigma, n), wt.res = e)
    }
    else {
        e[abs(e) < 100 * .Machine$double.eps * median(abs(e))] <- 0
        mqr <- .getQRlm(model)
        mqr$qr<-mqr$qr[!rownames(mqr$qr)%in%"",]
        n <- as.integer(nrow(mqr$qr))
        if (is.na(n)) 
            stop("invalid model QR matrix")
        if (NROW(e) != n) 
            stop("non-NA residual length does not match cases used in fitting")
        do.coef <- as.logical(do.coef)
        tol <- 10 * .Machine$double.eps
        C_infl <- get("C_influence",getNamespace("stats"))
        res <- .Call(C_infl, mqr, do.coef, e, tol)
        res$wt.res <- e # .Call above is not returning wt.res anymore?
        if (!is.null(model$na.action)) {
            hat <- naresid(model$na.action, res$hat)
            hat[is.na(hat)] <- 0
            res$hat <- hat
            if (do.coef) {
                coefficients <- naresid(model$na.action, res$coefficients)
                coefficients[is.na(coefficients)] <- 0
                res$coefficients <- coefficients
            }
            sigma <- naresid(model$na.action, res$sigma)
            sigma[is.na(sigma)] <- sqrt(deviance(model)/df.residual(model))
            res$sigma <- sigma
        }
    }
    res$wt.res <- naresid(model$na.action, res$wt.res)
    res$hat[res$hat > 1 - 10 * .Machine$double.eps] <- 1
    names(res$hat) <- names(res$sigma) <- names(res$wt.res)
    if (do.coef) {
        rownames(res$coefficients) <- names(res$wt.res)
        colnames(res$coefficients) <- names(coef(model))[!is.na(coef(model))]
    }
    res
}

#' Influence bayesglm object
#'
#' The influence function
#' @importFrom stats influence
#' @param model \code{bayesglm}
#' @param do.coef see \link{influence.glm}
#' @param ... ignored
#' @return see \link{influence.glm}
influence.bayesglm <- function (model, do.coef = TRUE, ...) 
{
    res <- bayesglm.influence(model, do.coef = do.coef, ...)
    pRes <- na.omit(residuals(model, type = "pearson"))[model$prior.weights != 
                                                        0]
    pRes <- naresid(model$na.action, pRes)
    names(res)[names(res) == "wt.res"] <- "dev.res"
    c(res, list(pear.res = pRes))
}

#' rstandard for bayesglm objects.
#'
#' rstandard bayesglm object S3 method
#' @importFrom stats rstandard
#' @param model \code{bayesglm}
#' @param infl see \link{rstandard}
#' @param type see \link{rstandard}
#' @param ... ignored
#' @return \code{numeric} residuals
rstandard.bayesglm <- function (model, infl = influence(model, do.coef = FALSE), type = c("deviance", "pearson"), ...)
{
    type <- match.arg(type)
    res <- switch(type, pearson = infl$pear.res, infl$dev.res)
    res <- res/sqrt(summary(model)$dispersion * (1 - infl$hat))
    res[is.infinite(res)] <- NaN
    res
}

#' @describeIn collectResiduals Standardized deviance residuals hook. Computes the sum of the standardized deviance residuals for the discrete and continuous models (scaled to have unit variance).  If the observation is zero then only the discrete component is used.
#' @export
deviance_residuals_hook<-function (x) 
{
    if (all(x@fitted)) {
        class(x@fitC) <- c("glm", "lm")
        class(x@fitD) <- c("bayesglm", "glm", "lm")
        cont.resid<-rstandard(x@fitC,type="deviance")
        disc.resid<-rstandard(x@fitD,type="deviance")
        cont.resid<-data.table(id=names(x@fitC$y),cont.resid)
        disc.resid<-data.table(id=names(x@fitD$y),disc.resid)
        resid<-merge(cont.resid,disc.resid,by="id",all=TRUE)
        namean <- function(x, y){
            nax <- is.na(x)
            nay <- is.na(y)
            (ifelse(nax, 0, x)+ifelse(nay,0, y))/sqrt((!nax)*1+(!nay)*1)
        }
        resid[, comb:=namean(cont.resid, disc.resid)]
        resid <- resid[data.table(id=rownames(x@modelMatrix)),,on='id']
        return(setNames(resid[,comb], resid[,id]))
    }
}

if(getRversion() >= "2.15.1") globalVariables(c('comb'))

#' @describeIn collectResiduals Hook to return p_hat, the predicted probability of expression.
#' @export
fitted_phat <- function(x){
    if(all(x@fitted)){
        class(x@fitC) <- c("glm","lm")
        class(x@fitD) <- c("bayesglm","glm","lm")
        fd <- fitted(x@fitD)
        fd
    }
}



safeCP <- function(x, y){
    cx <- complexifyNA(x)
    cy <- complexifyNA(y)
    res <- crossprod(cx, cy)
    uncomplexify(res)
}

#' @export
#' @describeIn collectResiduals Compute \eqn{Y_i-E(V_i|X_i, Z_0)E(U|X_i, Z_0)}, where \eqn{Z_0} is a  treatment effect (being left in) and \eqn{X_i} is a nuisance effect (being regressed out).
#' @param effectRegex a regular expression naming columns of the design corresponding to \eqn{Z_0}.
#' Generally these should be the treatment effects of interest.
partialScore <- function(x, effectRegex){
    MMall <- x@LMlike@modelMatrix
    effects <- colnames(MMall) %like% effectRegex
    MM <- MMall[,!effects,drop=FALSE]
    coefD <- coef(x, 'D')[,!effects, drop=FALSE]
    coefC <- coef(x, 'C')[,!effects,drop=FALSE]
    predC <- safeCP(t(coefC), t(MM))
    predD <- safeCP(t(coefD), t(MM))
    res <- assay(x@sca)-predC * invlogit(predD)
}
