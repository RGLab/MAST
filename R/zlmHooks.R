#' hHook to get the continuous residuals
#'
#' @param x the ZLMFit
#' @export
discrete_residuals_hook<- function(x){
    if(all(x@fitted["D"])){
        class(x@fitD) <- c("bayesglm","glm","lm")
        Rd <- residuals(x@fitD,type="response")
        Rd
    }
}

#' Hook to get the discrete residuals
#'
#' @param x the ZLMFit
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

#' Hook to get the combined residuals
#' @param x the ZLMFit
#'
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
    C_infl<-get("C_influence",getNamespace("stats"))
    res <- .Call(C_infl, mqr, do.coef, e, tol)
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

#' influence bayesglm object S3 method
#' @importFrom stats influence
#' @name influence.bayesglm
#' @title influence for bayesglm objects.
#' @export
R.methodsS3:::setMethodS3("influence","bayesglm",definition=function (model, do.coef = TRUE, ...) 
{
  res <- bayesglm.influence(model, do.coef = do.coef, ...)
  pRes <- na.omit(residuals(model, type = "pearson"))[model$prior.weights != 
                                                        0]
  pRes <- naresid(model$na.action, pRes)
  names(res)[names(res) == "wt.res"] <- "dev.res"
  c(res, list(pear.res = pRes))
})

#' rstandard bayesglm object S3 method
#' @importFrom stats rstandard
#' @name rstandard.bayesglm
#' @title rstandard for bayesglm objects.
#' @export
R.methodsS3:::setMethodS3("rstandard","bayesglm",definition=function (model, infl = influence(model, do.coef = FALSE), type = c("deviance", 
                                                                                                                                "pearson"), ...) 
{
  type <- match.arg(type)
  res <- switch(type, pearson = infl$pear.res, infl$dev.res)
  res <- res/sqrt(summary(model)$dispersion * (1 - infl$hat))
  res[is.infinite(res)] <- NaN
  res
})

#' Standardized deviance residuals hook
#' 
#' Computes the average of the standardized deviance residuals for the discrete and continuous models
#' @param x the ZLMFit
#' @export
deviance_residuals_hook<-function (x) 
{
  if (all(x@fitted)) {
    class(x@fitC) <- c("glm", "lm")
    class(x@fitD) <- c("bayesglm", "glm", "lm")
    cont.resid<-rstandard(x@fitC,type="deviance")
    disc.resid<-rstandard(x@fitD,type="deviance")
    cont.resid<-data.frame(id=names(x@fitC$y),cont.resid)
    disc.resid<-data.frame(id=names(x@fitD$y),disc.resid)
    resid<-merge(cont.resid,disc.resid,by="id",all=TRUE)
    resid<-data.frame(data.table(melt(resid))[,list(resid=mean(value,na.rm=TRUE)),id])
    rownames(resid)<-resid[,"id"]
    resid<-resid[,-1,drop=FALSE]
    resid<-resid[rownames(x@modelMatrix),] #ensure consistent ordering
    resid
  }
}

#' Used to compute module "scores"
#'
#' Hook to compute \eqn{ y_i-E(V_i)-\beta_{ngeneson} ng_i} for \eqn{y_i>0} and \eqn{E(U_i)(E(V_i)-\beta_{ngeneson} \times ng_i)} for \eqn{y_i=0}
#' @param x the ZLMFit
#' @export
score_hook <- function(x){
    if(all(x@fitted)){
        class(x@fitC) <- c("glm","lm")
        class(x@fitD) <- c("bayesglm","glm","lm")
        fd <- fitted(x@fitD)
        wh <- colnames(x@modelMatrix)%like%"cngeneson"
        wh2 <- names(coef(x@fitC))%like%"cngeneson"
        correction <- x@modelMatrix[,wh,drop=FALSE]%*%coef(x@fitC)[wh2,drop=FALSE]
        fc <- x@modelMatrix%*%coef(x@fitC)-correction
        score <- fd*fc
        score[x@response>0] <- ((x@response-correction))[x@response>0]
        score <- matrix(score,nrow=1)
        colnames(score) <- names(fd)
        score
    }
}

#' Hook to return p_hat from the model
#'
#' Used to weight the unobserved values in a module score calculation analogous to what is done in Shalek et.al.
#' @param x the ZLMFit
#' @export
shalek_weights <- function(x){
    if(all(x@fitted)){
        class(x@fitC) <- c("glm","lm")
        class(x@fitD) <- c("bayesglm","glm","lm")
        fd <- fitted(x@fitD)
        fd
    }
}



#' Hook to compute the combined residuals with treatment effects added back.
#'
#' @param x the ZLMFit
#' @export
score_2_hook <- function(x){
    if(all(x@fitted)){
        class(x@fitC) <- c("glm","lm")
        class(x@fitD) <- c("bayesglm","glm","lm")
        wh <- !colnames(x@modelMatrix)%like%"cngeneson"
        wh2 <- !names(coef(x@fitC))%like%"cngeneson"
        correction <- x@modelMatrix[,wh,drop=FALSE]%*%coef(x@fitC)[wh2,drop=FALSE] #treatment effect
        fc <- x@modelMatrix%*%coef(x@fitC)-correction #fitted minus treatment effect
        fd <- fitted(x@fitD) #discrete fitted effect
        R <- matrix((x@response-fc*fd),nrow=1) #residuals corrected for ngeneson in the continuous part
        colnames(R) <- names(residuals(x@fitD))
        R
    }
}

## # as with score_2 but we correct the discrete part for ngeneson
## score_2_discrete_hook <- function(x){
##     if(all(x@fitted)){
##         class(x@fitC) <- c("glm","lm")
##         class(x@fitD) <- c("bayesglm","glm","lm")
##         wh <- !colnames(x@modelMatrix)%like%"cngeneson"
##         wh2 <- !names(coef(x@fitD))%like%"cngeneson"
##         correction <- x@modelMatrix[,wh,drop=FALSE]%*%coef(x@fitD)[wh2,drop=FALSE]
##         fd <- x@modelMatrix%*%coef(x@fitD)-correction #remove only ngeneson effect
##         R <- matrix(((x@response>0)-invlogit(fd)),nrow=1)
##         colnames(R) <- names(residuals(x@fitD))
##         R
##     }
## }

## ngeneson_hook_D<- function(x){
##     if(all(x@fitted)){
##         class(x@fitC) <- c("glm","lm")
##         class(x@fitD) <- c("bayesglm","glm","lm")
##         wh <- colnames(x@modelMatrix)%like%"cngeneson"
##         wh2 <- names(coef(x@fitD))%like%"cngeneson"
##         ngd <- (x@modelMatrix[,wh,drop=FALSE]%*%coef(x@fitD)[wh2,drop=FALSE])
##         return(t(ngd))
##     }
## }

## ngeneson_hook_C<- function(x){
##     if(all(x@fitted)){
##         class(x@fitC) <- c("glm","lm")
##         class(x@fitD) <- c("bayesglm","glm","lm")
##         wh <- colnames(x@modelMatrix)%like%"cngeneson"
##         wh2 <- names(coef(x@fitC))%like%"cngeneson"
##         ngc <- x@modelMatrix[,wh,drop=FALSE]%*%coef(x@fitC)[wh2,drop=FALSE]
##         return(t(ngc))
##     }
## }
