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
