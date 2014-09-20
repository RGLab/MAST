## Methods for LMlike
setMethod('show',  signature=c(object='LMlike'), function(object){
    if(all(object@fitted)){
        cat('Fitted Continuous and discrete')
    } else if(object@fitted['C']){
        cat('Fitted Continuous')
    }else if(object@fitted['D']){
        cat('Fitted Discrete')
    }else {
        cat('Unfitted')
    }
    cat(':', as.character(object@formula), '\n')
    cat(class(object), ':', nrow(object@design), ' cases\n', sep='')
})

.fit <- function(object){
    frame <- sys.frame(-1)
    positive <- object@response>0
    object@fitted <- c(C=FALSE, D=FALSE)
    assign('pos', positive, pos=frame)
    assign('object', object, pos=frame)
    
    return(any(positive))
}

setMethod('fit', signature=c(object='LMlike', response='vector'), function(object, response, silent=TRUE, fitArgsC=list(), fitArgsD=list(), ...){
    object@response <- response
    object@fitArgsC <- fitArgsC
    object@fitArgsD <- fitArgsD
    ## validObject(object)
    fit(object, silent=silent, ...)
})


setMethod('coef', signature=c(object='LMlike'), function(object, which, singular=TRUE, ...){
    stopifnot(which %in% c('C', 'D'))
    co <- object@defaultCoef
    if(which=='C' & object@fitted['C']){
        co <- coef(object@fitC)}
    else if(object@fitted['D']){
        co <- coef(object@fitD)
    }
    if(!singular) co <- co[!is.na(co)]
    co
})



setMethod('summary', signature=c(object='LMlike'), function(object){
    print('===========Discrete============')
    print(coef(object, 'D'))
    print('==========Continuous===========')
    print(coef(object, 'C'))
})

setMethod('update', signature=c(object='LMlike'), function(object, formula., ...){
    object@formula <- update.formula(object@formula, formula.)
    model.matrix(object) <- model.matrix(object@formula, object@design, ...)
    object@fitC <- object@fitD <- numeric(0)
    object@fitted <- c(C=FALSE, D=FALSE)
    object
})

setMethod('model.matrix', signature=c(object='LMlike'), function(object){
    object@modelMatrix
})

setReplaceMethod('model.matrix', signature=c(object='LMlike'), function(object, value){
    qrm <- qr(value)
    est <- qrm$pivot[seq_len(qrm$rank)]
    if(length(est)<ncol(value)) warning('Coefficients ', paste(colnames(value)[setdiff(qrm$pivot, est)], collapse=', '), ' are never estimible and will be dropped.')
    MM <- value[,est, drop=FALSE]
    object@modelMatrix <- MM
    object@defaultCoef <- as.numeric(setNames(rep(NA, ncol(MM)), colnames(MM)))
    validObject(object)
    object
})





makeChiSqTable <- function(lambda, df, test){
    ## either data.frames or vectors
    stopifnot(all(names(lambda) == c('C', 'D')))
    stopifnot(all(names(df) == c('C', 'D')))
    ## which functions will we use for vectors vs data.frames
    if(inherits(lambda, 'data.frame')){
        Combine <- cbind
        Sum <- rowSums
        Glue <- function(...) abind(..., rev.along=0)
    }else{
        Combine <- c
        Sum <- sum
        Glue <- cbind
    }
    lambdaC <- setNames(Combine(lambda, Sum(lambda)), c('cont', 'disc', 'hurdle'))
    dfC <- setNames(Combine(df, Sum(df)), c('cont', 'disc', 'hurdle'))

    pchi <- pchisq(as.matrix(lambdaC), df=as.matrix(dfC), lower.tail=FALSE)
    tab <- Glue(lambda=lambdaC,
               df=dfC, 'Pr(>Chisq)'=pchi)
    structure(tab, test=test)
}

setMethod('waldTest', signature=c(object='LMlike', hypothesis='character'), function(object, hypothesis){
    if(object@fitted['C']){
            C <- car::linearHypothesis.default(object@fitC, hypothesis.matrix=hypothesis, test='Chisq', vcov.=vcov(object, which='C'), coef.=coef(object, which='C', singular=FALSE), singular.ok=TRUE)[2,c('Df', 'Chisq'), drop=TRUE]
    }else{
        C <- list(Chisq=0, Df=0)
    }
    if(object@fitted['D']){
        D <- car::linearHypothesis.default(object@fitD, hypothesis.matrix=hypothesis, test='Chisq', vcov.=vcov(object, which='D'), coef.=coef(object, which='D', singular=FALSE), singular.ok=TRUE)[2,c('Df', 'Chisq'), drop=TRUE]
    }else{
        D <- list(Chisq=0, Df=0)
    }
    makeChiSqTable(c(C[['Chisq']], D[['Chisq']]), c(C[['Df']],D[['Df']]),hypothesis)
})


## object1 full model (fitted)
## newMM is new model to be tested against
## returns chisqtable
.lrTest <- function(object1, newMM){        
    l1 <- logLik(object1)
    object0 <- object1
    model.matrix(object0) <- newMM
    object0 <- fit(object0)
    l0 <- logLik(object0)

    ## don't test when all coefficients are aliased in large model
    ## ...or maybe it's fine, though conservative
    ## testName <- setdiff(colnames(model.matrix(object1)), colnames(newMM))
    ## missingCoefC <- all(is.na(coef(object, which='C', singular=TRUE))[testName])
    ## missingCoefD <- all(is.na(coef(object, which='D', singular=TRUE))[testName])
    bothfitted <- object1@fitted & object0@fitted
    dl <- ifelse(bothfitted, -2*(l0-l1), c(0, 0))
    df <- ifelse(bothfitted, dof(object1) - dof(object0), c(0, 0))
    drop.terms <- setdiff(colnames(model.matrix(object1)), colnames(model.matrix(object0)))
    makeChiSqTable(dl, df, drop.terms)
}

setMethod('lrTest', signature=c(object='LMlike', hypothesis='character'), function(object, hypothesis){
    F <- update.formula(object@formula, formula(sprintf(' ~. - %s', hypothesis)))
    U <- update(object, F)
    .lrTest(object, U@modelMatrix)
})

.rotateMM <- function(object, contrast){
      ## from glmLRT in edgeR
    qrc <- qr(contrast)
    ncontrasts <- qrc$rank
    if(ncontrasts==0) stop("contrasts are all zero")
    testIdx <- 1:ncontrasts
    if(ncontrasts < ncol(contrast)) contrast <- contrast[,qrc$pivot[testIdx]]
    ## if(ncontrasts>1) {
    ##     coef.name <- paste("LR test of",ncontrasts,"contrasts")
    ## } else {
    ##     contrast <- drop(contrast)
    ##     i <- contrast!=0
    ##     coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
    ## }

    ## rotate design and drop columns
    design <- model.matrix(object)
    Dvec <- rep.int(1,ncol(design))
    Dvec[testIdx] <- diag(qrc$qr)[testIdx]
    ## what is Dvec??
    Q <- qr.Q(qrc,complete=TRUE,Dvec=Dvec)
    design <- design %*% Q
    ## Ok, so we rotated the design, and now the coefficients are arbitrary
    ## But might be needed for some subclasses (eg glmer)
    colnames(design) <- paste('X', seq_len(ncol(design)), sep='')
    structure(design, testIdx=testIdx)
}

setMethod('lrTest', c(object='LMlike', hypothesis='CoefficientHypothesis'), function(object, hypothesis){
    testIdx <- hypothesis@transformed
    .lrTest(object, model.matrix(object)[,-testIdx,drop=FALSE])
})


setMethod('lrTest', signature=c(object='LMlike', hypothesis='Hypothesis'), function(object, hypothesis){
    lrTest(object, hypothesis@transformed)    
})

setMethod('lrTest', signature=c(object='LMlike', hypothesis='matrix'), function(object, hypothesis){
    MM <- .rotateMM(object,  hypothesis)
    testIdx <- attr(MM, 'testIdx')
    .lrTest(object, MM[,-testIdx,drop=FALSE])
    ## drop tested contrast
    })

setMethod('residuals', signature=c(object='LMlike'), function(object, type='response', which, ...){
    which <- match.arg(which, c('Discrete', 'Continuous', 'Marginal'))
    RD <- residuals(object@fitD, type=type)
    RC <- residuals(object@fitC, type=type)
    if(which=='Discrete') return(RD)
    if(which=='Continuous') return(RC)
    if(which=='Marginal'){
        if(type != 'response') warning("Marginal residuals probably don't make sense unless predicting on the response scale")
        ## Zero inflated residuals
        PC <-  predict(object@fitC, newx=object@design, type=type)
        PD <- predict(object@fitD, newx=object@design, type=type)
        return(object@response-PC*PD)
    }
})
