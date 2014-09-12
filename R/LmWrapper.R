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
    validObject(object)
    fit(object, silent=silent, ...)
})

setMethod('coef', signature=c(object='LMlike'), function(object, which, singular=TRUE, ...){
    stopifnot(which %in% c('C', 'D'))
    co <- if(which=='C') coef(object@fitC) else coef(object@fitD)
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
    object@modelMatrix <- model.matrix(object@formula, object@design, ...)
    object@fitC <- object@fitD <- numeric(0)
    object@fitted <- c(C=FALSE, D=FALSE)
    object
})

setMethod('model.matrix', signature=c(object='LMlike'), function(object){
    object@modelMatrix
})

setReplaceMethod('model.matrix', signature=c(object='LMlike'), function(object, value){
    object@modelMatrix <- value
    validObject(object)
    object
})





makeChiSqTable <- function(lambda, df, test){
    stopifnot(all(names(lambda) == c('C', 'D')))
    stopifnot(all(names(df) == c('C', 'D')))
    lambdaC <- c(lambda, sum(lambda))
    dfC <- c(df, sum(df))
    tab <- cbind(lambda=lambdaC,
               df=dfC, 'Pr(>Chisq)'=pchisq(lambdaC, df=dfC, lower.tail=FALSE))
    row.names(tab) <- c('cont', 'disc', 'hurdle')
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
## object0 null model (possibly unfitted)
## returns chisqtable
.lrTest <- function(object1, object0){
    l1 <- logLik(object1)
    object0 <- fit(object0)
    l0 <- logLik(object0)
    bothfitted <- object1@fitted & object0@fitted
    dl <- ifelse(bothfitted, -2*(l0-l1), c(0, 0))
    df <- ifelse(bothfitted, dof(object1) - dof(object0), c(0, 0))
    makeChiSqTable(dl, df, drop.terms)
}

setMethod('lrTest', signature=c(object='LMlike', hypothesis='character'), function(object, hypothesis){
    F <- update.formula(object@formula, formula(sprintf(' ~. - %s', hypothesis)))
    U <- update(object, F)
    .lrTest(object, U)
})

setMethod('lrTest', signature=c(object='LMlike', hypothesis='CoefficientHypothesis'), function(object, hypothesis){
    object1 <- object
    mm <- object0@modelMatrix
    object0@modelMatrix <- mm[,setdiff(names(mm), hypothesis), drop=FALSE]

    ## Don't attempt to test any component with a missing coefficient
    missingCoefC <- setdiff(hypothesis, coef(object1, which='C', singular=FALSE))
    missingCoefD <- setdiff(hypothesis, coef(object1, which='D', singular=FALSE))
    if(length(missingCoefC)>0){
        warning('Not testing continuous.', paste(missingCoefC, ','), 'was missing')
        object1@fitted['C'] <- FALSE
}
        if(length(missingCoefD)>0){
        warning('Not testing discrete.', paste(missingCoefD, ','), 'was missing')
        object1@fitted['D'] <- FALSE
}
    .lrTest(object1, object0)
})


setMethod('lrTest', signature=c(object='LMlike', hypothesis='Hypothesis'), function(object, hypothesis){
    ## ## from glmLRT in edgeR
    contrast <- as.matrix(hypothesis@transformed)
    qrc <- qr(contrast)
    ncontrasts <- qrc$rank
    if(ncontrasts==0) stop("contrasts are all zero")
    coef <- 1:ncontrasts
    if(ncontrasts < ncol(contrast)) contrast <- contrast[,qrc$pivot[coef]]
    if(ncontrasts>1) {
        coef.name <- paste("LR test of",ncontrasts,"contrasts")
    } else {
        contrast <- drop(contrast)
        i <- contrast!=0
        coef.name <- paste(paste(contrast[i],coef.names[i],sep="*"),collapse=" ")
    }

    ## rotate design and drop columns
    Dvec <- rep.int(1,nlibs)
    Dvec[coef] <- diag(qrc$qr)[coef]
    ## what is Dvec??
    Q <- qr.Q(qrc,complete=TRUE,Dvec=Dvec)
    design <- model.matrix(object)
    design <- design %*% Q
    design0 <- design[,-coef,drop=FALSE]
    model.matrix(object) <- design
    coefname <- colnames(design[,-coef,drop=FALSE])

    ## Ok, so we rotated the design--now to see if we can estimate all of the columns that we need...
    object <- fit(object)
    ## now we have a CoefficientHypothesis defined in the rotated space
    hypothesis <- new('CoefficientHypothesis', .Data=coefname, transformed=coefname)
    lrTest(object,hypothesis)
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
