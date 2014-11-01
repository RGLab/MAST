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

setMethod('fit', signature=c(object='LMlike', response='vector'), function(object, response, silent=TRUE, fitArgsC=list(), fitArgsD=list(), quick=FALSE, ...){
    object@response <- response
    object@fitArgsC <- fitArgsC
    object@fitArgsD <- fitArgsD
    object@fitted <- c(C=FALSE, D=FALSE)
    object@fitC <- NULL
    object@fitD <- NULL
    if(!quick) validObject(object)      #save time in inner loop in zlm.SingleCellAssay
    fit(object, silent=silent, start=start, ...)
})


setMethod('coef', signature=c(object='LMlike'), function(object, which, singular=TRUE, ...){
    stopifnot(which %in% c('C', 'D'))
    co <- object@defaultCoef
    if(which=='C' & object@fitted['C']){
        co <- coef(object@fitC)
    }

    if(which=='D' & object@fitted['D']){
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
    object@defaultCoef <- setNames(as.numeric(rep(NA, ncol(MM))), colnames(MM))
    object@defaultVcov <- object@defaultCoef %o% object@defaultCoef
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
        Flatten <- function(x) x
    }else{
        Combine <- c
        Sum <- sum
        Glue <- cbind
        Flatten <- as.vector
    }
    lambdaC <- setNames(Combine(lambda, Sum(lambda)), c('cont', 'disc', 'hurdle'))
    dfC <- setNames(Combine(df, Sum(df)), c('cont', 'disc', 'hurdle'))
    pchi <- Flatten(pchisq(as.matrix(lambdaC), df=as.matrix(dfC), lower.tail=FALSE))
    tab <- Glue(lambda=lambdaC,
               df=dfC, 'Pr(>Chisq)'=pchi)
    structure(tab, test=test)
}

## this allows NAs to propagate via matrix multiplication, but still be killed when multiplied by zero
complexifyNA <- function(x){
    x[is.na(x)] <- 0+1i
    x
}
uncomplexify <- function(x){
    x[abs(Im(x))>.Machine$double.eps] <- NA
    nx <- as.numeric(x)
    dim(nx) <- dim(x)
    nx
}



## Takes coefficients (vector), contrast (matrix) and vcov (matrix), gets linear combination(s) defined by contrast matrix
## Get squared length of linear combinations
## This will be called by both ZlmFit and LmFit
.waldTest <- function(coefC, coefD, vcC, vcD, cm, fitted){
    ## linear combination of coefficients
    ## complexify NA allows NAs to propagate algebraically rather than symbolically (ie 0*NA = 0)
    contrC <- crossprod(cm, complexifyNA(coefC))
    ## covariance matrix of linear combination
    contrCovC <- crossprod(cm, complexifyNA(vcC) %*% cm)
    contrD <- crossprod(cm, complexifyNA(coefD))
    contrCovD <- crossprod(cm, complexifyNA(vcD) %*% cm)

    ## Don't consider comparisons when we failed to fit
    dof <- c(C=nrow(contrC), D=nrow(contrD))*(fitted*1)
    contrC <- uncomplexify(contrC)
    contrD <- uncomplexify(contrD)
    ## squared length of coefficient linear combination, using its covariance
    lambdaC <- tryCatch(contrC %*% solve(matrix(uncomplexify(contrCovC),ncol=ncol(contrCovC)), contrC), error=function(e){
        reraise(e, convertToWarning=TRUE)
        return(0)
    })

    lambdaD <- tryCatch(contrD%*%solve(matrix(uncomplexify(contrCovD),ncol=ncol(contrCovD)), contrD), error=function(e){
        reraise(e, convertToWarning=TRUE)
        return(0)
    })

    makeChiSqTable(c(C=lambdaC, D=lambdaD)*(fitted*1), dof,cm)
}

.makeContrastMatrixFromCoefficientHypothesis <- function(hypothesis, coefnames){
    h <- generateHypothesis(hypothesis, coefnames)
    testIdx <- h@transformed
    cm <- matrix(0, nrow=length(testIdx), ncol=length(coefnames), dimnames=list(contrast=h, coefnames))
    cm[cbind(seq_along(testIdx), testIdx)] <- 1
    t(cm)
}

setMethod('waldTest', signature=c(object='LMlike', hypothesis='CoefficientHypothesis'), function(object, hypothesis){
    cm <- .makeContrastMatrixFromCoefficientHypothesis(hypothesis,names(object@defaultCoef))
    waldTest(object, cm)
})


setMethod('waldTest', signature=c(object='LMlike', hypothesis='matrix'), function(object, hypothesis){
    .waldTest(coef(object, 'C', singular=TRUE),
                coef(object, 'D', singular=TRUE),
                vcov(object, 'C', singular=TRUE),
                vcov(object, 'D', singular=TRUE),
              hypothesis,
              object@fitted)
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
