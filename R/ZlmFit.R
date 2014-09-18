## static constructor of summary components
initSummaries <- function(genes, coefNames){
    p <- length(coefNames)
    ng <- length(genes)
    coefD <- coefC <- matrix(NA, nrow=ng, ncol=p, dimnames=list(primerid=genes, coef=coefNames))
    vcovC <- vcovD <- array(NA, dim=c(ng, p, p), dimnames=list(primerid=genes, coefNames, coefNames))
dispersionMLEC <- df.residD <- df.residC <- df.nullC <- df.nullD <- devianceC <- devianceD <- setNames(rep(NA, ng), genes)
    as.environment(list(coefC=coefC, vcovC=vcovC, df.residC=df.residC, df.nullC=df.nullC, devianceC=devianceC, dispersionMLEC=dispersionMLEC, coefD=coefD, vcovD=vcovD, df.residD=df.residD, df.nullD=df.nullD, devianceD=devianceD))
}

collectSummaries <- function(listOfSummaries){
    summaries <- list()
    summaries[['coefC']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'coefC'))
    summaries[['vcovC']] <- do.call(abind, lapply(listOfSummaries, '[[', 'vcovC'))
    summaries[['df.residC']] <- unlist(lapply(listOfSummaries, '[[', 'df.residC'))
    summaries[['df.nullC']] <- unlist(lapply(listOfSummaries, '[[', 'df.nullC'))
    summaries[['devianceC']] <- unlist(lapply(listOfSummaries, '[[', 'devianceC'))
    summaries[['dispersionMLEC']] <- unlist(lapply(listOfSummaries, '[[', 'devianceC'))

    summaries[['coefD']] <- do.call(rbind, lapply(listOfSummaries, '[[', 'coefD'))
    summaries[['vcovD']] <- do.call(abind, lapply(listOfSummaries, '[[', 'vcovD'), along=-1)
    summaries[['df.residD']] <- unlist(lapply(listOfSummaries, '[[', 'df.residD'))
    summaries[['df.nullD']] <- unlist(lapply(listOfSummaries, '[[', 'df.nullD'))
    summaries[['devianceD']] <- unlist(lapply(listOfSummaries, '[[', 'devianceD'))

    summaries
}

setMethod('lrTest', signature=c(object='ZlmFit', hypothesis='character'), function(object, hypothesis){


})

setMethod('lrTest', signature=c(object='ZlmFit', hypothesis='Hypothesis'), function(object, hypothesis){
    ## initialize junk
    testNames <- makeChiSqTable(c(0, 0), c(1, 1), '')
    vcovNames <- coefNames <- colnames(MM)

        for(h in seq_len(nhypo)){
        ltests[[h]] <- array(0, dim=c(ng, nrow(testNames), ncol(testNames)), dimnames=list(primerid=genes, test.type=row.names(testNames), metric=colnames(testNames)))
        ltests[[h]][,,'Pr(>Chisq)'] <- 1
}


})
