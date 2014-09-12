## static constructor of summary components
initSummaries <- function(genes, coefNames){
    p <- length(coefNames)
    ng <- length(genes)
    coefD <- coefC <- matrix(NA, nrow=ng, ncol=p, dimnames=list(primerid=genes, coef=coefNames))
    vcovC <- vcovD <- array(NA, dim=c(ng, p, p), dimnames=list(primerid=genes, coefNames, coefNames))
dispersionMLEC <- df.residD <- df.residC <- df.nullC <- df.nullD <- devianceC <- devianceD <- setNames(rep(NA, ng), genes)
    as.environment(list(coefC=coefC, vcovC=vcovC, df.residC=df.residC, df.nullC=df.nullC, devianceC=devianceC, dispersionMLEC=dispersionMLEC, coefD=coefD, vcovD=vcovD, df.residD=df.residD, df.nullD=df.nullD, devianceD=devianceD))
}
