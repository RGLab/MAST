lrtest <- function(w.x, w.y, x, y){
    ## w.x, w.y vectors of zeros/ones for expressed or not in each group
    ## x, y vectors of the positive observations (must be of length sum(w.x) and sum(w.y))

    e.x <- sum(w.x)
    e.y <-  sum(w.y)
    n.x <-  length(w.x)
    n.y <-  length(w.y)
    stopifnot(e.x == length(x) && e.y == length(y))


    p.0 <- (e.x+e.y)/(n.x + n.y)
    p.x <- e.x/n.x
    p.y <- e.y/n.y

    m0 <-  (sum(x)+sum(y))/(e.x+e.y)
    mu.x <-  mean(x)
    mu.y <-  mean(y)

    Tstar <-  1+e.x*e.y/(e.x+e.y)* (mu.x - mu.y)^2/(sum((mu.x - x)^2) + sum((mu.y-y)^2))

    if(!is.finite(Tstar)){
        Tstar <- 1
    }

    binom <- logProd(e.x, p.0/p.x) +
        logProd(e.y, p.0/p.y) +
        logProd(n.x-e.x, (1-p.0)/(1-p.x)) +
        logProd(n.y-e.y, (1-p.0)/(1-p.y))
    binomsign <- (p.y>p.x)*2 -1

    norm <- -(e.x+e.y)/2 * log(Tstar)
    normsign <- (mu.y>mu.x)*2-1

    logLR <- binom+norm

    maxsign <- c(binomsign, normsign)[which.min(c(binom, norm))]
    resultvec <- c(-2*binom, binomsign, pchisq(-2*binom, 1, lower.tail=FALSE),
                   -2*norm, normsign, pchisq(-2*norm, 1, lower.tail=FALSE),
                   -2*logLR, maxsign, pchisq(-2*logLR, 2, lower.tail=FALSE))
    result <- matrix(resultvec, nrow=3, ncol=3, dimnames=list(metric=c('lrstat', 'direction', 'p.value'), component=c('binom', 'norm', 'comb')))
}

logProd <- function(prod, logand){
    ifelse(prod==0, 0, prod*log(logand))
}


##' @rdname LRT
##' @param sca A \code{SingleCellAssay} class object
##' @param comparison A \code{character} specifying the factor for comparison
##' @param referent A \code{character} specifying the reference level of \code{comparison}.
##' @param groups A optional \code{character} specifying a variable on which to stratify the test.  For each level of \code{groups}, there will be a separate likelihood ratio test.
##' @param returnall A \code{logical} specifying if additional columns should be returned with information about the different components of the test.
##' @export
setMethod("LRT",signature=c("SingleCellAssay","character"),function(sca,comparison,referent=NULL,groups=NULL,returnall=FALSE){
    if(missing(groups))
        groups<-NULL
    if(missing(referent))
        referent <- NULL
    lrt(sca,comparison,referent,groups=groups,returnall=returnall)
})

lrt <- function(sca, comparison, referent=NULL, groups=NULL, returnall=TRUE){
    if (missing(comparison) || !checkGroups(sca, comparison))
        stop("'comparison' missing or incorrect")
    ## what happens if comparision has length >1?

    if(!is.null(groups)){
        checkGroups(sca, groups)
        ## we should check what happens if comparison has a different number of levels
        scL <- split(sca, groups)
        lapp <- lapply(scL, lrt, comparison=comparison, referent=referent, groups=NULL, returnall=returnall)
        ## fix
        retme<-do.call(rbind, lapp)
        nr<-lapply(lapp,nrow)
        nms<-names(lapp)
        retme<-rename(cbind(retme,groups=factor(do.call(c,lapply(seq_along(nr),function(i)rep(nms[i],nr[i]))))),c(groups=groups))
        return(retme)
    }

                                        #getMapping returns a list.. code expects a vector
    probeid <- 'primerid' 
    measure <- 'value'

    scadt <- as(sca, 'data.table')
    phenocol <- scadt[[comparison]]
    if(is.null(referent)){
        pheno.order <- factor(phenocol)
    } else{
        pheno.order <- factor(phenocol)
        pheno.order <- relevel(pheno.order, ref=referent)
    }
    nlev <- nlevels(pheno.order)

    ssca <- split(cbind(scadt[, c(measure, comparison), with=FALSE], pheno.order), scadt[,probeid,with=FALSE], drop=TRUE) #drop=TRUE: seems like the more reasonable default if probeid is a factor and unused levels are present (after subsetting, for example)
    lrout <- vapply(ssca, FUN.VALUE=array(0, dim=c(nlev-1, 3, 4)), FUN=function(x){
        res <- array(NA, dim=c(nlev-1, 3, 4))
        phenosplit <- split(x[[measure]], x$pheno.order, drop=FALSE)
        unstim <- phenosplit[[1]]
        if(any(is.na(unstim))){
            warning('dropping NA measurements')
            unstim <- unstim[!is.na(unstim)]
        }
        w.x <- (unstim>0)*1
        x <- unstim[w.x==1]

        for(i in seq(from=2, to=nlev)){
            stim <- phenosplit[[i]]
            if(any(is.na(stim))){
                warning('dropping NA measurements')
                stim <- stim[!is.na(stim)]
            }
            if (length(stim)==0){
                res[i-1,,] <- NA
                lrtmp <- lrtest(1, 1, 1, 1)        #needed to fill out dimnames of res
                                        #in case all groups had zero measurements
            } else{
                w.y <- (stim>0)*1
                y <- stim[w.y==1]
                lrtmp <- lrtest(w.x, w.y, x, y)
                res[i-1,,1:3] <- lrtmp
                tt <- t.test(2^unstim-1, 2^stim-1, var.equal=TRUE)
                res[i-1,1,4] <- tt$stat
                res[i-1,2,4] <- sign(tt$stat)
                res[i-1,3,4] <- tt$p.value
            }
        }
        dn <- dimnames(lrtmp)
        dn$component <- c(dn$component, 'zeroes')
        dimnames(res) <- c(list(geneid=names(phenosplit)[-1]), dn)
        res
    })

    m <- reshape2::melt(lrout)
    m <- rename(m, c('Var1'=comparison, 'Var2'='metric', 'Var3'='test.type', 'Var4'=probeid))
    if(returnall){
        return(m)
    }
    retme<-subset(m, test.type=='comb')
    return(dcast(rename(retme,c(metric="variable")), formula=...~variable))
}


##' Plot a likelihood ratio test object
##'
##' Constructs a forest-like plot of signed log10 p-values, possibly adjusted for multiple comparisons
##' \code{adjust} can be one of  "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
##' @param lr output from lrtest, with returnall=FALSE
##' @param adjust \code{character}, passed along to \code{p.adjust}, see below
##' @param thres \code{numeric} genes with adjusted pvalues above this value are not depicted
##' @param trunc \code{numeric} p values below this value are truncated at this value 
##' @param groups \code{character} grouping value.  If provided, must match groups argument passed to lrtest.  Plots done separately for each group.
##' @return Constructs a dotplot
##' @author andrew
##' @export
plotlrt <- function(lr, adjust='fdr', thres=.1, trunc=1e-6, groups=NULL){
    lr$adj <- p.adjust(lr$p.value, method=adjust)
    posgene <- suppressMessages(reshape2::dcast(lr[, c('gene', 'adj')], gene ~ ., fun.aggregate=function(x) any(x<thres)))
    posgene <- posgene[posgene[,2],]
    pvalue <- with(lr, pmin(p.value, trunc))
    if(length(posgene)>0){
        lattice::dotplot(gene ~ -log10(pvalue)*direction, lr, auto.key=TRUE, subset=gene %in% posgene$gene)
    } else{
        warning("No significant genes")
    }

}
