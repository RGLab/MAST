## Set up standard filter set
sigmaProportion = c(3, 5, 7, 9)
filtercrit = rbind(expand.grid(filter=T,
  nOutlier = c(1, 2),
  sigmaContinuous=c(3, 5, 7, 9),
  sigmaProportion=sigmaProportion,
  sigmaSum = NA,
  K=1.48
  ),
  expand.grid(filter=T,
              nOutlier = 0,
              sigmaContinuous = 0,
              sigmaProportion=sigmaProportion,
              sigmaSum = c(30, 35, 40, 45),
              K=1.48))


## Run filters, write summary output
tryFilters <- function(mysc, myhc_filt, filtercrit, groups=NULL, doR2=FALSE, ...){
  filtercritout = cbind(filtercrit, wss=Inf, ss=Inf, rc=0)
   nref <- nrow(mysc)

  ## Internal function called by tryFilters
.tryFilters <- function(i){
  thisfilt = filter(mysc, filt_control=filtercrit[i,])
  x = getConcordance(myhc_filt, thisfilt, groups)
  r2 <- NA
  ncells <- nref - nrow(thisfilt)
  if(doR2) r2 <- mean(getR2(thisfilt, ...), na.rm=TRUE)
  list(concord=x, r2=r2, ncells=ncells)
}

  if(any(grepl("parallel",loadedNamespaces()))){
   app <- mclapply(1:nrow(filtercrit), .tryFilters)
  } else{
    app <- lapply(1:nrow(filtercrit), .tryFilters)
}

  lo <- lapply(app, function(x) x$concord)
  filtercritout$r2 <- sapply(app, function(x) x$r2)
  filtercritout$ncells <- sapply(app, function(x) x$ncells)
  concord <- getConcordance(myhc_filt, mysc, groups)
  for(i in  seq_len(nrow(filtercritout))){
    filtercritout[i, 'wss'] <- getwss(lo[[i]], concord$nexp.comp)
    filtercritout[i, 'ss'] <- getss(lo[[i]])
    filtercritout[i, 'rc'] <- getrc(lo[[i]])
  }


  list(filtercritout, lo)

}


## R2 evaluation
FREQ_CUTOFF <- .7
getR2 <- function(sc, hkeepers=c('GAPDH', 'POLR2A'), what=function(x) summary(x)$r.square){
  therest <- setdiff(names(which(freq(sc)>FREQ_CUTOFF)), hkeepers)
  out <- rep(NA, length.out=length(therest))
  names(out) <- therest
  theexp <- exprs(sc)
  for(g in therest){
    form <-  sprintf("%s ~ %s + 1", g, paste(hkeepers, collapse=" + "))
    if(all(is.na(theexp[,g]))) next
    fit <- lm(form, as.data.frame(theexp))
    out[g] <- what(fit)
  }
  out
}

## Plotting helpers
as.orderedReverse <- function(x){lev <- sort(unique(x), decreasing=T); factor(x, levels=lev, ordered=T)}

getReferents <- function(fc.sel, baseline, lo.sel){
  nlev <- sqrt(nrow(fc.sel))
  sigmap <- sort(unique(fc.sel$sigmaProportion))
  sigmac <- sort(unique(fc.sel$sigmaContinuous), decreasing=TRUE)

  matc <- matco <- matrix(sigmac, nrow=nlev, ncol=nlev, byrow=TRUE) #it's easiest to assign the referents if we put everyone into matrices
  matp <- matpo <- matrix(sigmap, nrow=nlev, ncol=nlev) #then move things around 'ad-hoc'
  matc[,2:nlev] <- matc[, -nlev]
  matp[-nlev,1] <- matp[2:nlev,1]
  matp[nlev, 1] <- 'baseline'
  matc[nlev, 1] <- 'baseline'

                                        #  contKey <- fc.sel$sigmaContinuous
                                        #  contKey[contKey == max(sigmac)] = 'baseline' #set leftmost column to baseline sigmaCont
                                        #  contKey[!contKey=='baseline'] = sigmac[match(contKey[!contKey=='baseline'], sigmac)+1] #set other columns to neighbors to the left
                                        #  propKey <- ifelse(fc.sel$sigmaContinuous==max(sigmac), 'baseline', fc.sel$sigmaProportion)
  lo2<-c(lo.sel, list(baseline))
  fc.sel$sigmaContinuous <- as.character(fc.sel$sigmaC)
  fc.sel$sigmaProportion <- as.character(fc.sel$sigmaP)
  fc.sel2<- rbind.fill(fc.sel, data.frame(sigmaContinuous='baseline', sigmaProportion='baseline', stringsAsFactors=F))
  matmap <- with(fc.sel2, match(paste(matpo, matco), paste(sigmaProportion, sigmaContinuous))) #(matpo, matco) -> fc.sel2
  matmap2 <- order(matmap) # fc.sel2 -> (matmp, matco)
  sel <- with(fc.sel2, match(paste(matp, matc), paste(sigmaProportion, sigmaContinuous)))[matmap2]
  ref.sel <- do.call(rbind, lo2[sel])
  ref.sel
}

## tryFilterOut: list of length two from tryFilter
## concord: output from getConcordance
concatFilterOutput <- function(tryFilterOut, concord){
  lo <- tryFilterOut[[2]]
  filtercrit <- tryFilterOut[[1]]
  filtercrit.sel <- filtercrit[rep(1:length(lo), each=nrow(lo[[1]])),]
  filtercrit.sel$sigmaContinuous <- as.orderedReverse(filtercrit.sel$sigmaC)
  filtercrit.sel$sigmaProportion <- as.orderedReverse(filtercrit.sel$sigmaP)

  concord.sel <- do.call(rbind, lo)[, c('et.ref', 'et.comp')]
  concord.sel2 <- cbind(concord.sel, filtercrit.sel)

  ref.sel <- getReferents(filtercrit, concord, lo) #line up reference positions
  list(concord.sel2, ref.sel)
}

makePlots <- function(tf, concat, mysc){
  tf[[1]]$sigmaContinuous <- as.orderedReverse(tf[[1]]$sigmaContinuous)
    tf[[1]]$sigmaProportion <- as.orderedReverse(tf[[1]]$sigmaProportion)
  p <- levelplot(ncells~sigmaContinuous*sigmaProportion, tf[[1]], main='Number of cells filtered', col.regions=gray.colors)
print(p)
p <- levelplot(wss~sigmaContinuous*sigmaProportion, tf[[1]], main='Wighted Sum of Squares', col.regions=gray.colors)
print(p)

  p <- levelplot(ss~sigmaContinuous*sigmaProportion, tf[[1]], main='Sum of Squares', col.regions=gray.colors)
  print(p)

print(p)
p <- levelplot(R2~sigmaContinuous*sigmaProportion, tf[[1]], main='Average R2', col.regions=gray.colors)

print(p)
p <- xyplot(et.ref~et.comp|sigmaContinuous*sigmaProportion, concat[[1]], comparison=concat[[2]][, c('et.comp', 'et.ref')], panel=panel.shifts, pch='.', cex=2, strip=strip.custom(style=3, strip.names=c(T, T)))
print(p)
}

## Helper, rewrite character fields in dataframe as factors
## Lattice likes factors.
stringsAsFactors <- function(dataframe, otherfactors=NULL){
  for(i in seq_len(ncol(dataframe))){
    x <- dataframe[,i]
    if(is.character(x) || names(dataframe)[i] %in% otherfactors)
      dataframe[,i] <- as.factor(x)
  }
  dataframe
}

## Make a bunch of plots of filtering metrics
plotEtz <- function(mysc, groups, byGroup=FALSE, sigmaContinuous=c(3, 5, 7, 9)){
  etzmelt <- .getEtz(mysc, groups, byGroup, 5)
  p <- qqmath(~etz|gene, data=etzmelt, groups=eval(parse(text=groups)),
       panel=panel.qqall, pch="+", cex=1.3, ylim=c(-10, 10), layout=c(0, 48), main='Distribution of etz by group', alpha=.5)
  print(p)
  noutlier <- matrix(NA, nrow=nrow(mysc), ncol=length(sigmaContinuous), dimnames=list(getwellKey(mysc), sigmaContinuous))
  for(i in seq_along(sigmaContinuous))
    noutlier[,i] <- apply(abs(fl$z.exprs)>sigmaContinuous[i], 1, sum, na.rm=TRUE)

  moutlier <- melt(noutlier, varnames=c('wellkey', 'sigmaContinuous'))
  moutlier <- rename(moutlier, c('value'='noutlier'))
  p <- histogram(~noutlier|as.factor(sigmaContinuous), moutlier, subset=noutlier>0, breaks=c(0:10, 1000), type='count', main='Number of outliers at threshold', xlim=c(-1, 11))
  print(p)
  pet0zbind <- data.frame(cData(mysc), pet0z=fl$z.fet0)
  pet0zbind <- stringsAsFactors(pet0zbind)
  p <- qqmath(~pet0z|eval(groups), pet0zbind, tails.n=5,panel=function(...){
         panel.qqmath(...)
       panel.qqmathline(...)
       }, main='Distribution of pet0z by group')
  print(p)
}

.getEtz <- function(mysc, groups=NULL, byGroup, minThres){
 checkGroups(mysc, groups)
 conditionby <- NULL
  if(byGroup) conditionby <- groups
  latticeGroups <- paste(groups, collapse=':')
  fl <- filter(mysc, groups=conditionby, filt_control=list(filter=F), apply_filter=FALSE)
  if(byGroup && !is.null(conditionby)){
    z.exprs <- do.call(rbind, lapply(fl, function(x){x$z.exprs}))
  } else{
    z.exprs <- fl$z.exprs
  }
 
  allna <- apply(is.na(z.exprs), 2, all)
  z.exprs <- z.exprs[,!allna]
  anyOver <- apply(abs(z.exprs)>minThres, 1, any, na.rm=TRUE)
  flbind <- data.frame(cData(mysc), wellKey=getwellKey(mysc), anyOver)
  flbind <- merge(flbind, z.exprs, by=0)[,-1] #drop Row.names column
  etzmelt <- melt(flbind, id.vars=c('anyOver', names(cData(mysc)), 'wellKey'))
  etzmelt <- rename(etzmelt, c('value'='etz', 'variable'='gene'))
  etzmelt <- stringsAsFactors(etzmelt, c(groups, 'gene'))
  etzmelt
}

dotplotFilters <- function(mysc, groups, byGroup=FALSE, minThres=5){
   etzmelt <- .getEtz(mysc, groups, byGroup=byGroup, minThres)
   etzmelt$etztrunc <- with(etzmelt, sign(etz)*pmin(abs(etz), 10))
   grp <- with(etzmelt, as.factor(ifelse(anyOver, wellKey, 'clean')))
   nlev <- length(levels(grp))-1
               
   dotplot(~etztrunc|gene, groups=grp, etzmelt, jitter.y=TRUE, amount=.2, pch=c(rep(2:10, length.out=nlev), 1), cex=c(rep(.7, nlev), .3))
   
}



panel.qqall <- function(x, subscripts, groups, distribution=qnorm, ...){
  if(missing(subscripts) || is.null(subscripts)) subscripts <- TRUE
  ord <- order(x)
  sx <- sort(x)
  py <- ppoints(sx)
  grp <- groups[subscripts][ord]
  panel.xyplot(x=distribution(py), y=sx, groups=grp, subscripts=TRUE, ...)
  panel.qqmathline(x, y=x, distribution, groups=NULL, subscripts=TRUE, ...)
  panel.abline(h=c(-5, -3, 3, 5), col="grey")
}

##' @importFrom latticeExtra panel.smoother
pairsPlots <- function(mysc, FREQ_CUTOFF=.7, include=NULL){
#  require(latticeExtra)
genelist <- union(names(which(freq(mysc)>FREQ_CUTOFF)), include)
  theexp <- exprs(mysc)
  theexp[theexp==0] <- NA
  filt.crit <- apply(filter(mysc, apply_filter=F), 1, any)
   p <- splom(~theexp[, genelist], pch='+', groups=filt.crit, panel=function(x, y, groups, ...){
     panel.xyplot(x, y, groups=groups, ...)
     panel.smoother(x[!groups],y[!groups], method='lm', ...)
   } )
  print(p)
}

filtDiffs <- function(mysc, groups, filt_control=NULL){
  checkGroups(mysc, groups)
  if(is.null(filt_control)) filt_control <- list(filter=FALSE)
  fo <- filter(mysc, groups, filt_control, apply_filter=F)
  fglobal <- filter(mysc, filt_control=filt_control, apply_filter=F)
  fo <- c(fo, global=list(fglobal))
  fo.mad <- lapply(fo, '[[', f='et.mad')
  fo.med <- lapply(fo, '[[', f='et.med')
  fomelt <- melt(list(mad=fo.mad, med=fo.med))
  fomelt <- rename(fomelt, c('L2'=groups, 'L1'='metric'))
    fomelt$gene <- fData(mysc)[, getMapping(mysc,"primerid")]
  fomelt <- stringsAsFactors(fomelt)
  cb.lattice <- paste(groups, collapse=':')
 p <- barchart(~value|gene, groups=eval(parse(text=cb.lattice)), subset=metric=='med', data=fomelt, layout=c(0, 48), main='Filtering metrics by stratum (median)', auto.key=TRUE)
  print(p)
    p <- barchart(~value|gene, groups=eval(parse(text=cb.lattice)), subset=metric=='mad', data=fomelt, layout=c(0, 48), main='Filtering metrics by stratum (mad)', auto.key=TRUE)
  print(p)
}

burdenOfFiltering <- function(sc, groups, byGroup=FALSE, filt_control = NULL){
  checkGroups(sc, groups)
  conditionby <- NULL
  if(byGroup) conditionby <- groups
  filt <- filter(sc, groups=conditionby, filt_control=filt_control, apply_filter=FALSE)
  index <- apply(filt, 1, function(x){mx <- min(which(x)); if(!is.finite(mx)) mx <- 4; mx})
  outcome <- c(names(filt), 'none')[index]
  filt <- cbind(outcome, cData(sc))
  tab <- do.call(table, c(filt[, groups, drop=FALSE], filt[, 'outcome', drop=FALSE]))

  barchart(tab, main=sprintf('Burden of filtering by %s', paste(groups, collapse=" and ")), auto.key=TRUE)
}

#'Concordance plots of filtered single vs n-cell assays
#'
#'This will generate the types of plots shown in the bioinformatics manuscript
#'@param SCellAssay is a SingleCellAssay for the 1-cell per well assay
#'@param NCellAssay is a SingleCellAssay for the n-cell per well assay
#'@param filterCriteria is a list of filtering criteria to apply to the SCellAssay and NCellAssay
#'@param groups is a character vector naming the group within which to perform filtering. NULL by default.
#'@import ggplot2
#'@export
plotSCAConcordance<-function(SCellAssay, NCellAssay, filterCriteria=list(nOutlier = 2, sigmaContinuous = 9,sigmaProportion = 9), groups=NULL){
  if(!(inherits(SCellAssay,"SingleCellAssay")&inherits(NCellAssay,"SingleCellAssay"))){
    stop("SCellAssay and NCellAssay must inherit from SingleCellAssay");
  }
  filtered.sc<-filter(SCellAssay,filt_control=filterCriteria,groups=groups)
  filtered.hc<-filter(NCellAssay,filt_control=filterCriteria, groups=groups)
  concord.unfiltered<-getConcordance(SCellAssay,NCellAssay)
  concord.filtered<-getConcordance(filtered.sc,filtered.hc)
  toplot<-rbind(concord.filtered,concord.unfiltered)
  toplot<-data.frame(toplot,filter=c(rep("filtered",nrow(concord.filtered)),rep("unfiltered",nrow(concord.unfiltered))))
  
  g<-ggplot(toplot)+
    geom_point(data=subset(toplot,filter=="filtered"),aes(x=et.ref,y=et.comp),alpha=0.75)+
    theme_bw()+scale_x_continuous("Reference Assay")+scale_y_continuous("Comparison Assay")
 
  seg<-cast(melt(toplot,measure=c("et.ref","et.comp"),id=c("primerid","filter")),primerid~filter+variable)
  g<-g+geom_segment(data=seg,aes(x=unfiltered_et.ref,xend=filtered_et.ref,y=unfiltered_et.comp,yend=filtered_et.comp),alpha=0.5,lwd=0.5)
  
  wss.filt<-getwss(concord.filtered,concord.filtered$nexp.ref)
  wss.unfilt<-getwss(concord.unfiltered,concord.unfiltered$nexp.ref)
  message(sprintf("Sum of Squares before Filtering: %s\n After filtering: %s\n Difference: %s",round(wss.unfilt,2),round(wss.filt,2),round(wss.unfilt-wss.filt,2)))
  g
}
