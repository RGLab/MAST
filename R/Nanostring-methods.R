#'@exportMethod thresholdNanoString
setGeneric('thresholdNanoString', function(nsa, ...) standardGeneric('thresholdNanoString'))

#' Estimate thresholds for positive expression
#'
#' Estimates per-gene x unit thresholds for positive expression and truncates values below this threshold
#' Uncertain values (in terms of posterior probability of membership) can be set to NA or rounded left or right
#' Thresholds are estimated using a Gaussian mixture model with prior supplied by population estimates.
#' @name threshold
#' @param nsa NanostringAssay object
#' @param groups groups to apply thresholding
#' @param thresholds data.frame of thresholds  of groups x genes user may specify
#' @param posteriorprob Min posterior probability of cluster membership for an observation to be truncated 
#' @param clip Should values with uncertain posterior probs be clipped
#' @return modified nsa or list with elements nsa and debugging info regarding the clustering
setMethod('thresholdNanoString', signature='NanoStringAssay', function(nsa, include.primers, exclude.primers, posteriorprob, clip=c('left', 'right', 'NA'), debug=FALSE){
  layer(nsa) <- 'lCount'
  if(!('et' %in% dimnames(nsa)[[3]])) nsa <- addlayer(nsa, 'et')

  if(missing(include.primers)){
    include.primers <- fData(nsa)$primerid
  } else if(length(setdiff(include.primers, fData(nsa)$primerid))>0){
    warning('include.primers member ', setdiff(include.primers, fData(nsa)$primerid), ' not found.')
  }
  if(!missing(exclude.primers)) include.primers <- setdiff(include.primers, exclude.primers)
  m <- melt(nsa)
  dat <- subset(m, primerid %in% include.primers)$lCount
  
  fc <- flowClust::flowClust(dat, K=2, trans=0)
  prior<-flowClust::flowClust2Prior(fc,kappa=1,Nt=50) ## Lambda and Omega both scale by kappa*Nt.  w0 scales by Nt alone (but we don't use w0)
  prior$Omega0 <- prior$Omega0*5        #vague mean hyperprior
  prior$Lambda0 <- prior$Lambda0/2     #smaller variances
  prior$w0 <- c(5, 5)                     #vague cluster membership prior, but weighted towards noise cluster
  dats <- split(m$lCount, m$primerid)
  fitgene <- mclapply(dats, function(set){
    fc.tmp <- flowClust::flowClust(set, K=2, prior=prior, level=.9, usePrior='yes', trans=0, z.cutoff=.95)
    fc.tmp
  })
  means <- do.call(rbind, lapply(fitgene, function(x){t(x@mu)}))
  w.max <- apply(means, 1, which.max)
  p.signal <- do.call(c, lapply(1:length(fitgene), function(x){
    fitgene[[x]]@z[,w.max[x]]
  })) #genes in alpha-order, then each idvars in alpha order, same as melted nsa

  lab <- do.call(c, lapply(1:length(fitgene), function(x){
    if(w.max[x]==1){
      fitgene[[x]]@label
    } else if(w.max[x]==2){
      c(2, 1)[fitgene[[x]]@label]
    } else{
      stop('ruhroh, there should only be two clusters')
    }
  })) #genes in alpha-order, then each idvars in alpha order, same as melted nsa


  m <- cbind(m, ps=p.signal, clusterID=as.factor(ifelse(is.na(lab), 3, lab)))
  layer(nsa) <- 'et'
  exprs(nsa) <- ifelse(m$clusterID==1, m$lCount, 0) #fix rounding
  if(debug) return(list(m=m, nsa=nsa))
  else return(nsa)
})
