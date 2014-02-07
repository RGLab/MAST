##' Estimate thresholds for positive expression
##'
##' @param nsa NanostringAssay object
##' @return modified nsa
setGeneric('thresholdNanoString', function(nsa, ...) standardGeneric('thresholdNanoString'))


##' Estimate thresholds for positive expression
##'
##' Estimates per-gene x unit thresholds for positive expression and truncates values below this threshold
##' Uncertain values (in terms of posterior probability of membership) can be set to NA or rounded left or right
##' Thresholds are estimated using a Gaussian mixture model with prior supplied by population estimates.
##' @param nsa NanostringAssay object
##' @param include.primers primeris to use for population estimates.  If missing, then all primers are used.
##' @param exclude.primers primers to exclude from population estimates
##' @param posteriorprob currently ignored
##' @param clip currently ignored
##' @param debug should extra fields helpful for examining the thresholding be returned as a list?
##' @param location.strength scaling of prior on location
##' @param pseudo.counts total strength of prior vs data
##' @param hard.threshold location estimates less than this value will be treated as belonging to the "noise cluster"
##' @param startLayer which layer should be used to threshold.  Default 'lCount'.
##' @aliases thresholdNanoString,NanoStringAssay-method
##' @return object of class \code{ThresholdedNanoString} if \code{debug=TRUE}, otherwise returns \code{nsa} with the thresholded expression in layer \code{et}
##' @seealso ThresholdedNanoString-class
setMethod('thresholdNanoString', signature='NanoStringAssay', function(nsa, include.primers, exclude.primers, posteriorprob, clip=c('left', 'right', 'NA'), debug=FALSE, location.strength=1, pseudo.counts=5, hard.threshold=3, startLayer='lCount'){
  layer(nsa) <- startLayer
  if(!('et' %in% dimnames(nsa)[[3]])) nsa <- addlayer(nsa, 'et')
  
  if(missing(include.primers)){
    include.primers <- fData(nsa)$primerid
  } else if(length(setdiff(include.primers, fData(nsa)$primerid))>0){
    warning('include.primers member ', setdiff(include.primers, fData(nsa)$primerid), ' not found.')
  }
  if(!missing(exclude.primers)) include.primers <- setdiff(include.primers, exclude.primers)
  nsa <- sort(nsa)
  m <- melt(nsa)
  dat <- subset(m, primerid %in% include.primers)[, startLayer]
  
  fc <- flowClust::flowClust(dat, K=2, trans=0)
  prior<-flowClust::flowClust2Prior(fc,kappa=3,Nt=5) ## Lambda and Omega both scale by kappa*Nt.  w0 scales by Nt alone (but we don't use w0)
  prior$Omega0 <- prior$Omega0*5/location.strength        #vague mean hyperprior
  prior$Lambda0 <- prior$Lambda0/2     #smaller variances
  prior$w0 <- c(pseudo.counts, pseudo.counts)                     #cluster membership prior
  if(debug) cat('Omega0', prior$Omega0, ' Lambda0 ', prior$Lambda0)
  dats <- split(m[, startLayer], m$primerid)
  fitgene <- mclapply(dats, function(set){
    fc.tmp <- flowClust::flowClust(set, K=2, prior=prior, u.cutoff=NULL, level=1, usePrior='yes', trans=0, z.cutoff=0)
    fc.tmp
  })
  means <- do.call(rbind, lapply(fitgene, function(x){t(x@mu)}))
  props <- do.call(rbind, lapply(fitgene, function(x){t(x@w)}))
  w.max <- apply(means, 1, which.max)
  p.signal <- do.call(c, lapply(1:length(fitgene), function(x){
    fitgene[[x]]@z[,w.max[x]]
  })) #genes in alpha-order, then each idvars in alpha order, same as melted nsa

  lab <- do.call(c, lapply(1:length(fitgene), function(x){
      label <- flowClust::Map(fitgene[[x]])
    if(means[x,w.max[x]] < hard.threshold){
      rep(2, length(label)) #truncate all if both cluster means are less than hard.treshold
    } else if(w.max[x]==1){
      label
    } else if(w.max[x]==2){
      c(2, 1)[label]
    } else{
      stop('ruhroh, there should only be two clusters')
    }
  })) #genes in alpha-order, then each idvars in alpha order, same as melted nsa

  densities <- lapply(fitgene, function(x){
      mu <- x@mu
      sigma <- x@sigma
      w <- x@w
      f <- function(y){
          
          dnorm(y, mu[1], sigma[1])*w[1] + dnorm(y, mu[2], sigma[2])*w[2]
      }
      return(f)
  })

  m <- cbind(m, ps=p.signal, clusterID=as.factor(ifelse(is.na(lab), 3, lab)))
  layer(nsa) <- 'et'
  exprs(nsa) <- ifelse(m$clusterID==1, m[,startLayer], 0) #fix rounding
  if(debug) return(new('ThresholdedNanoString', melted=m, nsa=nsa, densities=densities, means=means, props=props, startLayer=startLayer))
  else return(nsa)
})

setMethod("show","ThresholdedNanoString",function(object){
  cat(class(object), ' of ', nrow(object@nsa), ' rows and ', ncol(object@nsa), ' columns')
  invisible(NULL)
})


##' Plot NanoString log Counts vs cluster ID
##'
##' Histogram of log Counts by primerid, colored by cluster
##' @param thresholdedNanoString output from thresholdNanoString, debug=TRUE
##' @param primerids character vector of primerids to plot
##' @return ggplot object
plot.threshold <- function(thresholdedNanoString, primerids){
    if(!inherits(thresholdedNanoString, 'ThresholdedNanoString')) stop('thresholdedNanoString must inherit from class ThresholdedNanoString')
    sub <- subset(thresholdedNanoString@melted, primerid %in% primerids)
    measure <-  thresholdedNanoString@startLayer                        #
sub <- ddply(sub, .(primerid), function(x){
    x$den.est <- thresholdedNanoString@densities[[x$primerid[1]]](x[, measure])
    x
})
p <- ggplot(sub, aes_string(x=measure, fill='clusterID')) + geom_histogram(aes(y=..density..)) + facet_wrap(~primerid) + geom_line(aes_string(x=measure, y='den.est'))+ylim(0, 1)
    p
}
