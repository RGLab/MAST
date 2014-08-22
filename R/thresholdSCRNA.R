
find_peaks<-function (x, y = NULL, num_peaks = NULL, adjust = 2, plot = FALSE, 
                      ...) 
{
  x <- as.vector(x)
  if (length(x) < 2) {
    warning("At least 2 observations must be given in 'x' to find peaks.")
    return(NA)
  }
  if (is.null(y)) {
    dens <- density(x, adjust = adjust, ...)
  }
  else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of 'x' and 'y' must be equal.")
    }
    dens <- list(x = x, y = y)
  }
  second_deriv <- diff(sign(diff(dens$y)))
  which_maxima <- which(second_deriv == -2) + 1
  which_maxima <- which_maxima[findInterval(dens$x[which_maxima], 
                                            range(x)) == 1]
  which_maxima <- which_maxima[order(dens$y[which_maxima], 
                                     decreasing = TRUE)]
  if (length(which_maxima) > 0) {
    peaks <- dens$x[which_maxima]
    if (is.null(num_peaks) || num_peaks > length(peaks)) {
      num_peaks <- length(peaks)
    }
    peaks <- peaks[seq_len(num_peaks)]
  }
  else {
    peaks <- NA
  }
  peaks <- data.frame(x = peaks, y = dens$y[which_maxima][seq_len(num_peaks)])
  if (plot) {
    plot(dens, main = paste("adjust =", adjust))
    points(peaks, , col = "red")
  }
  peaks
}
find_valleys<-function (x, y = NULL, num_valleys = NULL, adjust = 2, ...) 
{
  x <- as.vector(x)
  if (length(x) < 2) {
    warning("At least 2 observations must be given in 'x' to find valleys.")
    return(NA)
  }
  if (is.null(y)) {
    dens <- density(x, adjust = adjust, ...)
  }
  else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of 'x' and 'y' must be equal.")
    }
    dens <- list(x = x, y = y)
  }
  second_deriv <- diff(sign(diff(dens$y)))
  which_minima <- which(second_deriv == 2) + 1
  which_minima <- which_minima[findInterval(dens$x[which_minima], 
                                            range(x)) == 1]
  which_minima <- which_minima[order(dens$y[which_minima], 
                                     decreasing = FALSE)]
  if (length(which_minima) > 0) {
    valleys <- dens$x[which_minima]
    if (is.null(num_valleys) || num_valleys > length(valleys)) {
      num_valleys <- length(valleys)
    }
    valleys <- valleys[seq_len(num_valleys)]
  }
  else {
    valleys <- NA
  }
  valleys
}
between_interval<-function (x, interval) 
{
  x <- x[findInterval(x, interval) == 1]
  if (length(x) == 0) {
    x <- NA
  }
  x
}

#'Threshold a count matrix using an adaptive threshold.
#'
#'An adaptive threshold is calculated from the conditional mean of expression, based on 10 bins
#'of the genes with similar expression levels. Thresholds are chosen by estimating cutpoints in the bimodal density estimates of the
#'binned data. 
#'@param data \code{matrix} of counts
#'@param bin_by \code{character} "mean"
#'@param plot \code{logical}.
#'@return \code{list} of thresholded counts, thresholds, and bins
#'@export
thresholdSCRNACountMatrix<-function(data,bins=10,main="",bin_by="mean",plot=FALSE){
  arg<-match.arg(bin_by,c("proportion","mean"))
  if(bins!=10){
    message("bins is ignored, using default of 10")
    bins<-10
  }
  if(arg=="mean"){
    proportions<-apply(data,1,function(x)log1p(mean(x[x>0],na.rm=TRUE)))
    proportions[is.nan(proportions)]<-0
  }else{
    proportions<-rowMeans(data>0)
  }
  proportion_bins<-cut(proportions, bins)
  peaks<-by(data,proportion_bins,function(x)find_peaks(unlist(log1p(x)),adjust = 2))
  valleys<-by(data,proportion_bins,function(x)find_valleys(unlist(log1p(x)),adjust = 2))
  dens<-by(data,proportion_bins,function(x)density(unlist(log1p(x)),adjust=2))
  single_modes<-do.call(c,lapply(peaks,function(x)abs(diff(x[1:2,1]))))<1|(lapply(list(do.call(c,lapply(peaks,function(x)abs(diff(x[1:2,2]))))),function(x)(x-median(na.omit(x)))/mad(na.omit(x)))[[1]]>2)
  cutpoints<-lapply(seq_along(peaks),function(j){
    valleys[[j]][which(findInterval(valleys[[j]],sort(peaks[[j]][1:2,1]))==1)]
  })
  cutpoints[single_modes]<-NA
  cutpoints<-lapply(cutpoints,function(x)ifelse(length(x)==0,NA,max(x)))
  names(cutpoints)<-names(peaks)
  
  #impute cutpoints if NA
  for(i in 1:length(cutpoints)){
    if(is.na(cutpoints[[i]])){
      if(i!=length(cutpoints)){
        cutpoints[[i]]<-do.call(c,cutpoints)[min(which(!is.na(do.call(c,cutpoints))))]
      }else{
        cutpoints[[i]]<-do.call(c,cutpoints)[max(which(!is.na(do.call(c,cutpoints))))]
      }
    }
  }
  
  #ensure cutpoints are increasing
  for(i in length(cutpoints):2){
    if(cutpoints[[i-1]]>cutpoints[[i]]){
      cutpoints[[i-1]]<-cutpoints[[i]]
    }
  }
  
  data_threshold<-log1p(data)
  for(i in levels(proportion_bins)){
    data_threshold[proportion_bins%in%i,][data_threshold[proportion_bins%in%i,]<cutpoints[i]]<-0
  }
  if(plot){
    par(mfrow=c(2,5))
    for(i in 1:bins){
      xl<-c(0,max(dens[[i]]$x,unlist(cutpoints)))
      plot(dens[[i]],xlim=xl,xlab="Conditional Mean Expression",ylab="Density",main=paste0("Bin=",levels(proportion_bins)[[i]]));
      abline(v=cutpoints[[i]]) 
    }
  }
  nms<-names(cutpoints)
  cutpoints<-unlist(cutpoints)
  names(cutpoints)<-nms
  list(counts_threshold=expm1(data_threshold),cutpoint=cutpoints,bin=proportion_bins)
}

