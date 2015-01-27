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
    return( peaks )
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
    return( valleys )
}

between_interval<-function (x, interval)
{
    x <- x[findInterval(x, interval) == 1]
    if (length(x) == 0) {
        x <- NA
    }
    return( x )
}

apply_by<-function(x,by_idx,fun,...){
    res<-sapply(unique(by_idx),function(a){
        apply(x[,by_idx==a],1,fun,...)}
    )
    return(res)
}


#'Threshold a count matrix using an adaptive threshold.
#'
#'An adaptive threshold is calculated from the conditional mean of expression, based on 10 bins
#'of the genes with similar expression levels. Thresholds are chosen by estimating cutpoints in the bimodal density estimates of the
#'binned data.
#' @param data_all \code{matrix} of counts.  Rows are cells and columns are genes.
#' @param conditions Bins are be determined per gene and per condition.  Typically contrasts of interest should be specified.
#' @param cutbins \code{vector} of cut points.
#' @param nbins \code{integer} number of bins when cutbins is not specified.
#' @param bin_by \code{character} "median", "proportion", "mean"
#' @param qt when \code{bin_by} is "quantile", what quantile should be used to form the bins
#' @param min_per_bin minimum number of genes within a bin
#' @param absolute_min \code{numeric} giving a hard threshold below which everything is assumed to be noise
#' @param return_log return the logged expression matrix or not.  By default, returned expression matrix will be logged ( base 2 ).
#' @param G \code{integer} number of mixture model components to fit to the data above the threshold (default 2)
#' @param plot \code{logical} whether or not to plot thresholding results
#'@return \code{list} of thresholded counts (on natural scale), thresholds, bins, densities estimated on each bin, as well as null component weights for each observation from the mixture model fit, and the original data
#'@importFrom plyr ldply
#'@export
thresholdSCRNACountMatrix <-function( data_all              ,
                                      conditions  = NULL    ,
                                      cutbins     = NULL    ,
                                      nbins       = 10      ,
                                      bin_by      = "median",
                                      qt          = 0.975,
                                      min_per_bin = 50      ,
                                      absolute_min= 0.0     ,
                                      return_log  = TRUE,
                                      G=2,
                                      plot=FALSE
                                    )
{

    data_all <- t(data_all)   #internally we work with the data having rows as genes and columns as cells
                                        # when there is no condition to stratefy
    log_base <- 2
    if( is.null( conditions ) ){ 
        conditions <- rep( 1, dim( data_all )[2] ) 
    } else { 
        conditions <- as.character( conditions )
    }
    comp_zero_idx <- rowSums( log( data_all+1, base = log_base )> 0.0 ) == 0
    data          <- data_all[!comp_zero_idx,]
    uni_cond      <- unique( conditions )
    log_data      <- log( data + 1, base = log_base )

    arg <- match.arg( bin_by, c( "quantile", "proportion", "mean", "median","iqr" ) )
    if( arg == "median" ){
        cond_stat <- apply_by( log_data, conditions, function( x ) if( all( x <= 0.0 ) ){ 0 } else { median( x[x>0.0]) })
    } else if( arg == "mean" ){
        # mean is taken on the original scale than loged
        cond_stat <- apply_by( data, conditions, function( x ) if( all(x <= 0.0 ) ){ 0 } else { log( mean( x[x>0.0])+1, base = log_base) } )
    } else if( arg == "proportion" ){
        cond_stat <- apply_by( data > 0, conditions, mean )
    } else if( arg == "quantile" ){
        cond_stat <- apply_by( log_data, conditions, function(x) if( all( x <= 0.0 ) ){ 0 } else { quantile( x[x>0.0], qt ) } )
    } else if( arg == "iqr" ){
        cond_stat <- apply_by( log_data, conditions, function(x) if( all( x <= 0.0 ) ){ 0 } else { IQR( x[x>0.0]) } )
    } else {
        stop("choose bin_by from c('proportion', 'mean', 'median' )")
    }
    rprop     <- range( unlist(cond_stat) )
    # log bin
    if( is.null(cutbins) ){
        rprop     <- range( cond_stat )
        # log bin
        cutbins   <- exp( seq( log( rprop[1] + 1 ),log( rprop[2] + 1 ),length.out = nbins + 1 ) ) - 1
        cutbins[1]                  <- min( cutbins[1], rprop[1] - 1e-10 )
        cutbins[ length( cutbins ) ]<- max( cutbins[length( cutbins )], rprop[2] + 1e-10 )
    } else{
        if( !( min(cutbins)< min(cond_stat) & max(cutbins) > max(cond_stat) ) ){
            stop("range of cutbins must include the range of the log expression")
        }
    }
    cond_stat_bins <- cut( cond_stat, cutbins )
    t_prop_bins    <- table( cond_stat_bins )
    cond_stat_bins_array             <- array( cond_stat_bins, dim( cond_stat ) )
    dimnames( cond_stat_bins_array ) <- dimnames( cond_stat )

    # make sure each bin has at least min_per_bin genes
    while( any( t_prop_bins < min_per_bin ) ){
        i <- which.min( t_prop_bins )
        if( i == 1 & t_prop_bins[ i ] < min_per_bin ){
            cutbins <- cutbins[-( i + 1 )]
        } else if( i == length( t_prop_bins ) & t_prop_bins[i] < min_per_bin ){
            cutbins <- cutbins[ -i ]
        } else if( t_prop_bins[ i ] < min_per_bin ){
            cutbins <- cutbins[ -ifelse( t_prop_bins[i-1] < t_prop_bins[ i + 1 ], i, ( i + 1 ) ) ]
        }
        cond_stat_bins <- cut( cond_stat, cutbins )
        t_prop_bins    <- table( cond_stat_bins )
    }
    cond_stat_bins_array<-array(cond_stat_bins, dim(cond_stat))
    dimnames(cond_stat_bins_array)<-dimnames(cond_stat)

    # convert data into a list
    log_data_list <- vector("list", length(levels(cond_stat_bins)))
    names(log_data_list) <- levels(cond_stat_bins)
    for( i in levels(cond_stat_bins)){
        log_data_list[[i]]<-NULL
        for( j in uni_cond){
            x<-unlist(log_data[cond_stat_bins_array[,j]==i,conditions==j])
            log_data_list[[i]]<-c(log_data_list[[i]],x[x>0.0])
        }
    }
    #dens   <- lapply( log_data_list, function( x )      density(         x, adjust = 1 ) )
    #peaks  <- lapply(          dens, function( dd )  find_peaks( dd$x,dd$y, adjust = 1 ) )
    #valleys<- lapply(          dens, function( dd )find_valleys( dd$x,dd$y, adjust = 1 ) )
    dens   <- lapply( log_data_list, function( x )  { if(length(x)>2){ density( x, adjust = 1 ) } else{ NULL } })
    peaks  <- lapply(          dens, function( dd ) { if( is.null(dd) ){ data.frame( x=0, y=0 ) }else{ find_peaks( dd$x,dd$y, adjust = 1 )}} ) 
    valleys<- lapply(          dens, function( dd ) { if( is.null(dd) ){ list(0)} else{find_valleys( dd$x,dd$y, adjust = 1 ) } })

    single_modes<-do.call(c,lapply(peaks,function(x)abs(diff(x[1:2,1]))))<1|(lapply(list(do.call(c,lapply(peaks,function(x)abs(diff(x[1:2,2]))))),function(x)(x-median(na.omit(x)))/mad(na.omit(x)))[[1]]>2)
    #Check for single peaks found
    singles <- ldply( lapply( peaks,nrow ),function( x ) x == 1 )[,2]
    single_modes[singles] <- TRUE

    cutpoints<-lapply( seq_along(peaks), function( j ){
        valleys[[j]][which( findInterval( valleys[[j]], sort( peaks[[j]][1:2,1]) ) == 1 )]
    })
    cutpoints[single_modes]<-NA
    cutpoints       <- lapply( cutpoints, function(x) ifelse( length( x )==0, NA, max( x ) ) )
    names(cutpoints)<- names( peaks )

    #Check if all cutpoints are NA
    if( all( is.na( ldply( cutpoints )[,2] ) ) ){
        for( i in 1:length( cutpoints ) ){
            cutpoints[[i]] <- 0.0 #0
        }
    } else {
        #impute cutpoints if NA
        for(i in 1:length( cutpoints ) ){
            if( is.na( cutpoints[[i]] ) ){
                if( i != length( cutpoints ) ){
                    cutpoints[[i]]<-do.call(c,cutpoints)[min(which(!is.na(do.call(c,cutpoints))))]
                }else{
                    cutpoints[[i]]<-do.call(c,cutpoints)[max(which(!is.na(do.call(c,cutpoints))))]
                }
            }
        }
    }
    # ensure cutpoints are increasing and decreasing from the 75% of bins with 2 modes.
    # could be improved. 
    vals2    <- unlist(valleys[which(unlist(lapply(peaks,nrow))==2)])
    midindex <- which(names(peaks)==names(which.min(abs(vals2-quantile(vals2,c(0.75),na.rm=TRUE)))))
    if( length(midindex) > 0 ){
        for( i in midindex:2 ){
            if( cutpoints[[i-1]] > cutpoints[[i]] ){
                cutpoints[[i-1]] <- cutpoints[[i]]
            }
        }
        for( i in midindex:(length( cutpoints )-1) ){
            if( cutpoints[[i+1]] <= cutpoints[[i]] ){
                cutpoints[[i+1]] <- cutpoints[[i]]
            }
        }
    } else { # when no clear 2 peaked distribution exists, start from the top
        for( i in length( cutpoints ):2 ){
            if( cutpoints[[i-1]] > cutpoints[[i]] ){
                cutpoints[[i-1]] <- cutpoints[[i]]
            }
        }

    }
    cutpoints <- lapply( cutpoints, function(x) max( absolute_min ,x ) )

    if(return_log){
         data_threshold <- log_data
    } else {
         data_threshold <- data#log_data
    }
    
   #' Mixture model fit given cutpoints and bins
   Zout<-thresholdMMfit(log_data_list=log_data_list,cutpoints=cutpoints,plot=plot,G=G)
   data_null_weights<-matrix(1,nrow(data),ncol(data))
   dimnames(data_null_weights)<-dimnames(data)
   for( i in levels(cond_stat_bins)){
     log_data_list[[i]]<-NULL
     for( j in uni_cond){
       #grab the data originally in the bin
       x<-unlist(log_data[cond_stat_bins_array[,j]==i,conditions==j])
       #grab the weights in the same fashion and set them to the output for the non-zero elements.
       data_null_weights[cond_stat_bins_array[,j]==i,conditions==j][x>0.0]<-Zout[[i]][,1]
     }
   }
   
    for( j in uni_cond ){
        for( i in levels(cond_stat_bins) ){
            if(any(cond_stat_bins_array[,j]==i))
            data_threshold[cond_stat_bins_array[,j]==i,conditions==j][log_data[cond_stat_bins_array[,j]==i,conditions==j]<cutpoints[i]]<-0
        }
    }
   
    nms                <- names(  cutpoints )
    cutpoints          <- unlist( cutpoints )
    names( cutpoints ) <- nms
    print( cutpoints )
    data_threshold_all                  <- data_all*0
    #data_threshold_all[comp_zero_idx, ] <- 0
    data_threshold_all[!comp_zero_idx,] <- data_threshold #2^( data_threshold ) - 1
    bin_all       <- factor(cond_stat_bins_array)
    dim(bin_all)  <- dim(cond_stat_bins_array)
    res_obj       <- list( counts_threshold = t(data_threshold_all),
                          original_data     = t(log_data),
                          cutpoint          = cutpoints,
                          bin               =  bin_all,
                          conditions        = conditions,
                          density           = dens, 
                          peaks             = peaks, 
                          valleys= valleys,
                          null_weights=data_null_weights)
    class( res_obj )<- c( "list", "thresholdSCRNACountMatrix" )
    return( res_obj )

}

##' Plot cutpoints and densities for thresholding
##'
##' @param object output of \code{thresholdSCRNACountMatrix}
##' @param ask if TRUE then will prompt before displaying each plot
##' @param wait.time pause (in seconds) between each plot
##' @param type one or more of the following: 'bin' (plot the genes by the binning used for thresholding), or 'gene' (plot thresholding by gene -- see next argument)
##' @param indices if \code{type} is equal to 'gene', and is a integer of length 1, then a random sample of \code{indices} genes is taken.  If it is NULL, then 10 genes are sampled.  If it is a integer vector of length > 1, then it is interpreted as giving a list of indices of genes to be displayed.
##' @param ... further arguments passed to \code{plot}
##' @return displays plots
##' @export
plot.thresholdSCRNACountMatrix<-function(object, ask=FALSE, wait.time=0, type='bin', indices=NULL, ...)
{
    type <- match.arg(type, c('bin', 'gene'), several.ok=TRUE)
    op <- par(ask=ask)
    par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)
    if('bin' %in% type){
        ## plot by bins
        for(i in 1:length(object$density)){
            if(!is.null(object$density[[i]])){
                plot(object$density[[i]],main=names(object$cutpoint)[i], ...)
                abline(v=object$cutpoint[i],col="red",lty=2)
                Sys.sleep(wait.time)
            }
        }
    }

    if('gene' %in% type){
        ## plot by genes
        uni_cond <- unique(object$conditions)
        num_cond <- as.numeric(as.factor(object$conditions))
        if(is.null(indices)){
            indices <- 10L
        }
        if(is.finite(indices) && length(indices)==1){
            indices <- sample(ncol(object$original_data), indices)
        }
        if(any(indices < 1 | indices>ncol(object$original_data))) stop('`indices` must range between 1 and the number of columns of the count matrix.')
        for(i in indices){
            den <- density(object$original_data[,i])
            plot(den, main=paste0('Index ', i, '(', colnames(object$original_data)[i], ')'), ...)
            for(j in seq_along(uni_cond)){
                abline(v=with(object, cutpoint[bin[i,j]]), col=j, lty=2)
                with(object, rug(original_data[j==num_cond,i], col=j))
            }
        }
    }
    par(op)
}

##' Summarize the effect of thresholding
##'
##' Returns the proportion of (putative) expression, the variance of expressed cells, and -log10 shapiro-wilk tests for normality on the expressed cells
##' @param object a \code{thresholdSCRNACountMatrix}
##' @param ... currently ignored
##' @return a list of statistics on the original data, and thresholded data
##' @export
summary.thresholdSCRNACountMatrix <- function(object, ...){
    ## original <- with(object, data.table:::melt.data.table(data.table(conditions=conditions, type='original_data', original_data), id.vars=c('conditions', 'type')))
    ## threshold <- with(object, data.table:::melt.data.table(data.table(conditions=conditions, type='counts_threshold', counts_threshold), id.vars=c('conditions', 'type')))
    ## both <- rbind(original, threshold)
    ## summaries <- both[,list(zeroes=mean(value>0),
    ##            vars=var(value[value>0]),
    ##                         shapiro={
    ##                             pos <- value[value>0]
    ##                             if(length(pos)>3 & length(pos) < 5000)            -log10(shapiro.test(pos)$p.value) else NA_real_
    ##                         })
    ##                  ,keyby=list(conditions, type, variable)]
    
                 
    zeros <- lapply(object[c('original_data', 'counts_threshold')], function(o){
        apply(o>0, 2, mean)
    })
    vars <- lapply(object[c('original_data', 'counts_threshold')], function(o){
        o[o==0] <- NA
        apply(o, 2, var, na.rm=TRUE)
    })
    shapiro <- lapply(object[c('original_data', 'counts_threshold')], function(o){
        apply(o, 2, function(x){
           
        })
    })
    out <- list(zeros=zeros, vars=vars, shapiro=shapiro)
    class(out) <- c(class(out), 'summaryThresholdSCRNA')
    out
}

##' @export
##' @describeIn summary.thresholdSCRNACountMatrix prints five-number distillation of the statistics and invisibly returns the table used to generate the summary
print.summaryThresholdSCRNA <- function(object, ...){
    class(object) <- class(object)[-length(class(object))]
    m <- as.data.table(reshape::melt.list(object, na.rm=TRUE))
    m <- m[!is.na(value),]
    setnames(m, 'L1', 'metric')
    summ <- m[,list(stat=names(summary(value)), value=summary(value)),keyby=list(L2, metric)]
    summ$stat <- factor(summ$stat, levels=c('Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.'))
    dc <- dcast.data.table(summ, stat + L2  ~  metric)
    
    setnames(dc, c('stat', 'L2', 'vars', 'zeros'), c('Quantile', 'Data', 'variance', 'prop. expressed'))
    dcdf <- as.data.frame(dc)
    print(dcdf, right=FALSE, row.names=FALSE)
    invisible(dc)
}

#'thresholdMMfit
#'
#'Fit a gamma distribution to the data below the cutpoint. Fit a G component mixture to the data above the cutpoint.
#'Compute the full density and calculate the mixture weights for each observation.
#'@param log_data_list the list of binned data
#'@param cutpoints the vector of cutpoints
#'@param plot logical whether to show plots or not
#'@param G number of components
#'@author Greg Finak
#'@import mclust
#'@export
thresholdMMfit<-function(log_data_list=NULL,cutpoints=NULL,plot=FALSE,G=3){
  options(warn=-1)
  if(plot==TRUE){
    if(length(cutpoints)>4){
      x<-length(cutpoints)%%4
    }else{
      
    }
    if(x>4){
      x<-4
    }
    #4 x 4 grid at most
    par(mfrow=c(x,4))
  }
  Zout<-list()
  message("Fitting mixtures to binned genes")
  for(i in names(cutpoints)){
    dat<-log_data_list[[i]]
    cut<-cutpoints[[i]]
    
    #'3-component Gaussian fit to positive part
    fit_pos<-densityMclust(dat[which(dat>cut)],G=G)
    
    #'Gamma fit to negative part
    #'Gradient Descent fit, starting at method of moment estimates
    neg_data<-dat[dat<=cut]
    ll<-function(par,data){-sum(dgamma(data,log=TRUE,shape = par[1],scale=par[2]))}
    opt_fit<-optim(par=c(mean(neg_data)/(mean(neg_data)-mean(log(neg_data))),mean(neg_data)-mean(log(neg_data))),fn=ll,data=neg_data,method = "L-BFGS-B",lower = c(1e-10,1e-10))
    #'Check for convergence
    if(opt_fit$convergence!=0){
      #Try Nelder Mead
      opt_fit<-optim(par=c(mean(neg_data)/(mean(neg_data)-mean(log(neg_data))),mean(neg_data)-mean(log(neg_data))),fn=ll,data=neg_data,method="Nelder-Mead")
      if(opt_fit$convergence!=0){
        warning("Fitting gamma failed to converge during thresholding for bin ",i)
      }
    }
    
    #'Mixture proportions
    p<-length(neg_data)/(length(neg_data)+fit_pos$n)
    
    #'Density of negative and positive parts
    Dneg<-dgamma(dat,shape=opt_fit$par[1],scale=opt_fit$par[2])
    Dpos<-dens(fit_pos$modelName,data = dat,logarithm = FALSE,parameters = fit_pos$parameters)
    
    #'Weights
    z<-(p*Dneg)/((p*Dneg+(1-p)*Dpos))
    Z<-cbind(z,1-z)
    Z[order(dat),1]<-smooth(Z[order(dat),1])
    Z[,2]<-1-Z[,1]
    if(plot==TRUE){
      plot(sort(dat),(p*Dneg+(1-p)*Dpos)[order(dat)],type="l",xlab="log(tpm)",ylab="density",main=paste0("Bin ",i))
      hist(dat,prob=TRUE,add=TRUE,50)
      abline(v=cut,lwd=2,col="black")
      lines(x=sort(dat),Z[order(dat),1],type="l",col="red",lwd=1)
      lines(x=sort(dat),Z[order(dat),2],type="l",col="green",lwd=1)
    }
    Zout[[i]]<-Z
  }
  message("Done")
  names(Zout)<-names(cutpoints)
  options(warn=0)
  return(Zout)
}

