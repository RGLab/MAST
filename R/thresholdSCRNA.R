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
#' @param data_all \code{matrix} of counts
#' @param conditions ??
#' @param cutbins \code{vector} of cut points.
#' @param nbins \code{integer} number of bins when cutbins is not specified.
#' @param bin_by \code{character} "median", "proportion", "mean"
#' @param qt ??
#' @param min_per_bin ??
#' @param absolute_min ??
#' @param log_base ??
#' @param plot \code{logical}.
#'@return \code{list} of thresholded counts, thresholds, and bins
#'@importFrom plyr ldply
#'@export
thresholdSCRNACountMatrix <-function( data_all              ,
                                      conditions  = NULL    ,
                                      cutbins     = NULL    ,
                                      nbins       = 10      ,
                                      bin_by      = "median",
                                      qt          = 0.975,
                                      min_per_bin = 50      ,
                                      absolute_min= 1.0     ,
                                      log_base    = 2
                                    )
{

    # when there is no condition to stratefy
    if( is.null( conditions ) ) conditions <- rep( 1, dim( data_all )[2] )
    comp_zero_idx <- rowSums( log( data_all+1, base = log_base )> absolute_min ) == 0
    data          <- data_all[!comp_zero_idx,]
    uni_cond      <- unique( conditions )
    log_data      <- log( data + 1, base = log_base )

    arg <- match.arg( bin_by, c( "quantile", "proportion", "mean", "median","iqr" ) )
    if( arg == "median" ){
        cond_stat <- apply_by( log_data, conditions, function( x ) if( all( x <= absolute_min ) ){ 0 } else { median( x[x>absolute_min]) })
    } else if( arg == "mean" ){
        # mean is taken on the original scale than loged
        cond_stat <- apply_by( data, conditions, function( x ) if( all(x <= absolute_min ) ){ 0 } else { log( mean( x[x>absolute_min])+1, base = log_base) } )
    } else if( arg == "proportion" ){
        cond_stat <- apply_by( data > 0, conditions, mean )
    } else if( arg == "quantile" ){
        cond_stat <- apply_by( log_data, conditions, function(x) if( all( x <= absolute_min ) ){ 0 } else { quantile( x[x>absolute_min], qt ) } )
    } else if( arg == "iqr" ){
        cond_stat <- apply_by( log_data, conditions, function(x) if( all( x <= absolute_min ) ){ 0 } else { IQR( x[x>absolute_min]) } )
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
            log_data_list[[i]]<-c(log_data_list[[i]],x[x>absolute_min])
        }

    }
    dens   <- lapply( log_data_list, function( x )      density(         x, adjust = 1 ) )
    peaks  <- lapply(          dens, function( dd )  find_peaks( dd$x,dd$y, adjust = 1 ) )
    valleys<- lapply(          dens, function( dd )find_valleys( dd$x,dd$y, adjust = 1 ) )

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
            cutpoints[[i]] <- absolute_min #0
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
        #ensure cutpoints are increasing
        for( i in length( cutpoints ):2 ){
            if( cutpoints[[i-1]] > cutpoints[[i]] ){
                cutpoints[[i-1]] <- cutpoints[[i]]
            }
        }
        cutpoints <- lapply( cutpoints, function(x) max( absolute_min ,x ) )
    }
    data_threshold <- data#log_data
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
    data_threshold_all                  <- data_all
    data_threshold_all[comp_zero_idx, ] <- 0
    data_threshold_all[!comp_zero_idx,] <- data_threshold #2^( data_threshold ) - 1
    bin_all        <- cond_stat_bins_array
    res_obj       <- list( counts_threshold = data_threshold_all ,
                           cutpoint         = cutpoints          ,
                           bin              = factor( bin_all )  )
    class( res_obj )<- c( "list", "thresholdSCRNACountMatrix" )
    return( res_obj )

}


