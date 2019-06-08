#' @describeIn mast_filter plot the proportions of wells are filtered due to different criteria
#' @param byGroup in the case of \code{burdenOfFiltering} should the filter be stratified by groups, or only the plotting.
#' @examples
#' burdenOfFiltering(vbetaFA, groups='ncells', byGroup=TRUE)
#' burdenOfFiltering(vbetaFA, groups='ncells')
#'@export
burdenOfFiltering <- function(sc, groups, byGroup=FALSE, filt_control = NULL){
    checkGroups(sc, groups)
    conditionby <- NULL
    if(byGroup) conditionby <- groups
    filt <- mast_filter(sc, groups=conditionby, filt_control=filt_control, apply_filter=FALSE)
    index <- apply(filt, 1, function(x){suppressWarnings(mx <- min(which(x))); if(!is.finite(mx)) mx <- 4; mx})
    outcome <- c(names(filt), 'none')[index]
    filt <- cbind(as.data.frame(colData(sc)), outcome=outcome)
    tab <- do.call(table, c(filt[, groups, drop=FALSE], filt[, 'outcome', drop=FALSE]))

    lattice::barchart(tab, main=sprintf('Burden of filtering by %s', paste(groups, collapse=" and ")), auto.key=TRUE, horizontal=FALSE, stack=FALSE)
}

if(getRversion() >= "2.15.1") globalVariables('filter')

#'Concordance plots of filtered single vs n-cell assays
#'
#'Plot the average expression value of two subsets of the data.
#' Generally these might be 1 cell and multiple-cell replicates,
#' in which case if the \code{mcols} column \code{ncells} is
#' set then the averages will be adjusted accordingly.
#' But it could be any grouping.
#'@param SCellAssay is a FluidigmAssay for the 1-cell per well assay
#'@param NCellAssay is a FluidigmAssay for the n-cell per well assay
#'@param filterCriteria is a list of filtering criteria to apply to the SCellAssay and NCellAssay
#'@param groups is a character vector naming the group within which to perform filtering. NULL by default.
#' @param ... passed to \code{getConcordance}
#' @seealso getConcordance
#' @return printed plot
#'@examples
#' data(vbetaFA)
#' sca1 <- subset(vbetaFA, ncells==1)
#' sca100 <- subset(vbetaFA, ncells==100)
#' plotSCAConcordance(sca1, sca100)
#'@export
plotSCAConcordance<-function(SCellAssay, NCellAssay, filterCriteria=list(nOutlier = 2, sigmaContinuous = 9,sigmaProportion = 9), groups=NULL, ...){
    if(!(is(SCellAssay,"SingleCellAssay")&is(NCellAssay,"SingleCellAssay"))){
        stop("SCellAssay and NCellAssay must inherit from SingleCellAssay");
    }
    filtered.sc<-mast_filter(SCellAssay,filt_control=filterCriteria,groups=groups)
    filtered.hc<-mast_filter(NCellAssay,filt_control=filterCriteria, groups=groups)
    concord.unfiltered<-getConcordance(SCellAssay,NCellAssay, groups=groups, ...)
    concord.filtered<-getConcordance(filtered.sc,filtered.hc, groups=groups, ...)
    toplot<-rbind(concord.filtered,concord.unfiltered)
    toplot<-data.frame(toplot,filter=c(rep("filtered",nrow(concord.filtered)),rep("unfiltered",nrow(concord.unfiltered))))
    
    g<-ggplot2::ggplot(toplot)+
        ggplot2::geom_point(data=subset(toplot,filter=="filtered"),ggplot2::aes_string(x="et.ref",y="et.comp"),alpha=0.75)+
        ggplot2::theme_bw()+ggplot2::scale_x_continuous("Reference Assay")+ggplot2::scale_y_continuous("Comparison Assay")
    
    seg<-reshape2::dcast(melt(toplot,measure=c("et.ref","et.comp"),id=c("primerid","filter")),primerid~filter+variable)
    g<-g+ggplot2::geom_segment(data=seg,ggplot2::aes_string(x="unfiltered_et.ref",xend="filtered_et.ref",y="unfiltered_et.comp",yend="filtered_et.comp"),alpha=0.5,lwd=0.5)
    
    wss.filt<-getwss(concord.filtered,concord.filtered$nexp.ref)
    wss.unfilt<-getwss(concord.unfiltered,concord.unfiltered$nexp.ref)
    message(sprintf("Sum of Squares before Filtering: %s\n After filtering: %s\n Difference: %s",round(wss.unfilt,2),round(wss.filt,2),round(wss.unfilt-wss.filt,2)))
    g
}
