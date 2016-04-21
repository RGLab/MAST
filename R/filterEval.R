#'Burden of filtering
#'
#'what proportions of wells are filtered due to different criteria
#'
#'@param sc SingleCellAssay or derived class
#'@param groups the groups by which to filter
#'@param byGroup logical indicating whether to filter by group
#'@param filt_control a list of control parameters.
#'@export
burdenOfFiltering <- function(sc, groups, byGroup=FALSE, filt_control = NULL){
  checkGroups(sc, groups)
  conditionby <- NULL
  if(byGroup) conditionby <- groups
  filt <- filter(sc, groups=conditionby, filt_control=filt_control, apply_filter=FALSE)
  index <- apply(filt, 1, function(x){suppressWarnings(mx <- min(which(x))); if(!is.finite(mx)) mx <- 4; mx})
  outcome <- c(names(filt), 'none')[index]
  filt <- cbind(as.data.frame(colData(sc)), outcome=outcome)
  tab <- do.call(table, c(filt[, groups, drop=FALSE], filt[, 'outcome', drop=FALSE]))

  lattice::barchart(tab, main=sprintf('Burden of filtering by %s', paste(groups, collapse=" and ")), auto.key=TRUE, horizontal=FALSE, stack=FALSE)
}

#'Concordance plots of filtered single vs n-cell assays
#'
#'This will generate the types of plots shown in the bioinformatics manuscript
#'@param SCellAssay is a FluidigmAssay for the 1-cell per well assay
#'@param NCellAssay is a FluidigmAssay for the n-cell per well assay
#'@param filterCriteria is a list of filtering criteria to apply to the SCellAssay and NCellAssay
#'@param groups is a character vector naming the group within which to perform filtering. NULL by default.
#'@export
plotSCAConcordance<-function(SCellAssay, NCellAssay, filterCriteria=list(nOutlier = 2, sigmaContinuous = 9,sigmaProportion = 9), groups=NULL){
  if(!(is(SCellAssay,"SingleCellAssay")&is(NCellAssay,"SingleCellAssay"))){
    stop("SCellAssay and NCellAssay must inherit from SingleCellAssay");
  }
  filtered.sc<-filter(SCellAssay,filt_control=filterCriteria,groups=groups)
  filtered.hc<-filter(NCellAssay,filt_control=filterCriteria, groups=groups)
  concord.unfiltered<-getConcordance(SCellAssay,NCellAssay)
  concord.filtered<-getConcordance(filtered.sc,filtered.hc)
  toplot<-rbind(concord.filtered,concord.unfiltered)
  toplot<-data.frame(toplot,filter=c(rep("filtered",nrow(concord.filtered)),rep("unfiltered",nrow(concord.unfiltered))))
  
  g<-ggplot2::ggplot(toplot)+
    ggplot2::geom_point(data=subset(toplot,filter=="filtered"),ggplot2::aes_string(x="et.ref",y="et.comp"),alpha=0.75)+
    ggplot2::theme_bw()+ggplot2::scale_x_continuous("Reference Assay")+ggplot2::scale_y_continuous("Comparison Assay")
 
  seg<-cast(melt(toplot,measure=c("et.ref","et.comp"),id=c("primerid","filter")),primerid~filter+variable)
  g<-g+ggplot2::geom_segment(data=seg,ggplot2::aes_string(x="unfiltered_et.ref",xend="filtered_et.ref",y="unfiltered_et.comp",yend="filtered_et.comp"),alpha=0.5,lwd=0.5)
  
  wss.filt<-getwss(concord.filtered,concord.filtered$nexp.ref)
  wss.unfilt<-getwss(concord.unfiltered,concord.unfiltered$nexp.ref)
  message(sprintf("Sum of Squares before Filtering: %s\n After filtering: %s\n Difference: %s",round(wss.unfilt,2),round(wss.filt,2),round(wss.unfilt-wss.filt,2)))
  g
}
