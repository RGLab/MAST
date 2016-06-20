ell <- function(x,xse,y,yse,segments=20,radius){
	center=c(x,y)
	angles <- (0:segments) * 2 * pi/segments
   unit.circle <- cbind(cos(angles), sin(angles))
   shape=matrix(c(xse^2,0,0,yse^2),ncol=2)
   Q <- chol(shape, pivot = TRUE)
   order <- order(attr(Q, "pivot"))
   ellipse <- t(center + radius * t(unit.circle %*% Q[, order]))
   colnames(ellipse)=c("x","y")
   ellipse
}

StatEll <- ggproto("StatEll", Stat,
  compute_group = function(data, scales,level=0.95,invert=FALSE) {
    e=ell(x=data$x, xse=data$xse,y=data$y,yse=data$yse,radius=sqrt(qchisq(level,df=2)))
    if(invert==c("x")){
    	e[,1] =arm::invlogit(e[,1])
    }else if(invert==c("y")){
    	e[,2] =arm::invlogit(e[,2])
    }
    return(e)
  },
  required_aes = c("x","xse", "y","yse","alpha")
)

#' Plot confidence ellipse in 2D.
#'
#' @param mapping Set of aesthetic mappings created by aes or aes_. If specified and inherit.aes = TRUE (the default), it is combined with the default mapping at the top level of the plot. You must supply mapping if there is no plot mapping.
#' @param data The data to be displayed in this layer. There are three options: If NULL, the default, the data is inherited from the plot data as specified in the call to ggplot. A data.frame, or other object, will override the plot data. All objects will be fortified to produce a data frame. See fortify for which variables will be created. A function will be called with a single argument, the plot data. The return value must be a data.frame., and will be used as the layer data.
#' @param geom The geometric object to use display the data
#' @param position Position adjustment, either as a string, or the result of a call to a position adjustment function.
#' @param na.rm If FALSE (the default), removes missing values with a warning. If TRUE silently removes missing values.
#' @param show.legend logical. Should this layer be included in the legends? NA, the default, includes if any aesthetics are mapped. FALSE never includes, and TRUE always includes.
#' @param inherit.aes If FALSE, overrides the default aesthetics, rather than combining with them. This is most useful for helper functions that define both data and aesthetics and shouldn't inherit behaviour from the default plot specification, e.g. borders.
#' @param fill A color or aesthetic mapping to fill color. Defaults to NA for empty ellipses.
#' @param level The confidence level at which to draw an ellipse (default is level=0.95).
#' @param lty The linetype to use. Can map to a variable. Defaults to 2 (dashed line)
#' @param invert vector of length 1 that should either be \code{"x","y",or TRUE}. Specifies whether to plot the estimates from the discrete component on the inverse logit scale. invert specifies which axis to invert.
#' @param ... other arguments passed on to layer. These are often aesthetics, used to set an aesthetic to a fixed value, like color = "red" or size = 3. They may also be parameters to the paired geom/stat.
#'
#' @export
#'
#' @importFrom ggplot2 layer
#' @importFrom ggplot2 ggproto
#' @importFrom ggplot2 Stat
stat_ell = function(mapping = NULL, data = NULL, geom = "polygon", position = "identity", na.rm = FALSE, show.legend = NA,inherit.aes = TRUE,fill=NA, level=0.95,lty=2,invert=FALSE,...) {
  ggplot2::layer(
    stat = StatEll, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm,level=level,lty=lty,fill=fill,invert=invert, ...)
  )
}


