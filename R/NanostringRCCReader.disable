
#'RCC Nanostring Reader
#'
#'Reads RCC Nanostring Lane files
#'and constructs a SingleCellAssay object
#'
#'Open each lane file and read line by line constructing the different components of the file
#'
#'Metadata is present in the files. One lane per file.
#'@param x is a \code{character} vector of full path names to RCC files
#'@return a list
#' @importFrom stringr str_trim
#' @import foreach
#'@export
readNanoStringLanes<-function(x){
  nfiles<-length(x)
  #lanes<-vector("list",nfiles)a
  lanes<-suppressWarnings(foreach(i=1:length(x),.combine=c) %dopar% {
    R<-NULL
    this.file<-x[i]
    con<-suppressWarnings({
    capture.output(r<-fread(this.file,sep="\n",sep2=",",header=FALSE,autostart=1,verbose=FALSE))
    r}
    )
    #    con<-str_trim(readLines(this.file))
    R$header<-.readSection(con, this.file, 'Header', header=FALSE, as.is=TRUE)
    R$sample_attributes<-.readSection(con, this.file, 'Sample_Attributes', header=FALSE, as.is=TRUE)
    R$lane_attributes<-.readSection(con, this.file, 'Lane_Attributes', header=FALSE, as.is=TRUE)
    R$code_summary<-.readSection(con, this.file, 'Code_Summary', recast=FALSE, header=TRUE, as.is=TRUE)
    list(R)
  })
  return(lanes)
}



##' Read a section of a Nanostring RCC file
##'
##' The section is surrounded by <type> </type>
##' @rdname readSection
##' @param x character vector of lines from RCC file, newlines stripped
##' @param file name 
##' @param type the section of the file to be read in
##' @param recast should the section be transposed?
##' @param ... argument spassed to read.csv
##' @return data.frame of section
.readSection <- function(x, file, type, recast=TRUE, ...){
  #START <- sprintf('<%s>', type)
  #END <- sprintf('</%s>', type)
  #si <- match(START, x)
  #ei <- match(END, x)
  #if(is.na(si) || is.na(ei) || ei<si){
  #  stop('Malformed ', type, ' in ', file)
  #}
  wh<-which(x[,grepl(type,V1),])
  if(length(wh)!=2){
    stop("Malformed ", type, ' in ', file)
  }
  secind<-wh+c(1,-1)
  section<-x[eval(as.call(c(as.name("seq"),as.list(secind)))),str_split(V1,",")]
  if(recast){
    names<-section[1,]
    section<-section[2:nrow(section),]
    setnames(section,as.character(names))
  }else{
    names<-t(section[,1,with=FALSE])
    section<-data.table(t(section[,2:ncol(section),with=FALSE]))
    setnames(section,names)
  }
  #  section<-read.csv(textConnection(x[(si+1):(ei-1)]),...)
  #if(recast){
  #  wide <- data.frame(t(section), stringsAsFactors=FALSE)[-1,]
  #wide <- data.table(t(section))[-1L,]
  #setnames(wide,old=colnames(wide),new=t(section)[1,])
  #  names(wide) <- t(section)[1,]
  #} else{
  wide <- section
  wide
}

#'Merge NanoString lanes with a key file
#'
#'Reads in a key file and maps to nanostring lanes
#'@param rcc is a list returned from readNanoStringLanes
#'@param file is a character vector name of the key file
#'@return mapped rcc list
#'@importFrom reshape rename
#'@export
mergeWithKeyFile<-function(rcc,file){
  capture.output(f<-fread(file,header=F))
  key<-setnames(f,c("V1","V2"),c("Name","GeneID"))
  setkey(key,"Name")
  for(i in 1:length(rcc)){
    setkey(rcc[[i]]$code_summary,"Name")
  }
  #key<-data.frame(rename(read.csv(file,header=FALSE),c("V1"="Name","V2"="GeneID")))
  return(lapply(rcc,function(x){x$code_summary<-merge(x$code_summary,key,by="Name");x}))
}


#Could write a constructor that takes a post-processing function...
##' @exportClass NanoStringAssay
setClass('NanoStringAssay', contains='FluidigmAssay', prototype=prototype(mapNames=SingleCellAssay:::FluidigmMapNames,mapping=new("Mapping",mapping=SingleCellAssay:::FluidigmMap)),validity=SingleCellAssayValidity)


##' Constructor for a NanoStringAssay
##'
##' Constructs a NanoStringAssay object. Differs little from the FluidigmAssay constructor. Accepts a function for post-processing.
##' mapping argument has been removed as well for simplicity.
##' @title NanoString Assay Constructor
##' @param rccfiles A list of rcc files with full paths
##' @param keyfile A full path to a keyfile
##' @param idvars See \code{\link{SingleCellAssay}}
##' @param primerid See \code{\link{SingleCellAssay}}
##' @param measurement The measurement column for raw data. Will be placed in a variable named \code{raw} in the environment storing the data. The transformed and thresholded data will be placed in data.
##' @param ncells A \code{character} specifying the column which gives the number of cells per well
##' @param geneid See \code{\link{SingleCellAssay}}
##' @param id An identifier for the resulting object. Should be a meaningful name
##' @param cellvars See \code{\link{SingleCellAssay}}
##' @param featurevars See \code{\link{SingleCellAssay}}
##' @param phenovars See \code{\link{SingleCellAssay}}
##' @param post.process.function function applied to \code{data.frame} of all rcc files, before the NanostringAssay object is constructed.
##' @param ltrans Should the counts be log2 + 1  transformed? Default is FALSE. Set to TRUE when constructing from RAW data.
##' @param ... Additional parameters passed to \code{SingleCellAssay} constructor
##' @return A FluidigmAssay object
##' @author Andrew McDavid and Greg Finak
##' @import data.table
##' @export NanoStringAssay
NanoStringAssay<-function(rccfiles=NULL,keyfile=NULL,idvars,primerid,measurement, ncells=NULL, geneid=NULL,id=NULL, cellvars=NULL, featurevars=NULL, phenovars=NULL, post.process.function=NULL,ltrans=FALSE,...){
  
    mapping<-try(get("mapping",list(...)),silent=TRUE)
  if(inherits(mapping,"try-error")){
    mapping<-new("Mapping",mapping=SingleCellAssay:::FluidigmMap) 
  }
  #transform the counts and rename the mapping for measurement to the new name
  #TODO fix the duplication of column names in a more reasonable way
  if(!is.null(ncells)){
    mapping<-addMapping(mapping,list(ncells=ncells))
  }
  if(!is.null(rccfiles)){
    rcclist<-readNanoStringLanes(rccfiles);
    if(!is.null(keyfile)) rcclist<-mergeWithKeyFile(rcc=rcclist,keyfile)
    suppressMessages(
      for(i in 1:length(rcclist)){
        #cbind(rcclist[[i]]$code_summary,rcclist[[i]]$lane_attributes,rcclist[[i]]$sample_attributes, stringsAsFactors=FALSE)  
        rcclist[[i]]<-data.table(rcclist[[i]]$code_summary,rcclist[[i]]$lane_attributes,rcclist[[i]]$sample_attributes, stringsAsFactors=FALSE)  
      })
    dataframe<-do.call(data.table:::.rbind.data.table,rcclist)
    if(!is.null(post.process.function)){
      dataframe<-post.process.function(dataframe) 
    }
  }else if(inherits(dataframe<-get("dataframe",list(...)),"try-error")){
    stop("must provide either RCC files or a dataframe for NanoStringAssay")
  }

  #the explicit call to .rbind.data.table is necessary to resolve even though it's imported.
  ## dataframe <- lapply(dataframe, function(col){
  ##   suppressWarnings(numTry <- as.numeric(col))
  ##   ## non-numeric character vectors return NA
  ##   ## if all entries that weren't originally NA are not NA, then we suppose everything was a numeric
  ##   if( all(!is.na(numTry[!is.na(col)])) ){
  ##     return(numTry)
  ##   }
  ##   return(col)
  ## })

  #reconstruct this for the next call
  #transform with log+1
  #make sure the measurement column is numeric
  if(is.null(measurement)){
    measurement<-getMapping(mapping,"measurement")[[1]]
  }
  nummeas<-as.numeric(as.character(dataframe[,eval(as.name(substitute(measurement)),envir=.SD)]))
  dataframe[,eval(substitute(measurement)):=nummeas]
  #need to programmatically construct the transform and then apply it. ugh.
  #specify flag whether to transform
  if(ltrans){
    f<-(substitute(lCount:=log(eval(as.name(substitute(measurement)),envir=.SD)+1)))
    dataframe[,eval(f)]
    raw<-measurement
    measurement<-"lCount"
  }
  

  this.frame <- as.list(environment())
  #remove unneeded variables
  this.frame$post.process.function<-NULL
  this.frame$rcclist<-NULL
  this.frame$rccfiles<-NULL
  this.frame$keyfile<-NULL
  
  ## Factor out code that builds the mapping
  if(!is.null(ncells)&&!is.null(cellvars)){
    this.frame$cellvars <- c(ncells, cellvars)
  }
  sc <- do.call(SingleCellAssay, this.frame)
  if(is.null(getMapping(sc,"raw"))&ltrans){
    #best guess given the above
    sc@mapping<-addMapping(getMapping(sc),list(raw="Count"))
  }
  #TODO need a cast for NanoStringAssay
  sc<-as(sc, 'FluidigmAssay')
  sc<-as(sc,'NanoStringAssay')
  return(sc)
}

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
#' @return modifies nsa in place
NULL
## setMethod('threshold', signature='NanostringAssay', function(nsa, groups, thresholds, posteriorprob, clip=c('left', 'right', 'NA'){
##   exprs.all <- log2(melt(nsa)[, getMapping(nsa, 'raw')]+1)
##   out <- Mclust(exprs.all, G=1:3, modelNames='V')
##   if(out$G != 2)
##     stop("Uhoh, didn't find two clusters in complete data set")
##   means <- out$param$mean
##   scales <- out$param$variance$scale

##   exprs.split <- split(exprs.all, melt(nsa)[, c(groups, getMapping(nsa, 'primerid'))])

##   lapp <- lapply(exprs.split, function(x){
##     out <- Mclust(x, G=1:2, modelNames='V', prior=priorControl(shrinkage=.03, mean=means, scale=scales))
##     print(out$G)
##     if(out$G==1){
##       print(out$param$mean)
##     }
##     out
##   })

## })

setAs('FluidigmAssay', 'NanoStringAssay', function(from)  new("NanoStringAssay",env=from@env,mapping=addMapping(from@mapping,list(ncells=NULL)),id=from@id, wellKey=from@wellKey, featureData=from@featureData, phenoData=from@phenoData, cellData=from@cellData, description=from@description))
