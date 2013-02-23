library(SingleCellAssay)

#develop a reader for nanostring data

datafiles<-list.files(pattern="RCC",path="/Users/gfinak/Documents/Projects/Nanostring/H9Cells/",recursive=TRUE,full=TRUE)
keyfiles<-"/Users/gfinak/Documents/Projects/Nanostring/h9key.csv"







process<-function(full=NULL){
#Date
full$Date<-as.Date(as.character(full$Date),format="%Y%m%d")

#Plate to Cell Cycle mapping
full<-merge(full,data.frame(Plate=factor(c(1,2,3,4,5,6)),CellCycle=c("G0/G1","G0/G1","S","S","G2/M","G2/M")))

#Number of cells is constant (1 cell per lane)
full$ncells<-1

#parse the comments in the sample info to get the cell type, plate, and row ID
full<-cbind(full,data.frame(CellType=substring(as.character(full$Comments),1,2),
                            Plate=substring(as.character(full$Comments),3,3),
                            RowID=substring(as.character(full$Comments),4,4)))
full
}
F<-NanoStringAssay(rccfiles=datafiles[1:10],keyfile=keyfiles,idvars=c("Plate","CartridgeID","ID"),geneid="GeneID",measurement="Count",phenoVars=c("CellType","CellCycle"),ncells="ncells",primerid="GeneID",cellvars=c("ncells","CellType"),id="nanostring",post.process.function=process)

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
##' @param ... Additional parameters passed to \code{SingleCellAssay} constructor
##' @return A FluidigmAssay object
##' @author Andrew McDavid and Greg Finak
##' @export NanoStringAssay
NanoStringAssay<-function(rccfiles=NULL,keyfile=NULL,idvars,primerid,measurement, ncells=NULL, geneid=NULL,id=NULL, cellvars=NULL, featurevars=NULL, phenovars=NULL, post.process.function=NULL,...){
  mapping<-new("Mapping",mapping=SingleCellAssay:::FluidigmMap)
  
  #transform the counts and rename the mapping for measurement to the new name
  #TODO fix the duplication of column names in a more reasonable way
  if(!is.null(ncells)){
    mapping<-addMapping(mapping,list(ncells=ncells))
  }
  
  rcclist<-readNanoStringLanes(rccfiles);
  rcclist<-mergeWithKeyFile(rcc=rcclist,keyfile)
  dataframe<-do.call(rbind,lapply(1:length(rcclist),function(i){
    cbind(rcclist[[i]]$code_summary,rcclist[[i]]$lane_attributes,rcclist[[i]]$sample_attributes)  
  }))
  if(!is.null(post.process.function)){
    dataframe<-post.process.function(dataframe) 
  }
  #reconstruct this for the next call
  #transform with log+1
  dataframe$lCount<-log(get(measurement,dataframe)+1)
  
  #reame the mapping 
  raw<-measurement
  measurement<-"lCount"
  this.frame <- as.list(environment())
  #remove unneeded variables
  this.frame$post.process.function<-NULL
  this.frame$rcclist<-NULL
  this.frame$rccfiles<-NULL
  this.frame$keyfile<-NULL
  
  ## Factor out code that builds the mapping
  this.frame$cellvars <- c(ncells, cellvars)
  sc <- do.call(SingleCellAssay, this.frame)
  sc@mapping<-addMapping(getMapping(sc),list(raw="Count"))
  #TODO need a cast for NanoStringAssay
  sc<-as(sc, 'FluidigmAssay')
  sc<-as(sc,'NanoStringAssay')
  return(sc)
}

setAs('FluidigmAssay', 'NanoStringAssay', function(from)  new("NanoStringAssay",env=from@env,mapping=addMapping(from@mapping,list(ncells=NULL)),id=from@id, wellKey=from@wellKey, featureData=from@featureData, phenoData=from@phenoData, cellData=from@cellData, description=from@description))

