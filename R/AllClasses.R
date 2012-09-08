
setOldClass("ncdf")

##' @import reshape
##' @import BiocGenerics
##' @importFrom plyr rbind.fill
##' @importClassesFrom flowCore ncdfHandler NcdfOrMatrix
setClass("SCA",
         representation=representation(exprs="NcdfOrMatrix",
           "VIRTUAL"),
         prototype=list(exprs=matrix(numeric(0),
                          nrow=0,
                          ncol=0)
           ))


##' Vbeta Data Set
##' @docType data
##' @name vbeta
##' @rdname vbeta-dataset
##' @format a data frame with 11 columns
NULL

### SingleCellAssay class
SingleCellAssayValidity <- function(object){
  if(!inherits(get("data",envir=object@env),"data.frame")){
    message("Argument `dataframe` should be a data.frame.")
    return(FALSE)
  }
  if(!any(getMapping(object,"idvars") %in% colnames(object@env$data))){
    warning("Invalid idvars column name. Not in data.frame")
    return(FALSE)
  }
  if(! ( getMapping(object,"geneid")%in%colnames(object@env$data))){
    warning("Invalid geneid column name. Not in data.frame")
    return(FALSE)
  }
  if(!(getMapping(object,"primerid")%in%colnames(object@env$data))){
    warning("Invalid primerid column name. Not in data.frame")
    return(FALSE)
  }
  if(!(getMapping(object,"measurement")%in%colnames(object@env$data))){
    warning("Invalid measurement column name. Not in data.frame")
    return(FALSE)
  }
  return(TRUE)
}


Mandatory_Fields <- c('primerid', 'geneid', 'measurement', 'idvars')
SingleCellAssayNames <- c(Mandatory_Fields, 'cellvars', 'featurevars', 'phenovars')
SingleCellAssayMap<-vector('list',length(SingleCellAssayNames))
names(SingleCellAssayMap)<-SingleCellAssayNames
                           

#Mapping Class Related Functions and Methods and the Class Definition
.isValidNamedList<-function(mylist){
  #if the named list is empty, we don't add anything
  if(length(mylist)==0){
    warning("namedlist should not be empty")
    return(TRUE)
  }
  if(is.null(names(mylist))){
    warning("namedlist must contain names")
    return(FALSE)
  }
  if(any(gsub(" ","",names(mylist))%in%"")){
    warning("namedlist names cannot be empty")
    return(FALSE)
  }
  if(!all(do.call(c,lapply(mylist,class))%in%c("character","NULL"))){
    warning("namedlist must be a named list of character vectors")
    return(FALSE)
  }
  return(TRUE)
}

setClass("Mapping",representation(mapping="list"),prototype=list(mapping=list("_empty_mapping_"=NULL)),validity=function(object){
  return(.isValidNamedList(object@mapping))
})

##' Accessor for mapNames
##'
##' This returns the mapNames, which are fields recognized in the names of mapping
##' @title getMapNames
##' @param object An object with a \code{Mapping} or \code{mapNames} slot
##' @return character vector of map names
##'
##' @exportMethod getMapNames
##' @docType methods
##' @rdname getMapNames-methods
##' @keywords accessor
setGeneric('getMapNames', function(object) standardGeneric('getMapNames'))

##' ##' @rdname getMapNames-methods
##' @aliases getMapNames,Mapping-method
setMethod("getMapNames","Mapping",function(object){
  return(names(object@mapping))
})

##' Accessor for mapping
##'
##' This returns the mapping, which is a named list.  The names of the list are keywords, while the contents of the list are column names in the melted \code{SingleCellAssay}.  The \code{mapping} identifies columns in the melted  which trigger special object behavior. Names recognized by the class are contained in \link[=getMapNames]{mapNames}, but additional names may be defined.
##' @title getMapping
##' @usage \code{getMapping(object)}
##' @usage \code{getMapping(object, mapnames)}
##' @param object A \code{Mapping} object or a \code{SingleCellAssay} or object inheriting from \code{SingleCellAssay}
##' @param mapnames A \code{character} vector of map names. Can be empty
##' @return A \code{Mapping} object if \code{object} is a \code{SingleCellAssay}. A character vector of mapped columns if \code{mapnames} is provided.
##' 
##' @export
##' @docType methods
##' @rdname getMapping-methods
##' @keywords accessor
setGeneric("getMapping",function(object,mapnames){standardGeneric("getMapping")})

##' @rdname getMapping-methods
##' @aliases getMapping,Mapping,missing-method
##' @export
setMethod("getMapping",c("Mapping","missing"),function(object){
  return(object@mapping)
})

##' @aliases getMapping,Mapping,character-method
##' @rdname getMapping-methods
##' @export
setMethod("getMapping",c("Mapping","character"),function(object,mapnames){
  return(object@mapping[mapnames])
})

setGeneric("isEmpty",function(object){standardGeneric("isEmpty")})
setMethod("isEmpty","Mapping",function(object){
  if(length(getMapNames(object))==0){
    return(TRUE)
  }
  if(any(getMapNames(object)%in%"_empty_mapping_")&length(getMapNames(object))==1){
    return(TRUE)
  }else{
    return(FALSE)
  }
 })

##' @rdname addMapping-methods
##' @docType methods
##' @export
##' @title addMapping
##' @name addMapping
setGeneric("addMapping",function(object,namedlist,...){standardGeneric("addMapping")})

##' @rdname addMapping-methods
##' @docType methods
##' @export
##' @aliases addMapping,Mapping,list-method
setMethod("addMapping",c("Mapping","list"),function(object,namedlist,replace=FALSE){
  if(!.isValidNamedList(namedlist)){
    stop("namedlist is not a valid named list for a Mapping")
  }
  if(isEmpty(object)){
    object@mapping<-namedlist
    return(object)
  }else if(!replace){
    inmapping<-getMapNames(object)
    idx<-!names(namedlist)%in%inmapping
    object@mapping<-c(object@mapping,namedlist[idx])
    append<-namedlist[!idx]
    for(nm in names(append)){
      object@mapping[[nm]]<-unique(c(object@mapping[[nm]],append[[nm]]))
    }
  }else{
    inmapping<-getMapNames(object)
    idx<-!names(namedlist)%in%inmapping
    object@mapping<-c(object@mapping,namedlist[idx])
    repl<-namedlist[!idx]
    for(nm in names(repl)){
      object@mapping[[nm]]<-repl[[nm]]
    }
    
  }
    return(object)
})

setGeneric("removeMapping",function(object,namedlist){standardGeneric("removeMapping")})
setMethod("removeMapping",c("Mapping","character"),function(object,namedlist){
  idx<-!getMapNames(object)%in%namedlist
  object@mapping<-object@mapping[idx]
  if(length(object@mapping)==0){
    object@mapping<-list("_empty_mapping_"=NULL)
  }
  return(object)
})

##' mappingIntersection checks the intersection of two maps in a \code{Mapping} object.
##' @param object A \code{Mapping} object
##' @param name1 An optional \code{character} map name
##' @param name2 Another optional \code{character} map name to be compared against \code{name1}
##' @return A matrix of intersecitons between different mappings if \code{name1} and \code{name2} are omitted, or a numeric value otherwise.
##' @usage \code{mappingIntersection(object)
##' @usage \code{mappingIntersection(object,name1,name2)
##' @rdname mappingIntersection-methods
##' @docType methods
##' @export
##' @title mappingIntersection
##' @name mappingIntersection-methods
##' @aliases mappingIntersection
setGeneric("mappingIntersection",function(object,name1,name2){standardGeneric("mappingIntersection")})

##' @rdname mappingIntersection-methods
##' @export
##' @aliases mappingIntersection,Mapping,missing,missing-method
setMethod("mappingIntersection",c("Mapping","missing","missing"),function(object){
  mapping<-getMapping(object)
  M<-matrix(0,nrow=length(mapping),ncol=length(mapping))
  colnames(M)<-names(mapping)
  rownames(M)<-names(mapping)
  for(i in 1:length(mapping)){
    for(j in 1:length(mapping)){
      M[i,j] <- length(intersect(mapping[[i]],mapping[[j]]))
    }
  }
  return(M)
})

##' @rdname mappingIntersection-methods
##' @aliases mappingIntersection,Mapping,character,character-method
##' @export
setMethod("mappingIntersection",c("Mapping","character","character"),function(object,name1,name2){
  return(length(intersect(getMapping(object,name1)[[1]],getMapping(object,name2)[[1]])))
})

##' show method for Mapping
##' @rdname show-methods
##' @aliases show,Mapping-method
setMethod("show","Mapping",function(object){
  if(isEmpty(object)){
    cat("An Empty Mapping\n")
  }else{
    cat("Mapping of size ",length(object@mapping),"\n")
    cat("With maps: ",getMapNames(object),"\n")
  }
})

##Not to be called directly by the user
##' SingleCellAssay represents an arbitrary single cell assay
##' It is meant to be flexible and can (and will) be subclassed to represent specific assay
##' types like Fluidigm and others. It should be constructed using the \code{SingleCellAssay} constructor, or ideally the \code{SCASet} constructor.
##' Required mapNames for the SingleCellAssay class are in the object \code{SingleCellAssayMapNames}
##' Required mapNames for the FluidigmAssay class are in the object \code{FluidigmMapNames}
##' }
##' \section{Slots}{
##' \describe{
##'   \item{featureData}{an \code{AnnotatedDataFrame} that describes feature-level metadata (i.e. genes)}
##'   \item{phenoData}{an \code{AnnotatedDataFrame} that describes the phenotype-level metadata (i.e. subject or experimental unit)}
##'   \item{cellData}{an \code{AnnotatedDataFrame} that describes the cell-level metadata (i.e. per individual cell)}
##'   \item{mapNames}{a \code{character} vector that describes some mandatory fields used by the class to map data from the raw file to the object. These are defined in the package, class definintion and subclasses.}
##'   \item{mapping}{a named \code{character} vector that maps \code{mapNames} to column names in the raw data file or data frame. This provides some flexibility for changing file formats and future assay types.}
##'   \item{description}{a \code{data.frame}}
##'   \item{wellKey}{A unique key that identifies the well, which USUALLY contains a single cell being measured.}
##'   \item{id}{a vector of type \code{character} that identifies the set of columns acting as a primary key to uniquely identify a single-cell or single-well across all wells / cells / assays / subjects / conditions in the data set.}
##'   \item{env}{an environment that will hold the data}
##' }
##' @name SingleCellAssay-class
##' @docType class 
##' @aliases SingleCellAssay-class
##' @aliases FluidigmAssay-class
##' @rdname SingleCellAssay-class
##' @exportClass SingleCellAssay
setClass("SingleCellAssay",
         representation=representation(featureData="AnnotatedDataFrame",
           phenoData="AnnotatedDataFrame",
           cellData="AnnotatedDataFrame",
           mapNames='character',
           mapping="Mapping",
           description='data.frame',
           cellKey='numeric',
           id="ANY",                    #experiment descriptor
           env="environment"),		   
         prototype=list(exprs=matrix(numeric(0),
                          nrow=0,
                          ncol=0),
           phenoData=new("AnnotatedDataFrame"),
           featureData=new("AnnotatedDataFrame"),
           cellData=new("AnnotatedDataFrame"),
           mapNames=SingleCellAssayNames,
           mapping=new("Mapping",mapping=SingleCellAssayMap),
           description=data.frame(),
           cellKey=numeric(0),
           id=numeric(0),
           env=new.env()),validity=SingleCellAssayValidity)


## Same as SingleCellAssay, but with additional mapNames
FluidigmMapNames <- c(SingleCellAssayNames, 'ncells')
FluidigmMap<-vector('list',length(FluidigmMapNames))
names(FluidigmMap)<-FluidigmMapNames

##' @exportClass FluidigmAssay
setClass('FluidigmAssay', contains='SingleCellAssay', prototype=prototype(mapNames=FluidigmMapNames,mapping=new("Mapping",mapping=FluidigmMap)))

##' @aliases getMapping,SingleCellAssay,missing-method
##' @rdname getMapping-methods
##' @export
setMethod("getMapping",c("SingleCellAssay","missing"),function(object){
  return(object@mapping)
})

##' @aliases getMapping,SingleCellAssay,character-method
##' @rdname getMapping-methods
##' @export
setMethod("getMapping",c("SingleCellAssay","character"),function(object,mapnames){
  mp<-getMapping(getMapping(object),mapnames)
  if(length(mp)>1){
    return(mp)
  }else{
    return(mp[[1]])
  }
})

##'SCASet is a set of SingleCellAssay objects or objects of its subclasses (i.e. FluidigmAssay)
##'The constructor \code{SCASet} should be used to make objects of this class.
##' }
##' \section{Slots}{
##' \describe{
##' \item{set}{A \code{list} of \code{SingleCellAssays} or its subclasses}
##' }
##' 
##' @rdname SCASet-class
##' @docType class
##' @name SCASet-class
##' @exportClass SCASet
##' @aliases SCASet-class
setClass("SCASet",
         representation=list(set="list"),validity=function(object){
           if(all(names(object@set)!=unlist(lapply(object@set,function(x) x@id),use.names=FALSE))){
             warning("Names of the SCASet don't match the SingleCellAssay id's. Plese use the SingleCellAssay() constructor.")
             return(FALSE)
           }
           return(TRUE)
         })


##' SingleCellAssay: A constructor for an object of type SingleCellAssay.
##'
##' This is the constructor for the class. This class intends to ease the analysis of single cell assays, in which multiple, exchangible, cells from an experimental unit (patient, or organism) are assayed along several (or many) dimensions, such as genes. A few examples of this might be Fluidigm gene expression chips, or single cell sequencing experiments.  The chief functionality is to make it easy to keep cellular-level metadata linked to the measurements through \code{cellData} and \code{phenoData}.  There are also subsetting and splitting measures to coerce between a SingleCellAssay, and a \link{SCASet}.
##' @param dataframe A 'flattened' data.frame containing columns giving cell and feature identifiers and  a measurement column
##' @param idvars character vector naming columns that uniquely identify a cell
##' @param primerid character vector of length 1 that names the column that identifies what feature (i.e. gene) was measured
##' @param measurement character vector of length 1 that names the column containing the measurement 
##' @param geneid character vector of length 1 that names a 'gene' column.  This could be placed into ...?
##' @param id optional numeric (not sure what this is supposed to do)
##' @param mapping named list.  Names are identifiers used by find special columns in the dataframe.  This shouldn't be modified directly by the user
##' @param cellvars Character vector naming columns containing additional cellular metadata
##' @param featurevars Character vector naming columns containing additional feature metadata
##' @param phenovars Character vector naming columns containing additional phenotype metadata
##' @param ... Additional keywords to be added to mapping
##' @export SingleCellAssay
##' @aliases SingleCellAssay
##' @name SingleCellAssay
##' @rdname SingleCellAssay-constructor
##' @docType methods
##' @return SingleCellAssay object
SingleCellAssay<-function(dataframe=NULL,idvars=NULL,primerid=NULL,measurement=NULL,geneid=NULL,id=NULL, mapping=NULL, cellvars=NULL, featurevars=NULL, phenovars=NULL, ...){
  ### Add pheno key
  ### throw error if idvars isn't disjoint from geneid, probeid

  env<-new.env()
  if(!is.null(mapping)){
    if(class(mapping)!="Mapping")
      stop("mapping argument must be a Mapping object")
  }else{
    mapping<-new("Mapping")
  }
  ## BEGIN: place into a validObject method!
  mapping.args<-list(idvars,cellvars,primerid,measurement,geneid,featurevars,phenovars)
  notnull<-unlist(lapply(mapping.args,function(x)!is.null(x)),use.names=FALSE)
  mapping<-addMapping(mapping,list(idvars=idvars,cellvars=cellvars,primerid=primerid,measurement=measurement,geneid=geneid,featurevars=featurevars,phenovars=phenovars)[notnull])
  if(! all(Mandatory_Fields %in% getMapNames(mapping)) )
      stop(paste('Mapping must contain at least ', paste(Mandatory_Fields, sep=', '), collapse=''))
  ## END: place into validObject method!
  #Update cellvars and featurevars with current mapping
  mapping<-addMapping(mapping,list(cellvars=unique(c(getMapping(mapping,"cellvars")[[1]],getMapping(mapping,"idvars")[[1]],getMapping(mapping,"phenovars")[[1]]))))
  mapping<-addMapping(mapping,list(featurevars=unique(c(getMapping(mapping,"featurevars")[[1]],getMapping(mapping,"primerid")[[1]],getMapping(mapping,"geneid")[[1]]))))
  ## mapping <- within(mapping, {
  ##        cellvars <- unique(c(cellvars, idvars, phenovars))
  ##        featurevars <- unique(c(featurevars, primerid, geneid))
  ##        })
  
    if(mappingIntersection(mapping,"cellvars","featurevars")>0)
      stop("'cellvars', 'idvars' must be disjoint from 'featurevars', 'primerid', 'geneid'")
    if(any(c(mappingIntersection(mapping,"phenovars","featurevars"),mappingIntersection(mapping,"phenovars","primerid"),mappingIntersection(mapping,"phenovars","geneid"))>0))
      stop("'phenovars' must be disjoint from 'featurevars', 'primerid', 'geneid'")
   if(nrow(unique(dataframe[, getMapping(mapping,"featurevars")[[1]], drop=FALSE])) != nrow((unique(dataframe[, getMapping(mapping,"primerid")[[1]], drop=FALSE]))))
       stop("'featurevars' must be keyed by 'primerid'")
   if(nrow(unique(dataframe[, getMapping(mapping,"cellvars")[[1]], drop=FALSE])) != nrow((unique(dataframe[, getMapping(mapping,"idvars")[[1]], drop=FALSE]))))
       stop("'cellvars' must be keyed by 'idvars'")
  
  ##check if idvars exists in dataframe
  ##check if probeid exists in dataframe
  ##check if geneid exists in dataframe
  ##check if measurement exists in dataframe
  
  cellCounts <- table(do.call(paste, dataframe[,getMapping(mapping,"idvars")[[1]]]))
  incomplete <- !all(cellCounts == cellCounts[1])
  
  if(incomplete){
    message("dataframe appears incomplete, attempting to complete it with NAs")
    skeleton <- expand.grid.df(unique(dataframe[,getMapping(mapping,"featurevars")[[1]], drop=FALSE]), unique(dataframe[, getMapping(mapping,"cellvars")[[1]], drop=FALSE]))
    dataframe <- merge(skeleton, dataframe, all.x=TRUE, by=c(getMapping(mapping,"featurevars")[[1]], getMapping(mapping,"cellvars")[[1]]))
    cellCounts <- table(do.call(paste, dataframe[,getMapping(mapping,"idvars")[[1]]]))
  }
  
  ord <- do.call(order, dataframe[, c(getMapping(mapping,"primerid")[[1]], getMapping(mapping,"idvars")[[1]])])
  dataframe <- dataframe[ord,]
  assign("data",dataframe,envir=env)
  cellKey <- seq_along(cellCounts)
  env$data$`__cellKey` <- rep(cellKey, times=cellCounts[1])
  protoassay <- new("SingleCellAssay",env=env,mapping=mapping,id=id,cellKey=cellKey)
  
    cell.adf  <- new("AnnotatedDataFrame")
    pData(cell.adf)<-melt(protoassay)[1:nrow(protoassay), getMapping(mapping,"cellvars")[[1]], drop=FALSE]
    sampleNames(cell.adf) <- getcellKey(protoassay)

    ##pheno.adf <- new('AnnotatedDataFrame')
    ##need a phenokey into the melted data frame for this to make sense
    f.adf <- new('AnnotatedDataFrame')
    pData(f.adf) <- unique(melt(protoassay)[,getMapping(mapping,"featurevars")[[1]], drop=FALSE])
    sampleNames(f.adf) <- unique(melt(protoassay)[,getMapping(mapping,"primerid")[[1]]])
    protoassay@cellData<-cell.adf
    protoassay@featureData <- f.adf
    

  
  return(protoassay)
}

FluidigmAssay<-function(dataframe,idvars,primerid,measurement, ncells, geneid=NULL,id=NULL, mapping=NULL, cellvars=NULL, featurevars=NULL, phenovars=NULL, ...){
  ## Factor out code that builds the mapping
  this.frame <- as.list(environment())
  this.frame$cellvars <- c(ncells, cellvars)
  sc <- do.call(SingleCellAssay, this.frame)
  as(sc, 'FluidigmAssay')
}

##' Constructs a SCASet, which is a list of SingleCellAssays
##'
##' FIXME
##' @title SCASet constructor
##' @param dataframe flat data.frame ala SingleCellAssay
##' @param splitby either a character vector naming columns or a factor or a list of factors used to split dataframe into SingleCellAssays
##' @param idvars character vector naming columns that uniquely identify a cell
##' @param primerid character vector of length 1 that names the column that containing what feature was measured
##' @param measurement character vector of length 1 that names the column containing the measurement
##' @param ... passed up to SingleCellAssay
##' @return SCASet
##' @aliases SCASet
##' @rdname SCAset-methods
##' @export
SCASet<-function(dataframe,splitby,idvars=NULL,primerid=NULL,measurement=NULL,...){
  if(is.character(splitby) && all(splitby %in% names(dataframe))){
  spl<-split(dataframe,dataframe[, splitby])
} else if(is.factor(splitby) || is.list(splitby) || is.character(splitby)){
  spl <- split(dataframe, splitby)
} else{
  stop("Invalid 'splitby' specification")
}

  set<-vector("list",length(spl))
  names(set)<-names(spl)
  for(i in seq_along(set)){
    set[[i]]<-SingleCellAssay(dataframe=spl[[i]],idvars=idvars,primerid=primerid,id=names(spl)[[i]], measurement=measurement,...)
  }
  new("SCASet",set=set)
}




