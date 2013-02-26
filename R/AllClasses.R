##' @import Biobase
##' @import BiocGenerics
##' @importFrom plyr rbind.fill

NULL

##' Mapping class for SingleCellAssay package
##'
##' A class that represents a mapping of columns in a raw data file to cell-leve, feature-level, and phenotype-level metadata, as well as unique identifiers for individual cells.
##' mapNames for the SingleCellAssay class are in the object \code{SingleCellAssay:::SingleCellAssayMapNames}
##' mapNames for the FluidigmAssay class are in the object \code{SingleCellAssay:::FluidigmMapNames}
##' }
##' \section{Slots}{
##' \describe{
##' \item{mapping}{A named list providing the mapping from required fields in a SingleCellAssay or FluidigmAssay to column names in the data file}
##' }
##' 
##' @docType class
##' @name Mapping-class
##' @aliases Mapping
##' @aliases Mapping-class
##' @rdname Mapping-class
##' @title Column Mapping for SingleCellAssays
##' @seealso \code{\link{SingleCellAssay}},\code{\link{FluidigmAssay}}, \code{\link{getMapping}},\code{\link{addMapping}},\code{\link{getMapNames}},\code{\link{isEmpty}},\code{\link{mappingIntersection}}
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

##' Is this Mapping Empty
##'
##' Tests if a Mapping is empty
##' @param object The \code{Mapping}
##' @return A \code{logical}. \code{TRUE} if empty. \code{FALSE} otherwise.
##' @title iEmpty
##' @name isEmtpy
##' @docType methods
##' @rdname isEmpty-methods
##' @exportMethod isEmpty
##' @aliases isEmpty,Mapping-method
##' @aliases isEmpty
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
##' @exportMethod addMapping
##' @title addMapping
##' @name addMapping
setGeneric("addMapping",function(object,namedlist,...){standardGeneric("addMapping")})

##' @rdname addMapping-methods
##' @docType methods
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

setOldClass("ncdf")


setClass("SCA",
         representation=representation(env="environment",
           "VIRTUAL",mapping="Mapping",mapNames="character"),
         prototype=list(env=new.env(),mapping=new("Mapping"),mapNames=""))


##' Vbeta Data Set
##' @docType data
##' @name vbeta
##' @rdname vbeta-dataset
##' @format a data frame with 11 columns
NULL

### SingleCellAssay validity method
##'
##' Function to check the validity of SingleCellAssay objects. More specific functionality can wrap this function.
##' @export
##' @title SingleCellAssayValidity
##' @param object A \code{SingleCellAssay} object or object inheriting from SingleCellAssay
##' @return A \code{logical}. \code{TRUE} if the object is valid, \code{FALSE} otherwise.
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
  ##are the required mappings present...
  if(!all(object@mapNames%in%getMapNames(object)))
    return(FALSE)
  
  return(TRUE)
}


Mandatory_Fields <- c('primerid', 'geneid', 'measurement', 'idvars')
SingleCellAssayMapNames <- c(Mandatory_Fields, 'cellvars', 'featurevars', 'phenovars')
SingleCellAssayMap<-vector('list',length(SingleCellAssayMapNames))
names(SingleCellAssayMap)<-SingleCellAssayMapNames
                           

#Mapping Class Related Functions and Methods and the Class Definition
##'@export
.isValidNamedList<-function(mylist){
  #if the named list is empty, we don't add anything
  if(length(mylist)==0){
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



##' SingleCellAssay class
##' 
##' SingleCellAssay represents an arbitrary single cell assay
##' It is meant to be flexible and can (and will) be subclassed to represent specific assay
##' types like Fluidigm and others. It should be constructed using the \code{SingleCellAssay} constructor, or ideally the \code{SCASet} constructor.
##' mapNames for the SingleCellAssay class are in the object \code{SingleCellAssay:::SingleCellAssayMapNames}
##' mapNames for the FluidigmAssay class are in the object \code{SingleCellAssay:::FluidigmMapNames}
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
setClass("SingleCellAssay",contains="SCA",
         representation=representation(featureData="AnnotatedDataFrame",
           phenoData="AnnotatedDataFrame",
           cellData="AnnotatedDataFrame",
           description='data.frame',
           wellKey='character',
           id="ANY",                    #experiment descriptor
           env="environment"),		   
         prototype=prototype(phenoData=new("AnnotatedDataFrame"),
           featureData=new("AnnotatedDataFrame"),
           cellData=new("AnnotatedDataFrame"),
           description=data.frame(),
           wellKey=NA_character_,
           id=numeric(0),
           env=new.env()),validity=SingleCellAssayValidity)


## Same as SingleCellAssay, but with additional mapNames
FluidigmMapNames <- c(SingleCellAssayMapNames, 'ncells')
FluidigmMap<-vector('list',length(FluidigmMapNames))
names(FluidigmMap)<-FluidigmMapNames

##' @exportClass FluidigmAssay
setClass('FluidigmAssay', contains='SingleCellAssay', prototype=prototype(mapNames=FluidigmMapNames,mapping=new("Mapping",mapping=FluidigmMap)),validity=SingleCellAssayValidity)

##' @aliases getMapping,SingleCellAssay,missing-method
##' @aliases getMapping,SCA,missing-method
##' @rdname getMapping-methods
##' @export
setMethod("getMapping",c("SCA","missing"),function(object){
  return(object@mapping)
})

##' @aliases getMapping,SCA,character-method
##' @aliases getMapping,SingleCellAssay,character-method
##' @rdname getMapping-methods
##' @export
setMethod("getMapping",c("SCA","character"),function(object,mapnames){
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
##' @name SingleCellAssay-constructor
##' @rdname SingleCellAssay-constructor
##' @docType methods
##' @return SingleCellAssay object
##' @importFrom digest digest
##' @importFrom plyr ddply
##' @importFrom reshape expand.grid.df
SingleCellAssay<-function(dataframe=NULL,idvars=NULL,primerid=NULL,measurement=NULL,geneid=NULL,id=NULL, mapping=NULL, cellvars=NULL, featurevars=NULL, phenovars=NULL, ...){
  ### Add pheno key
  ### throw error if idvars isn't disjoint from geneid, probeid
  #if geneid == primerid make a new primerid column ensuring it is unique
  ## if(geneid==primerid||is.null(primerid)){
  ##   #creates a new column called primerid
  ##   mkunique<-function(x,G){
  ##     cbind(x,primerid=make.unique(as.character(get(G,x))))
  ##   }
  ##   dataframe<-ddply(dataframe,idvars,mkunique,G=geneid)
  ##   primerid<-"primerid"
  ## }else{
  ##   ddply(dataframe,idvars,mkunique,G=primerid)
  ## }
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
  
  cellCounts <- table(do.call(paste, dataframe[,getMapping(mapping,"idvars")[[1]], drop=FALSE]))
  incomplete <- !all(cellCounts == cellCounts[1])
  
  if(incomplete){
    message("dataframe appears incomplete, attempting to complete it with NAs")
    skeleton <- expand.grid.df(unique(dataframe[,getMapping(mapping,"featurevars")[[1]], drop=FALSE]), unique(dataframe[, getMapping(mapping,"cellvars")[[1]], drop=FALSE]))
    dataframe <- merge(skeleton, dataframe, all.x=TRUE, by=c(getMapping(mapping,"featurevars")[[1]], getMapping(mapping,"cellvars")[[1]]))
    cellCounts <- table(do.call(paste, dataframe[,getMapping(mapping,"idvars")[[1]]]))
  }

   primerCounts <- table(do.call(paste, dataframe[,getMapping(mapping,"primerid")[[1]], drop=FALSE]))
  if(!all(primerCounts==primerCounts[1])){
    stop('Some primers appear more often than others.  Either your data is incomplete or you have duplicate primerid')
  }
  
  ord <- do.call(order, dataframe[, c(getMapping(mapping,"primerid")[[1]], getMapping(mapping,"idvars")[[1]])])
  dataframe <- dataframe[ord,]
  assign("data",dataframe,envir=env)
  #wellKey <- seq_along(cellCounts)
  wellKey <- sapply(names(cellCounts),digest)
  env$data$`__wellKey` <- rep(wellKey, times=cellCounts[1])
  protoassay <- new("SingleCellAssay",env=env,mapping=mapping,id=id,wellKey=wellKey)
  
    cell.adf  <- new("AnnotatedDataFrame")
    pData(cell.adf)<-melt(protoassay)[1:nrow(protoassay), getMapping(mapping,"cellvars")[[1]], drop=FALSE]
    sampleNames(cell.adf) <- getwellKey(protoassay)

    ##pheno.adf <- new('AnnotatedDataFrame')
    ##need a phenokey into the melted data frame for this to make sense
    f.adf <- new('AnnotatedDataFrame')
    pData(f.adf) <- unique(melt(protoassay)[,getMapping(mapping,"featurevars")[[1]], drop=FALSE])
    sampleNames(f.adf) <- unique(melt(protoassay)[,getMapping(mapping,"primerid")[[1]]])
    protoassay@cellData<-cell.adf
    protoassay@featureData <- f.adf
    

  
  return(protoassay)
}

##' Constructor for a FluidigmAssay
##'
##' Constructs a FluidigmAssay object. Differs little from the SingleCellAssay constructor. Only the \code{ncells} parameter is additionally required.
##' Mapping argument has been removed for simplicity.
##' @title Fluidigm Assay Constructor
##' @param dataframe A data frame containing the raw data
##' @param idvars See \code{\link{SingleCellAssay}}
##' @param primerid See \code{\link{SingleCellAssay}}
##' @param measurement See \code{\link{SingleCellAssay}}
##' @param ncells A \code{character} specifying the column which gives the number of cells per well
##' @param geneid See \code{\link{SingleCellAssay}}
##' @param id An identifier for the resulting object. Should be a meaningful name.
##' @param cellvars See \code{\link{SingleCellAssay}}
##' @param featurevars See \code{\link{SingleCellAssay}}
##' @param phenovars See \code{\link{SingleCellAssay}}
##' @param ... Additional parameters passed to \code{SingleCellAssay} constructor
##' @return A FluidigmAssay object
##' @author Andrew McDavid and Greg Finak
##' @export FluidigmAssay
FluidigmAssay<-function(dataframe,idvars,primerid,measurement, ncells=NULL, geneid=NULL,id=NULL, cellvars=NULL, featurevars=NULL, phenovars=NULL, ...){
  
   mapping<-try(get("mapping",list(...)),silent=TRUE)
    if(inherits(mapping,"try-error")){
      #no mapping provided so construct one
      mapping<-new("Mapping",mapping=SingleCellAssay:::FluidigmMap)
    }
   
  if(!is.null(ncells)){
    mapping<-addMapping(mapping,list(ncells=ncells))
  }
  this.frame <- as.list(environment())
  ## Factor out code that builds the mapping
  this.frame$cellvars <- c(ncells, cellvars)
  sc <- do.call(SingleCellAssay, this.frame)
  as(sc, 'FluidigmAssay')
}

##' Constructs a SCASet
##'
##' An SCASet is a list of SingleCellAssays or objects inheriting from SingleCellAssay. The type of constructor called is determined by the value of contentClass, which should be the class of the SCA inheriting object contained in this SCASet. Both the class and the constructor should exist and have the same name. The code dynamically looks to see if the a function with the same name exists, and ASSUMES it is the constructor for the class.
##' @title SCASet constructor
##' @param dataframe flat data.frame ala SingleCellAssay
##' @param splitby either a character vector naming columns or a factor or a list of factors used to split dataframe into SingleCellAssays
##' @param idvars character vector naming columns that uniquely identify a cell
##' @param primerid character vector of length 1 that names the column that containing what feature was measured
##' @param measurement character vector of length 1 that names the column containing the measurement
##' @param contentClass a character, the name of the class being constructed within this SCASet. Defaults to SingleCellAssay. Other methods may pass in other classes, i.e. FluidigmAssay.
##' @param ... passed up to SingleCellAssay or other dynamically called constructor.
##' @return SCASet
##' @note The dynamic lookup of the constructor could be made more robust. 
##' @aliases SCASet
##' @rdname SCAset-methods
##' @export
##' @TODO SCASet constructor should perhaps take a SingleCellAssay class or FluidigmClass rather than a dataframe. Then we can learn the class type for construction.
SCASet<-function(dataframe,splitby,idvars=NULL,primerid=NULL,measurement=NULL,contentClass="SingleCellAssay",...){
  if(is.character(splitby) && all(splitby %in% names(dataframe))){
  spl<-split(dataframe,dataframe[, splitby])
} else if(is.factor(splitby) || is.list(splitby) || is.character(splitby)){
  spl <- split(dataframe, splitby)
} else{
  stop("Invalid 'splitby' specification")
}
  ## mapping<-try(get("mapping",list(...)),silent=TRUE)
  ## if(inherits(mapping,"try-error")){
  ##   #Mapping wasn't passed to the constructor.. 
  ## }else{
    
  ## }
  set<-vector("list",length(spl))
  names(set)<-names(spl)
  for(i in seq_along(set)){
    ##construct a call using contentClass
    F <- try(getFunction(contentClass),silent=TRUE)
    if(inherits(F,"try-error"))
      message("Can't construct a class of type ",contentClass[[1]],". Constructor of this name doesn't exist")
    cl<-as.call(list(as.name(contentClass[[1]]),dataframe=spl[[i]],idvars=idvars,primerid=primerid,id=names(spl)[[i]], measurement=measurement,...))
    set[[i]]<-eval(cl)
#    set[[i]]<-SingleCellAssay(dataframe=spl[[i]],idvars=idvars,primerid=primerid,id=names(spl)[[i]], measurement=measurement,...)
  }
  new("SCASet",set=set)
}




