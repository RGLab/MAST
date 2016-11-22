#' Convert a MASTClassic SingleCellAssay
#' 
#' @description Convert a SingleCellAssay object created with the MASTClassic package to an object recognized by the new MAST package
#' @param object of class \code{SingleCellAssay} created by MASTClassic
#' @details The function will extract the relevant information from the attributes of the old object and construct a new
#' SingleCellAssay that is recognized by MAST. This function checks that the object is a MASTClassic SingleCellAssay object. It will stop 
#' if it is not a SingleCellAssay, return a converted SingleCellAssay if object was created by MASTClassic, and return the original
#' object if the object is already compatible.
#' @note Type checking for old object is not performed.
#' @return A MAST SingleCellAssay object. 
#' @export
#' @examples
#' data(vbetaFA)
#' convertMASTClassicToSingleCellAssay(vbetaFA)
convertMASTClassicToSingleCellAssay = function(object=NULL){
	if(is.null(object)){stop("object must not be NULL")	}
	if(!inherits(object,"SingleCellAssay")){stop("object must be a SingleCellAssay created by MASTClassic")}
	if("elementMetadata"%in%names(attributes(vbetaFA))){message("object does not appear to be created by MASTClassic. Returning original.");return(object)}
	attr(object,"class") = "array" #MASTClassic SingleCellAssay inherits from array
	nassays = dim(object)[3] #How many assay layers we want to create.
	phenodata = attr(object,"cellData") #this will be our phenoData
	featuredata = attr(object,"featureData") #these will be our row data
	message("Found ",nassays," layers in the old object.")
	dim_attr = attr(object,"dim")
	dimnames_attr = attr(object,"dimnames")
	attributes(object)=NULL
	attr(object,"dim")=dim_attr
	attr(object,"dimnames")=dimnames_attr
	object = aperm(object,c(2,1,3))
	DataFrame(phenodata)
	FromMatrix(exprsArray = object,cData = DataFrame(phenodata@data),fData = DataFrame(featuredata@data),class = "SingleCellAssay")
}