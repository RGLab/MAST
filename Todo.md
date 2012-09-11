Main remaining issues
=====================
*phenoData needs to be added via phenovars. 
*Testing that subsetting, combining and splitting do the right subset operations on the various AnnotatedDataFrames
*Rewrite the test suite using vbeta data. The public repository doesn't have the full set of data sets.

Class issues: Fixed
============
*Add mapping class.  This could be a good way to encapsulate what keys are supported in a mapping.
Done  
*SCASet constructor needs to be able to construct any subclass of SCA  
Done

The above issues are fixed as of 09/12/2012. I've added a mapping class and fixed the SCASet constructor as well as split,combine, etc.

SingleCellAssay issues
=======================
*SingleCellAssay could be refactored, it's getting a little long.   
This could still use some work.

*Validity of arguments is inconsistently verified.  :Fixed
I've made SingleCellAssay and FluidigmAssay inherit from SCA and all subclasses can now use the SingleCellAssayValidity function. It could be wrapped in another validity method for more specific checking.
*Constructors still do some validity checking. This could be improvedd.

*Add phenotype key to env$data  

*Rename cellData->wellData  
Half done. I've renamed cellKey to wellKey and updated the accessors for these as well as fixed import statements for cData and pData and fData etc.. I'm tempted to leave cellData named as is since it's consistent with the Biocbase generic method.


*Add replacement methods for cellData, featureData and phenoData.  Can someone just cbind additional columns onto cellData, featureData?  
Still to do

*Fix combine method.  Only the first two arguments are checked for alike mappings.  And maybe we don't need the mappings to totally align?  :Fixed
Combine works, i.e. checks for same mappings and also checks that the classes of the incoming objects is the same. I agree that we may not need mappings to totally align. We'll fix this as we see use cases.

*cellKey should be wellKey (since it identifies wells):Fixed
Done

SCASet
=======
*This might be a SingleCellAssay issue, but we would like to be able to construct SCASet of subclasses of SingleCellAssay by doing something   clever with an internal constructor or by post-construction coercion (likely slow). : Fixed
The above is now fixed. SCASet takes an additonal contentClass argument that is the name of the class it's meant to hold. It will check if the constructor for the class exists and call it dynamically. The same approach is used for combine and [[ subsetting of SCASet

*Add phenoData?   
Still need to do this. Not clear to me how to key this.. 

FluidigmAssay
=============
*Slow to construct because we first create a SingleCellAssay using its constructor, then pull out the slots

*No validity checking at the moment?  
Cleaned up the inheritance, and now it uses the default SingleCellAssayValidity function. The mapNames slots of these objects contain the "required" maps. The Mapping is all available maps. getMapNames returns all available maps. There is no accessor for @mapNames, since it doesn't really need public access. Default maps for subclasses of SCA are also available for construction and have been documented.

