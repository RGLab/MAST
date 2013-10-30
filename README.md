Note
======
Single Cell Assay is not working with lmer 1.0
You can download an older version [here](http://cran.r-project.org/src/contrib/Archive/lme4/lme4_0.999999-2.tar.gz)
while we work on a fix.

SingleCellAssay
===============

Infrastructure and Tools for Single Cell Assay Analysis


Installation Instructions
------------
     install.packages('devtools')
     library(devtools)
     install_github('SingleCellAssay', 'RGLab')
     vignette('SingleCellAssay-intro')


New Features 
------------
- Migrated underlying data storage DataLayer
- Added parallel support for reading Nanostring RCC files using foreach and dopar
- Hurdle Model implemented in zlm.SingleCellAssay
- Thresholding support for NanoString in thresholdNanoString

Bug Fixes
----------
- Fixed indexing of SingleCellAssay with an empty i index using [[

