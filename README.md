SingleCellAssay
===============

Infrastructure and Tools for Single Cell Assay Analysis


Installation Instructions
------------
     install.packages('devtools')
     library(devtools)
     install_github('SingleCellAssay', 'RGLab')
     vignette('SingleCellAssay-intro')

**Needs roxygen2 4.0.0.99** or higher (available on github) if the documentation is regenerated--otherwise the NAMESPACE file will not be correct.

Changes
------------
New interface for zlm.SingleCellAssay to specify hypothesis.  See ?Hypothesis


New Features 
------------
- Support tests of arbitrary contrasts using LRT/zlm.SingleCellAssay

![doi/10.5281/zendoo.9810](http://zenodo.org/badge/doi/10.5281/zenodo.9810.png)
