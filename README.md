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

New Features in version .85
------------
- The interface to zlm.SingleCellAssay has changed somewhat, and is not completely backwards compatible
    * But now it is easier to add hurdle models with different modeling functions.  Try method='glmer' or method='bayesglm' for mixed effects or bayesian methods.
    * Multiple hypothesis can be tested by passing a list of hypothesis (either terms to drop for Likelihood Ratio Tests, or hypotheses formated *a la* car::lht
- Thresholding support for NanoString in thresholdNanoString

![doi/10.5281/zendoo.9810](http://zenodo.org/badge/doi/10.5281/zenodo.9810.png)
