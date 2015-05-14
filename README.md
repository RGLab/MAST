MAST
===============

Model-based Analysis of Single-cell Transcriptomics


Installation Instructions
------------
     install.packages('devtools')
     library(devtools)
     install_github('MAST', 'RGLab')
     # *or* if you don't have a working latex setup
     install_github('MAST', 'RGLab', build_vignettes=FALSE)
     vignette('MAST-intro')


New Features 
------------
- `gseaAfterBoot` for competitive geneset analysis under variance inflation
- Support tests of arbitrary contrasts using LRT/zlm.SingleCellAssay

![doi/10.5281/zendoo.9810](http://zenodo.org/badge/doi/10.5281/zenodo.9810.png)
