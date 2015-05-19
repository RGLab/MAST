MAST
===============

Model-based Analysis of Single-cell Transcriptomics


Installation Instructions
------------
If you have previously installed the package `SingleCellAssay` you will want to remove it as `MAST` superscedes `SingleCellAssay`.  (If both `MAST` and `SingleCellAssay` are attached, opaque S4 dispatch errors will result.)  Remove it with:

     remove.packages('SingleCellAssay')

Then you may install or update `MAST` with:
     install.packages('devtools')
     library(devtools)
     install_github('RGLab/MAST')
     # *or* if you don't have a working latex setup
     install_github(RGLab/'MAST', build_vignettes=FALSE)
     vignette('MAST-intro')


New Features 
------------
- `gseaAfterBoot` for competitive geneset analysis under variance inflation
- Support tests of arbitrary contrasts using LRT/zlm.SingleCellAssay

![doi/10.5281/zendoo.9810](http://zenodo.org/badge/doi/10.5281/zenodo.9810.png)
