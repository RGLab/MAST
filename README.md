MAST: Model-based Analysis of Single-cell Transcriptomics
===============
MAST fits two-part, generalized linear models that are specially adapted for bimodal and/or zero-inflated single cell gene expression data.

Examples and vignettes
------------
MAST supports:

*  Easy importing, subsetting and manipulation of expression matrices
*  Filtering of low-quality cells
*  Adaptive thresholding of background noise
*  Tests for univariate differential expression, with adjustment for covariates
*  Gene set enrichment analysis, corrected for covariates and gene-gene correlations
*  Exploration of gene-gene correlations and co-expression

A [range of examples](https://github.com/RGLab/MASTdata/tree/master/vignettes) are available from the paper published in [Genome Biology](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5).
A vignette (in need of updating) is also available in the package via `vignette('MAST-intro')`.

New Features and announcements
------------
- MAST will soon be ported to use `SummarizedExperiment` under the hood, and submitted to Bioconductor.  See the summarizedExpt branch.

Installation Instructions
------------
If you have previously installed the package `SingleCellAssay` you will want to remove it as `MAST` supercedes `SingleCellAssay`.  (If both `MAST` and `SingleCellAssay` are attached, opaque S4 dispatch errors will result.)  Remove it with:

     remove.packages('SingleCellAssay')

Then you may install or update `MAST` with:

     install.packages('devtools')
     library(devtools)
     install_github('RGLab/MAST')
     # *or* if you don't have a working latex setup
     install_github(RGLab/'MAST', build_vignettes=FALSE)
     vignette('MAST-intro')



![doi/10.5281/zendoo.9810](http://zenodo.org/badge/doi/10.5281/zenodo.9810.png)
