[![Build Status](https://travis-ci.org/RGLab/MAST.svg?branch=master)](https://travis-ci.org/RGLab/MAST)

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


Vignettes are available in the package via `vignette('MAITAnalysis')` or `vignette('MAST-intro')`.

New Features and announcements
------------
- MAST has been ported to use `SummarizedExperiment` under the hood, and is in [Bioconductor](http://bioconductor.org/packages/release/bioc/html/MAST.html).
The main difference is that the data container is now transposed to follow bioconductor standards.
The older version will remain accessible under branch *MASTClassic*

Getting Help
----------------
For bug reports (something seems broken): open a bug report [here](https://github.com/RGLab/MAST/issues).  For
general questions, please submit a question to the [bioconductor support
site](https://support.bioconductor.org/t/MAST/) so that others can
benefit from the discussion.

Citation
----------------
If you find MAST useful in your work, please consider citing the
paper: [MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5)
G Finak, A McDavid, M Yajima, J Deng, V Gersuk, AK Shalek, CK Slichter
et al
Genome biology 16 (1), 278


Installation Instructions
------------
If you have previously installed the package `SingleCellAssay` you will want to remove it as `MAST` supercedes `SingleCellAssay`.  (If both `MAST` and `SingleCellAssay` are attached, opaque S4 dispatch errors will result.)  Remove it with:

     remove.packages('SingleCellAssay')

Then you may install or update `MAST` with:

    source("https://bioconductor.org/biocLite.R")
    biocLite("MAST")

Converting old MASTClassic SingleCellAssay objects
--------

If you have data analyzed using MASTClassic, starting with MAST package version 1.0.4 you can convert
objects from MASTClassic format to the new format based on SummarizedExperiment using
`convertMastClassicToSingleCellAssay()`.

