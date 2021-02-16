
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

Vignettes are available in the package via `vignette('MAITAnalysis')`, `vignette('MAST-Intro')` or `vignette('MAST-interoperability')`.

New Features and announcements
------------
- MAST has been ported to use `SingleCellExperiment` under the hood, and is in [Bioconductor](http://bioconductor.org/packages/release/bioc/html/MAST.html).
- We now make an effort to track assay contents (counts vs log counts).  This should facilitate interaction with Scater and SCRAN.

Getting Help
----------------
For general questions, please submit a question to the [bioconductor support
site](https://support.bioconductor.org/?tag=MAST) so that others can
benefit from the discussion.

For bug reports (something seems broken): open a bug report [here](https://github.com/RGLab/MAST/issues).

Installation Instructions
------------
**This version available here on github may only properly function if you are running Bioconductor Devel, which is not something you will want to run for existing analyses! Instead follow instructions below.**

You may install or update `MAST` with:

    install.packages("BiocManager") # Needed to install all Bioconductor packages
    BiocManager::install("MAST")
    
Citation
----------------
If you find MAST useful in your work, please consider citing the
paper: [MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5)
G Finak, A McDavid, M Yajima, J Deng, V Gersuk, AK Shalek, CK Slichter
et al
Genome biology 16 (1), 278

The version that was used in the Genome Biology paper is accesible under the branch `MASTClassic`.


Converting old MASTClassic SingleCellAssay objects
--------

If you have data analyzed using MASTClassic,  you can convert
objects from MASTClassic format to the new format based on SingleCellExperiment using
`convertMastClassicToSingleCellAssay()`.

