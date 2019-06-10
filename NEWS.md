## 1.5.3 ##
Transition to use `SingleCellExperiment` and manage log counts slot.

### Version 1.7.3
* Fixed a bug were weighted residuals weren't being returned, causing the residuals hook to fail.
* Fixed documentation and vignettes.
* Added CovFromBoots() API to return inter-gene covariance matrices for model coefficients from GSEA bootstraps.
* Added some files and directories to .Rbuildignore to address notes and warnings.

## 1.11.2 ##

- Deprecated functions in version 1.8.0 are now Defunct.  `filter` (defunct since 1.8.0) has been removed.  Use mast_filter now.
- Refactored the `FromMatrix` constructor so that list-like assays don't make a roundtrip through an array.   This allows sparse-matrix and HDF5-backed assays to be provided to the constructor.
- zlm gains an `exprs_value` argument to allow the assay to be selected that will be tested.
