\docType{methods}
\name{combine}
\alias{combine}
\alias{combine,SingleCellAssay,SingleCellAssay-method}
\title{Combine two SingleCellAssay or derived classes}
\description{
  combines two single cell assays provided they share a
  common mapping TODO combining based on a common mapping
  may be too restrictive. This may change depending on
  needs.

  Combines two Single Cell-like objects provided they have
  the same number of Features and Layers. The union of
  columns from featureData will be taken The union (padded
  if necessary with NA) will be taken from cellData.
}

