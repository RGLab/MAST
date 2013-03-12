SingleCellAssay
===============

Infrastructure and Tools for Single Cell Assay Analysis


New Features 
------------
- Migrated underlying data storage to data.table
- Added parallel support for reading Nanostring RCC files using foreach and dopar

Bug Fixes
----------
- Fixed indexing of SingleCellAssay with an empty i index using [[
- Fixed non-standard behaviour of SingleCellAssay indexing, foo[[TRUE]] would return the full set but now returns just the first well. To return the full set, do foo[[]] or foo[[1:nrow(foo)]]

