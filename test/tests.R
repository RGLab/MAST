require(SingleCellAssay)
data(vbeta)
vbeta <- computeEtFromCt(vbeta)
vbeta.fa <- FluidigmAssay(vbeta, idvars = c("Subject.ID", "Chip.Number", "Well"),
  primerid = "Gene", measurement = "Et", ncells = "Number.of.Cells", geneid = "Gene",
  cellvars = c("Number.of.Cells", "Population"), phenovars = c("Stim.Condition",
    "Time"), id = "vbeta all")
show(vbeta.fa)

## cData replace method
cDat <- cData(vbeta.fa)
cDat$foo <- "bar"
cData(vbeta.fa) <- cDat
