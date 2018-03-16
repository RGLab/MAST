geneid="Gene"
primerid='Gene'
measurement='et'
idvars=c('Subject.ID', 'Chip.Number', 'Stim.Condition', 'Population', 'Well', 'Number.of.Cells')
phenovars=NULL
cellvars='Experiment.Number'
featurevars=NULL
ncells <- 'Number.of.Cells'

data(vbeta)
vbeta$et <- ifelse(is.na(vbeta$Ct), 0, 40-vbeta$Ct)
VBeta = vbeta
fd <- FromFlatDF(vbeta, idvars=idvars, primerid=primerid, measurement=measurement,cellvars=cellvars, geneid=geneid, ncells='Number.of.Cells', class='FluidigmAssay')
