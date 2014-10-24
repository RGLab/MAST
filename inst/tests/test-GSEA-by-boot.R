context('GSEA')
data(vbetaFA)
vb1 = subset(vbetaFA, ncells==1)
vb1 = vb1[,freq(vb1)>.1]
zf = zlm.SingleCellAssay(~Stim.Condition, vb1)
boots = bootVcov1(zf, 10)
sets=list(A=1:5, B=3:10, C=15, D=1:5)
gsea=gseaAfterBoot(zf, boots, sets, CoefficientHypothesis('Stim.ConditionUnstim'))
calcZ(gsea)
test_that('equal sets yield equal results', {
    expect_equal(gsea['A',,,],gsea['D',,,])
})

test_that('Singletons agree with coefficients', {
    expect_equal(gsea['C','cont','stat','test'], coef(zf, 'C')[15,'Stim.ConditionUnstim'])
    expect_equal(gsea['C','disc','stat','test'], coef(zf, 'D')[15,'Stim.ConditionUnstim'])
})
