fd@keep.names <- FALSE
fd2 <- fd[, 1:20]

context("ZlmFit")
test_that('zlm.SingleCellAssay works', {
  zzinit <<- zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2)
  expect_that(zzinit, is_a('ZlmFit'))
})

context("Likelihood ratio tests")
test_that('Three flavors of lrt on ZlmFit', {
    zz <- lrTest(zzinit, CoefficientHypothesis('PopulationVbetaResponsive:Stim.ConditionUnstim'))
    expect_equal(dim(zz)[1], 20)
 
    zz2 <- lrTest(zzinit, 'Population:Stim.Condition')
    CM <- matrix(0, nrow=ncol(zzinit@coefD), ncol=1, dimnames=list(colnames(zzinit@coefD), NULL))
    CM['PopulationVbetaResponsive:Stim.ConditionUnstim',] <- 1
    zz3 <- lrTest(zzinit, CM)
    expect_equivalent(zz, zz2)
    expect_true(all.equal(zz, zz3, check.attributes=FALSE, tol=1e-5))
})


context("Wald tests")
test_that('Three flavors of wald on ZlmFit', {
    CM <- matrix(0, nrow=ncol(zzinit@coefD), ncol=1, dimnames=list(colnames(zzinit@coefD), NULL))
    CM['PopulationVbetaResponsive:Stim.ConditionUnstim',] <- 1
    zz3 <- waldTest(zzinit, CM)
    zz <- waldTest(zzinit, CoefficientHypothesis('PopulationVbetaResponsive:Stim.ConditionUnstim'))
    expect_equal(dim(zz)[1], 20)
 
    zz2 <- waldTest(zzinit, Hypothesis('`PopulationVbetaResponsive:Stim.ConditionUnstim`'))
    expect_equivalent(zz, zz2)
    expect_true(all.equal(zz, zz3, check.attributes=FALSE, tol=1e-5))



})

context('Log fold change calcs')
test_that('log fold changes match zero-inflated regression', {
    zzsimple <<- zlm.SingleCellAssay( ~ Population, fd2)
    lfc <- logFC(zzsimple)
    expect_equal(nrow(lfc$logFC), ncol(zzsimple@sca))
    expect_true(all(lfc$varLogFC>0, na.rm=TRUE))
    zlfc <- lm(exprs(zzsimple@sca)~ .+0, data=as.data.frame(model.matrix(zzsimple@LMlike)))
    diff <- sum((coef(zlfc)[-1,]-t(lfc$logFC))^2, na.rm=TRUE)
    expect_less_than(diff, 1e-6)

    vzlfc <- diag(vcov(zlfc))
    ## drop intercept terms
    vzlfcNoIntercept <- vzlfc[-seq(1, to=length(vzlfc), by=nrow(coef(zlfc)))]
    expect_more_than(cor(as.numeric(t(lfc$var)), vzlfcNoIntercept, use='pairwise'), .85)    
})
