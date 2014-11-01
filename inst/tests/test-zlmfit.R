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

context("Bootstrap")
test_that("Only return coef works", {
    zzinit2 <- zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2, onlyCoef=TRUE)
    expect_that(zzinit2, is_a('array'))
    expect_equal(dim(zzinit2)[1], ncol(fd2))
})

test_that("Bootstrap", {
    zf <- zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2)
    boot <- bootVcov1(zf, R=10)
    expect_is(boot, 'array')
    ## rep, genes, coef, comp
    expect_equal(dim(boot),c(10, dim(coef(zf, 'D')), 2))
})


