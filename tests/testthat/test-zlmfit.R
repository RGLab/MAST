fd2 <- fd[1:20,]

context("ZlmFit")
test_that('zlm.SingleCellAssay works', {
  zzinit <<- hushWarning(zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2, parallel=FALSE, method='glm', ebayes=FALSE), 'never estimible')
  expect_that(zzinit, is_a('ZlmFit'))
})

context("makeChiSqTable")
test_that("Zero dof gives unit P values", {
    ct <- makeChiSqTable(c(C=1, D=1), df=c(C=0, D=0), test='N')
    expect_equal(ct[,'Pr(>Chisq)'], c(1, 1, 1), check.attributes=FALSE)
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
    zz <- hushWarning(waldTest(zzinit, CoefficientHypothesis('PopulationVbetaResponsive:Stim.ConditionUnstim')), 'Some levels contain symbols')
    expect_equal(dim(zz)[1], 20)
 
    zz2 <- hushWarning(waldTest(zzinit, Hypothesis('`PopulationVbetaResponsive:Stim.ConditionUnstim`')), 'Some levels contain symbols')
    expect_equivalent(zz, zz2)
    expect_true(all.equal(zz, zz3, check.attributes=FALSE, tol=1e-5))



})

context('Log fold change calcs')
test_that('log fold changes match zero-inflated regression', {
    zzsimple <<- zlm.SingleCellAssay( ~ Stim.Condition, fd2, method='glm', ebayes=FALSE)
    lfc <- logFC(zzsimple)
    expect_equal(nrow(lfc$logFC), nrow(zzsimple@sca))
    expect_true(all(lfc$varLogFC>0, na.rm=TRUE))
    zlfc <- lm(exprs(zzsimple@sca)~ .+0, data=as.data.frame(model.matrix(zzsimple@LMlike)))
    diff <- sum((coef(zlfc)[-1,]-t(lfc$logFC))^2, na.rm=TRUE)
    expect_lt(diff, 1e-6)

    vzlfc <- diag(vcov(zlfc))
    ## drop intercept terms
    vzlfcNoIntercept <- vzlfc[-seq(1, to=length(vzlfc), by=nrow(coef(zlfc)))]
    expect_gt(cor(as.numeric(t(lfc$var)), vzlfcNoIntercept, use='pairwise'), .85)    
})

test_that('log fold change via contrasts', {
    lfc <- hushWarning(logFC(zzinit, contrast1=Hypothesis('`(Intercept)`+PopulationVbetaResponsive + `PopulationVbetaResponsive:Stim.ConditionUnstim`')), 'Some levels contain symbols')

    expect_equal(nrow(lfc$logFC), nrow(zzsimple@sca))
    expect_true(all(lfc$varLogFC>0, na.rm=TRUE))
    zlfc <- coef(lm(exprs(zzsimple@sca)~ .+0, data=as.data.frame(model.matrix(zzinit@LMlike))))
    lfc.lm <- colSums(zlfc[c('PopulationVbetaResponsive', '`PopulationVbetaResponsive:Stim.ConditionUnstim`'),])
    diff <- mean((lfc.lm   -lfc$logFC)^2/lfc.lm, na.rm=TRUE)
    expect_lt(diff, .04)
})



context('summary')
test_that('summary works', {
    zzs <- summary(zzsimple, logFC=TRUE, doLRT=FALSE)
    expect_is(zzs, 'summaryZlmFit')
    zzd <- zzs$datatable
    expect_equivalent(unique(as.character(zzd$contrast)), c('(Intercept)', 'Stim.ConditionUnstim'))
    expect_equivalent(sort(c('C', 'D', 'logFC', 'S')), sort(unique(as.character(zzd$component))))
    expect_equal(nrow(zzd), nrow(fd2)*(2*4-1)) #primerid, contrast (no intercept for logFC), component
    expect_output(print(zzs, n=2), c("Fitted zlm with top 2 genes per contrast:
( log fold change Z-score )
 primerid Stim.ConditionUnstim
 CCR7        7.3*             
 CD27        5.2*"), fixed=TRUE)
})
