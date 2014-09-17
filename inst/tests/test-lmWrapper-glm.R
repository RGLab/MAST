source('common-fixtures.R')
obj <- new('GLMlike', design=cData(fd), formula=~Stim.Condition)
obj <- fit(obj, response=exprs(fd)[,2])
objD <- glm(obj@response>0 ~ Stim.Condition, data=obj@design, family='binomial')
objC <- glm(obj@response ~ Stim.Condition, data=obj@design, subset=obj@response>0)

source('common-lmWrapper-tests.R', local=TRUE)

context('GLM construction')


test_that('Design is invariant to updates', {
        obj2 <- update(obj, ~ Stim.Condition+ncells)
        expect_true(length(obj2@fitC)==0)
        expect_equivalent(obj2@design, obj@design)
        expect_true(length(setdiff(colnames(obj@modelMatrix), colnames(obj2@modelMatrix)))==0)
        expect_true(length(setdiff(colnames(obj2@modelMatrix), colnames(obj@modelMatrix)))==1)
})

test_that('Ebayes shrinkage corner cases', {
    ## zero DOF
    obj@priorDOF <- 0
    obj <- fit(obj)
    expect_equal(obj@fitC$dispersion, summary(objC)$dispersion)
    ## Complete shrinkage
    obj@priorDOF <- 1e8
    obj@priorVar <- 999
    obj <- fit(obj)
    expect_equal(obj@fitC$dispersion, 999, tol=1e-5)
    expect_equal(obj@fitC$dispersionNoShrink, summary(objC)$dispersion)
})
