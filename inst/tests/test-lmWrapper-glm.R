source('common-fixtures.R')
obj <- new('GLMlike', design=cData(fd), formula=~Stim.Condition+ncells)
source('common-lmWrapper-tests.R', local=TRUE)

context('GLM construction')
objglmD <- glm(obj@response>0 ~ .+0, data=as.data.frame(obj@modelMatrix), family='binomial')
objglmC <- glm(obj@response ~ .+0, data=as.data.frame(obj@modelMatrix), subset=obj@response>0)


test_that('Can get variance/cov', {
    expect_equivalent(vcov(obj, 'C'), vcov(objglmC))
    expect_equivalent(vcov(obj, 'D'), vcov(objglmD))
})

test_that('Design is invariant to updates', {
        obj2 <- update(obj, ~ Stim.Condition)
        expect_true(length(obj2@fitC)==0)
        expect_equivalent(obj2@design, obj@design)
        expect_true(length(setdiff(colnames(obj@modelMatrix), colnames(obj2@modelMatrix)))==1)
        expect_true(length(setdiff(colnames(obj2@modelMatrix), colnames(obj@modelMatrix)))==0)
})


