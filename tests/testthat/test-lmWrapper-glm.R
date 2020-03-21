obj <- new('GLMlike', design=colData(fd), formula=~Stim.Condition)
obj <- fit(obj, response=assay(fd)[2,])
objD <- glm(obj@response>0 ~ Stim.Condition, data=obj@design, family='binomial')
objC <- glm(obj@response ~ Stim.Condition, data=obj@design, subset=obj@response>0)
context('GLMlike')
source('common-lmWrapper-tests.R', local=TRUE)

context('GLM construction')


test_that('Design is invariant to updates', {
        obj2 <- update(obj, ~ Stim.Condition+ncells)
        expect_true(length(obj2@fitC)==0)
        expect_equivalent(obj2@design, obj@design)
        expect_true(length(setdiff(colnames(obj@modelMatrix), colnames(obj2@modelMatrix)))==0)
        expect_true(length(setdiff(colnames(obj2@modelMatrix), colnames(obj@modelMatrix)))==1)
})

source('common-lmWrapper-glm-tests.R', local=TRUE)
