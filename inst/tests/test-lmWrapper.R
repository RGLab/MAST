source('common-fixtures.R')
T <- new('GLMlike', design=cData(fd), formula=~Stim.Condition+ncells)
context('Testing GLMLike construction')
test_that('Can construct GLMLike', {
    expect_is(T, 'GLMlike')
})

context('GLM construction')
T <- fit(T, response=exprs(fd)[,2])
TglmD <- glm(T@response>0 ~ .+0, data=as.data.frame(T@modelMatrix), family='binomial')
TglmC <- glm(T@response ~ .+0, data=as.data.frame(T@modelMatrix), subset=T@response>0)
test_that('Can fit', {
    expect_is(coefC(T), 'numeric')
    expect_is(coefD(T), 'numeric')    
})

test_that('Can get variance/cov', {
    expect_equivalent(vcovC(T), vcov(TglmC))
    expect_equivalent(vcovD(T), vcov(TglmD))
})

test_that('Design is invariant to updates', {
        T2 <- update(T, ~ Stim.Condition)
        expect_true(length(T2@fitC)==0)
        expect_equivalent(T2@design, T@design)
        expect_true(length(setdiff(colnames(T@modelMatrix), colnames(T2@modelMatrix)))==1)
        expect_true(length(setdiff(colnames(T2@modelMatrix), colnames(T@modelMatrix)))==0)
})


test_that('log likelihood is increasing in model complexity', {
    l1 <- logLik(T)
    T2 <- update(T, ~ Stim.Condition)
    T2 <- fit(T2)
    l0 <- logLik(T2)
    expect_true(l0<l1)
})

context('Post hoc testing')
test_that('LRT For Glm', {
 atest <- lrTest(T, 'Stim.Condition')
 expect_is(atest, 'data.frame')
 ## Test non-existant continuous/etc
})

test_that('Wald For Glm', {
 atest <- waldTest(T, 'Stim.ConditionUnstim')
 expect_is(atest, 'data.frame')
})

test_that('Handle 100% expression', {
    T2 <- fit(T, rnorm(nrow(fd))+20)
    expect_is(coefC(T2), 'numeric')
    expect_false(T2@fitted['D'])
})
