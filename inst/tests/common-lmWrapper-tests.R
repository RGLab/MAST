context('Testing LMLike construction')
test_that('Can construct LMLike', {
    expect_is(obj, 'LMlike')
})

obj <- fit(obj, response=exprs(fd)[,2])
test_that('Can fit', {
    expect_is(coef(obj, 'C'), 'numeric')
    expect_is(coef(obj, 'D'), 'numeric')    
})

test_that('Handle 0% expression', {
    obj2 <- fit(obj, rep(0, nrow(fd)))
    expect_false(any(obj2@fitted))
})

test_that('Handle Singular Designs', {
    obj2 <- update(obj, ~ . + Stim.Condition*Population)
    obj2 <- fit(obj2)
    expect_false(any(is.na(vcov(obj2, 'C'))))
    expect_false(any(is.na(vcov(obj2, 'D'))))
    expect_false(any(is.na(coef(obj2, which='C', singular=FALSE))))
    expect_false(any(is.na(coef(obj2, which='D', singular=FALSE))))
})

test_that('Handle 100% expression', {
    obj2 <- fit(obj, rnorm(nrow(fd))+20)
    expect_is(coef(obj2, 'C'), 'numeric')
    expect_false(obj2@fitted['D'])
})

test_that('Handle NA', {
    resp <- obj@response
    resp[1] <- NA
    expect_error(obj2 <- fit(obj, resp), 'NA')
})

## test_that('Handle expressions in formulae', {
##     obj2 <- update(obj, (~ . +cut(Experiment.Number, 2)))
##     obj2 <- fit(obj2)
##     expect_is(obj2, 'LMlike')
## })

context('Testing fit summaries')

test_that('log likelihood is increasing in model complexity', {
    l1 <- logLik(obj)
    obj2 <- update(obj, ~ . -Stim.Condition)
    obj2 <- fit(obj2)
    l0 <- logLik(obj2)
    expect_true(all(l0<=l1))
})



context('Post hoc testing')
test_that('LRT For Glm', {
 atest <- lrTest(obj, 'Stim.Condition')
 expect_is(atest, 'matrix')

 obj2 <- fit(update(obj, ~ .+Stim.Condition*Population))
 btest <- lrTest(obj2, 'Stim.Condition')
 expect_true(all(btest[,'df']==0))
 btest <- lrTest(obj2, 'Stim.Condition:Population')
 expect_equal(btest['cont','df'],1)
 
})

test_that('Wald For Glm', {
 atest <- waldTest(obj, 'Stim.ConditionUnstim')
 expect_is(atest, 'matrix')
})
