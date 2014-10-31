context('test prior construction')
test_that('default prior', {
    dp <- .defaultPrior(character(0))
    expect_equal(length(dp),0)
    dp <- .defaultPrior('a')
    dpn <- dimnames(dp)
    expect_equal(dpn[[1]],c('loc', 'scale', 'df'))
})


source('common-fixtures.R')
 obj <- new('BayesGLMlike', design=cData(fd), formula=~Stim.Condition + Population)
obj@coefPrior <- .defaultPrior(colnames(model.matrix(obj)))
test_that('update prior with model', {
    expect_equal(dim(obj@coefPrior)[3], ncol(model.matrix(obj))-1)
    model.matrix(obj) <- model.matrix(obj)[,-4]
    expect_equal(dim(obj@coefPrior)[3], ncol(model.matrix(obj))-1)
})

test_that('Strong prior changes estimates', {
    obj <- fit(obj, response=exprs(fd)[,2])
    obj@coefPrior['loc',,] <- 3
    obj@coefPrior['scale',,] <- .1
    obj@coefPrior['df',,] <- Inf
    obj2 <- fit(obj)
    diff <- coef(obj2, 'D')-coef(obj, 'D')
    expect_true(all(diff[names(diff) != '(Intercept)']>0))
    expect_true(diff['(Intercept)']<0)
})

