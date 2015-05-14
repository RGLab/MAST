
context('test prior construction')
test_that('default prior', {
    dp <- defaultPrior(character(0))
    expect_equal(length(dp),0)
    dp <- defaultPrior('a')
    dpn <- dimnames(dp)
    expect_equal(dpn[[1]],c('loc', 'scale', 'df'))
})


source('common-fixtures.R')
 obj <- new('BayesGLMlike', design=cData(fd), formula=~Stim.Condition + Population)
coefPrior <- defaultPrior(colnames(model.matrix(obj)))
obj@coefPrior <- coefPrior
test_that('update prior with model', {
    expect_equal(dim(obj@coefPrior)[3], ncol(model.matrix(obj))-1)
    model.matrix(obj) <- model.matrix(obj)[,-4]
    expect_equal(dim(obj@coefPrior)[3], ncol(model.matrix(obj))-1)
})

strongPrior <- coefPrior
strongPrior['loc',,] <- 3
strongPrior['scale',,] <- .1
strongPrior['df',,] <- Inf
test_that('Strong prior changes estimates', {
    obj <- fit(obj, response=exprs(fd)[,2])
    obj@coefPrior <- strongPrior
    obj2 <<- fit(obj)
    diff <- coef(obj2, 'D')-coef(obj, 'D')
    expect_true(all(diff[names(diff) != '(Intercept)']>0))
    expect_true(diff['(Intercept)']<0)
})

context('Fit args')
test_that('Same result with fit args', {
    obj3 <- new('BayesGLMlike', design=cData(fd), formula=~Stim.Condition + Population, fitArgsD=list(prior.mean=strongPrior['loc',1,], prior.scale=strongPrior['scale',1,], prior.df=strongPrior['df',1,]))
    obj4 <- fit(obj3, response=exprs(fd)[,2])
    expect_equal(obj3@fitArgsD, obj4@fitArgsD)
    expect_equal(coef(obj4, 'D'), coef(obj2, 'D'))
})

obj <- new('BayesGLMlike', design=cData(fd), formula=~Stim.Condition)
obj <- fit(obj, response=exprs(fd)[,2])
objC <- glm(obj@response ~ Stim.Condition, data=obj@design, subset=obj@response>0)

## Not really applicable since we've diverged from the arm codebase
## objD <- suppressWarnings(.bayesglm.fit(x=model.matrix(obj), y=obj@response>0, family=binomial()))
## ## bayesglm doesn't initialize this correctly, or consistently, in any case
## objD$aic <- deviance(objD)+2*objD$rank

## source('common-lmWrapper-tests.R', local=TRUE)


source('common-lmWrapper-glm-tests.R', local=TRUE)
