fd@keep.names <- FALSE
x <- matrix(runif(1000), 100)
colnames(x) <- paste('X', 1:10, sep='')
y <- x[, 1]*10 + rnorm(100)
y[y<median(y) & rbinom(100, prob=.8, size=1)>0] <- 0
y[y>median(y) & rbinom(100, prob=.1, size=1)>0] <- 0
context('testing zlm')
dat <- data.frame(y, x1=cut(x[,1], 4), x2=x[,2])
  disc <- glm(y>0 ~ x1 + x2, dat, family='binomial')
  cont <- lm(y ~ x1 + x2, subset=y>0, dat)

test_that('zlm throws meaningful error with matrix', {
  expect_error(zlm( y ~ x2, as.matrix(dat[, -2])), 'data.frame')
})

test_that('zlm throws error on NA', {
  dat2 <- dat
  dat2$y[1] <- NA
  expect_error(zlm( y ~ x2, dat2), 'NA')
})

test_that('zlm can run linear regression', {
  out <- zlm(y ~ x1 + x2, dat)
  expect_equivalent(coef(disc), coef(out$disc))
  expect_equivalent(coef(cont), coef(out$cont))
})

test_that('zlm accepts expressions in formulae', {
    out <- zlm(y ~ cut(x2, 3) + x1, dat)
})

  fd@keep.names <- FALSE
  fd2 <- fd[, 1:20]

if(require('lme4')){
  m <- melt(fd2)
  m$Subject.ID <- factor(m$Subject.ID)
  m$Stim.Condition <- factor(m$Stim.Condition)
test_that('zlm can run lmer', {
    lrout2 <- suppressWarnings(zlm(value ~ Population + (1|Subject.ID:Stim.Condition), data=m, method='lmer'))
      expect_is(lrout2$cont, c('mer','lmerMod','glmerMod'))
    expect_is(lrout2$disc, c('mer','lmerMod','glmerMod'))
})
    }
 

test_that('zlm.SingleCellAssay works', {
  zz <- suppressWarnings(zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2, hypothesis='PopulationVbetaResponsive'))
  expect_that(zz, is_a('array'))
  expect_equal(dim(zz)[1], 20)

  zz2 <- suppressWarnings(zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2, hypothesis='Population:Stim.Condition', type='LRT'))
  expect_that(zz, is_a('array'))
  expect_equal(dim(zz)[1], 20)
  
})

## test_that('Guessing coefficients in zlm.SingleCellAssay works', {
##     implicitCoefs <- zlm.SingleCellAssay(value ~ Subject.ID+Stim.Condition, fd2, hypo.terms=c('Subject.ID', 'Stim.Condition'), .drop=TRUE)
##     explicitCoefs <- zlm.SingleCellAssay(value ~ Subject.ID+Stim.Condition, fd2, hypo.contrasts=c('Subject.IDSub02', 'Stim.ConditionUnstim'), .drop=TRUE, keep.zlm=TRUE)
##     expect_equal(implicitCoefs, explicitCoefs$tests)

##     ## TODO
##     ## lrt <- zlm.SingleCellAssay(value ~ Subject.ID+Stim.Condition, fd2, hypothesis.matrix=c('Subject.ID', 'Stim.Condition'), type='LRT', .drop=TRUE, keep.zlm=TRUE)
##     ## lrt1.d <- drop1(lrt$models[[1]]$disc, ~Subject.ID + Stim.Condition, test='Chisq')
##     ## expect_equal(lrt$tests[1,'Chisq', 'disc'], lrt1.d[3, 'LRT'])
    
## })

test_that("zlm.SingleCellAssay doesn't die on 100% expression", {
  ee <- exprs(fd2)
  ee[,1] <- rnorm(nrow(fd))+20
  exprs(fd2) <- ee
  zz <- zlm.SingleCellAssay(~ Population*Stim.Condition, fd2, hypothesis='PopulationVbetaResponsive')
  expect_that(zz, is_a('array'))
  expect_equal(dim(zz)[1], 20)

  w.resp <- which(cData(fd2)$Population=='VbetaResponsive')
  ee[,1][w.resp] <- rbinom(length(w.resp), 1, .1)
  exprs(fd2) <- ee
  zz <- suppressWarnings(zlm.SingleCellAssay( ~ Population, fd2, hypothesis='PopulationVbetaResponsive'))
  expect_that(zz, is_a('array'))
  expect_equal(dim(zz)[1], 20)  
})

if(require('lme4')) detach('package:lme4')

context('Empirical Bayes')
if(require('numDeriv')){
test_that('Gradients match analytic', {
    set.seed(12345)
    rNg <- 101:200
    SSg <- rgamma(100, 3, 1)
    fn <- getMarginalHyperLikelihood(rNg, SSg)
    Grad <- getMarginalHyperLikelihood(rNg, SSg, deriv=TRUE)
    pts <- as.matrix(expand.grid(a0=c(.05, 1, 10), b0=c(.05, 1, 10)))
    for(i in nrow(pts))
        expect_equivalent(grad(fn, pts[i,]), Grad(pts[i,]))
})
}

test_that('Empirical Bayes works', {
     zz <- zlm.SingleCellAssay( ~ Population, fd2, hypothesis='PopulationVbetaResponsive', ebayes=TRUE)
})

context('Test error handling')
test_that('Give up after 5 errors', {
     expect_error(zlm.SingleCellAssay(value ~ Population*Stim.Condition, fd2, hypothesis='foo', force=FALSE), 'problems')

})

context('Test hooks')
test_that('Identity Hook', {
     zz <- zlm.SingleCellAssay(value ~ Population*Stim.Condition, fd2, hypothesis='Population:Stim.Condition', type='LRT', hook=function(x) x)
     expect_is(revealHook(zz)[[1]], 'GLMlike')
})

test_that('Residuals Hook', {
     zz <- zlm.SingleCellAssay(value ~ Population*Stim.Condition, fd2, hypothesis='Population:Stim.Condition', type='LRT', hook=residualsHook)
     fd3 <- collectResiduals(zz, fd2)
     expect_is(fd3, 'SingleCellAssay')
})

if(require('arm')){
context('zlm and bayesglm')
test_that('Can fit using bayesglm', {
    zz <- zlm.SingleCellAssay(~Population, fd2, hypothesis='Population', ebayes=FALSE, method='bayesglm', type='LRT', silent=FALSE)
    expect_is(zz, 'array')
})

test_that('Can do ebayes shrinkage using bayesglm', {
    zz <- zlm.SingleCellAssay(~Population, fd2, hypothesis='Population', ebayes=TRUE, method='bayesglm', type='LRT', silent=FALSE)
    expect_is(zz, 'array')
})

detach('package:arm')
}
