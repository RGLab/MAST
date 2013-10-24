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
  expect_equal(coef(disc), coef(out$disc))
  expect_equal(coef(cont), coef(out$cont))
})

test_that('zlm accepts expressions in formulae', {
    out <- zlm(y ~ cut(x2, 3) + x1, dat)
})

if(require('lme4')){
  m <- melt(fd)
  m$Subject.ID <- factor(m$Subject.ID)
  m$Stim.Condition <- factor(m$Stim.Condition)
      lrout2 <- suppressWarnings(zlm(value ~ Population + (1|Subject.ID:Stim.Condition), data=m, lm.fun=glmer))
test_that('zlm can run lmer', {
  m$Subject.ID <- factor(m$Subject.ID)
  m$Stim.Condition <- factor(m$Stim.Condition)
      expect_is(lrout2$cont, c('mer','lmerMod','glmerMod'))
    expect_is(lrout2$disc, c('mer','lmerMod','glmerMod'))
})
    }

test_that('test.zlm works', {
    out <- zlm(y ~ x1 + x2, dat)
    test.zlm(out, car::matchCoefs(out$disc, 'x1'))

    if(require('lme4')){
      test.zlm(lrout2, car::matchCoefs(lrout2$disc, 'Population'))
    }
    
})

test_that("test.zlm handles 0-DoF cases gracefully", {
    out <- zlm(y~x1 + x2, dat[1:3,])
    var <- car::matchCoefs(out$disc, 'x1')
    suppressWarnings(test.zlm(out, var, type='LRT'))
})

  fd@keep.names <- FALSE
  fd2 <- fd[, 1:20]

test_that('zlm.SingleCellAssay works', {
  zz <- suppressWarnings(zlm.SingleCellAssay(value ~ Population*Stim.Condition, fd2, hypothesis.matrix='PopulationVbetaResponsive', .drop=TRUE))
  expect_that(zz, is_a('array'))
  expect_equal(dim(zz)[1], 20)
})

test_that('Guessing coefficients in zlm.SingleCellAssay works', {
    implicitCoefs <- zlm.SingleCellAssay(value ~ Subject.ID+Stim.Condition, fd2, hypo.terms=c('Subject.ID', 'Stim.Condition'), .drop=TRUE)
    explicitCoefs <- zlm.SingleCellAssay(value ~ Subject.ID+Stim.Condition, fd2, hypo.contrasts=c('Subject.IDSub02', 'Stim.ConditionUnstim'), .drop=TRUE, keep.zlm=TRUE)
    expect_equal(implicitCoefs, explicitCoefs$tests)

    ## TODO
    ## lrt <- zlm.SingleCellAssay(value ~ Subject.ID+Stim.Condition, fd2, hypothesis.matrix=c('Subject.ID', 'Stim.Condition'), type='LRT', .drop=TRUE, keep.zlm=TRUE)
    ## lrt1.d <- drop1(lrt$models[[1]]$disc, ~Subject.ID + Stim.Condition, test='Chisq')
    ## expect_equal(lrt$tests[1,'Chisq', 'disc'], lrt1.d[3, 'LRT'])
    
})

test_that("zlm.SingleCellAssay doesn't die on 100% expression", {
  ee <- exprs(fd2)
  ee[,1] <- rnorm(nrow(fd))+20
  exprs(fd2) <- ee
  zz <- zlm.SingleCellAssay(value ~ Population*Stim.Condition, fd2, hypothesis.matrix='PopulationVbetaResponsive', .drop=TRUE)
  expect_that(zz, is_a('array'))
  expect_equal(dim(zz)[1], 20)

  w.resp <- which(cData(fd2)$Population=='VbetaResponsive')
  ee[,1][w.resp] <- rbinom(length(w.resp), 1, .1)
  exprs(fd2) <- ee
  zz <- suppressWarnings(zlm.SingleCellAssay(value ~ Population, fd2, hypothesis.matrix='PopulationVbetaResponsive', .drop=TRUE, keep.zlm=TRUE))
  expect_that(zz$tests, is_a('array'))
  expect_equal(dim(zz$tests)[1], 20)  
})

if(require('lme4')) detach('package:lme4')
