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

test_that('zlm can run linear regression', {
  out <- zlm(y ~ x1 + x2, dat)
  expect_equal(coef(disc), coef(out$disc))
  expect_equal(coef(cont), coef(out$cont))
})

if(require('lme4')){
  m <- melt(fd)
  m$Subject.ID <- factor(m$Subject.ID)
  m$Stim.Condition <- factor(m$Stim.Condition)
      lrout2 <- zlm(value ~ Population + (1|Subject.ID:Stim.Condition), data=m, lm.fun=glmer)
test_that('zlm can run lmer', {
  m$Subject.ID <- factor(m$Subject.ID)
  m$Stim.Condition <- factor(m$Stim.Condition)
      expect_is(lrout2$cont, 'mer')
    expect_is(lrout2$disc, 'mer')
})
    }

test_that('test.zlm works', {
    out <- zlm(y ~ x1 + x2, dat)
    test.zlm(out, matchCoefs(out$disc, 'x1'))

    if(require('lme4')){
      test.zlm(lrout2, matchCoefs(lrout2$disc, 'Population'))
    }
    
})

test_that('zlm.SingleCellAssay works', {
  zlm.SingleCellAssay(value ~ Population*Stim.Condition, fd, hypothesis.matrix='PopulationVbetaResponsive')
})
