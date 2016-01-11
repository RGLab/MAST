fd@keep.names <- FALSE
set.seed(1234)
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
    skip_on_cran()
  m <- melt(fd2)
  m$Subject.ID <- factor(m$Subject.ID)
  m$Stim.Condition <- factor(m$Stim.Condition)
  test_that('zlm can run lmer', {
    lrout2 <- suppressWarnings(zlm(value ~ Population + (1|Subject.ID:Stim.Condition), data=m, method='lmer'))
    expect_is(lrout2$cont, c('mer','lmerMod','glmerMod'))
    expect_is(lrout2$disc, c('mer','lmerMod','glmerMod'))
})
    test_that('zlm.SingleCellAssay can run lmer', {
        z <- zlm.SingleCellAssay(~Population + (1|Subject.ID), fd2, method='lmer')
        expect_is(z, 'ZlmFit')
        expect_equal(nrow(z@df.null), 20)
        expect_equal(dim(z@vcovC)[[3]], 20)
    })
} else{
    message('Install lme4')
}
 

test_that('zlm.SingleCellAssay works', {
  zzinit <<- suppressWarnings(zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2))
  expect_that(zzinit, is_a('ZlmFit'))
  expect_equal(rownames(zzinit@coefC), fData(fd2)$primerid)
})

test_that("zlm.SingleCellAssay doesn't die on 100% expression", {
    fd3 <- fd2[,1:5]
  ee <- exprs(fd3)
  ee[,1] <- rnorm(nrow(fd3))+20
  exprs(fd3) <- ee
  zz <- zlm.SingleCellAssay( ~ Population, fd3)
    expect_that(zz, is_a('ZlmFit'))
  ## Really should be testing for linear separation or fitted values = 0 or 1 I guess..
  ##expect_less_than(zz@df.resid[1,'D'], 1)

        if(suppressPackageStartupMessages(require(arm))){
            zz3 <- zlm.SingleCellAssay( ~ Population, fd3, method='bayesglm')
            expect_that(zz3, is_a('ZlmFit'))
            expect_true(zz3@converged[1,'D'])
            detach('package:arm')
        } else{
            message('install arm')
        }

    w.resp <- which(cData(fd3)$Population=='VbetaResponsive')
    ee[,1] <- 0
    ee[,1][w.resp] <- rbinom(length(w.resp), 1, .2)
    exprs(fd3) <- ee
    zz2 <- zlm.SingleCellAssay( ~ Population, fd3)
    expect_that(zz2, is_a('ZlmFit'))
    expect_true(zz2@converged[1,'D'])
    
})

try(detach('package:lme4'), silent=TRUE)

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
} else{
    message('Install numDeriv')
}

test_that('Empirical Bayes works', {
     zz <- zlm.SingleCellAssay( ~ Population, fd2,ebayes=TRUE)
     expect_that(zz@dispersion, not(is_equivalent_to(zz@dispersionNoshrink)))
})

context('Test error handling')
test_that('Give up after 5 errors', {
     expect_error(zlm.SingleCellAssay(~ Population1234*Stim.Condition, fd2, force=FALSE), 'Population1234')
})

test_that('No holes in output', {
    ee <- exprs(fd2)
    ee[1,2] <- NA
    exprs(fd2) <- ee
    zze <- zlm.SingleCellAssay(~Stim.Condition, fd2)
    expect_equal(nrow(zze@coefD), ncol(fd2))
    expect_true(all(is.na(zze@coefD[2,])))
    expect_equal(dim(zze@vcovD)[3], ncol(fd2))
    expect_true(all(is.na(zze@vcovC[,,2])))
    expect_equal(nrow(zze@dispersion), ncol(fd2))
    expect_true(all(is.na(zze@dispersion[2,])))
})

context('Test hooks')
test_that('Identity Hook', {
     zz <- zlm.SingleCellAssay(value ~ Population, fd2, hook=function(x) x)
     expect_is(revealHook(zz)[[1]], 'GLMlike')
})

test_that('Residuals Hook', {
     zz <- zlm.SingleCellAssay(value ~ Population, fd2, hook=residualsHook)
     fd3 <- collectResiduals(zz, fd2)
     expect_is(fd3, 'SingleCellAssay')
})

if(suppressPackageStartupMessages(require('arm'))){
context('zlm and bayesglm')

test_that('Can fit using bayesglm', {
    zzinit <<- zlm.SingleCellAssay(~Population, fd2, ebayes=FALSE, method='bayesglm', silent=FALSE)
    expect_is(zzinit, 'ZlmFit')
})

test_that('Can do ebayes shrinkage using bayesglm', {
    zzinitshrink <- zlm.SingleCellAssay(~Population, fd2,  ebayes=TRUE, method='bayesglm', silent=FALSE)
    expect_that(zzinit@dispersion, not(is_equivalent_to(zzinitshrink@dispersion)))
    expect_equal(zzinit@dispersion, zzinitshrink@dispersionNoshrink)
})

try(detach('package:arm'), silent=TRUE)
}

