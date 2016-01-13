skip_on_cran()
getX <- function(groups, minsize, N){
    gid <- LETTERS[seq_len(groups)]
    fidMin <- rep(gid, each=minsize)
    fid <- sample(gid[-1], size=N-length(fidMin), replace=TRUE)
    fid <- data.frame(group=c(fid, fidMin))
    X <- model.matrix(~group, data=fid)
    structure(X, group=fid)
}


simYs <- function(m, X, beta, rho, sigma, p){
    ## Simulate a matrix of continuous expression values
    ## m is number of genes
    ## X is design
    ## beta is m * nrow(X) matrix of coefficients
    ## rho is covariance between cells
    ## sigma is variance per gene
    ## p is vector of marginal frequency of expression
    ## returns expression, design and errors
    
    ## genes-major order
    eps <- rnorm(nrow(X)*m)*sigma
    Z <- rnorm(nrow(X))*rho
    ## cells X genes
    err <- matrix(eps+rep(Z, each=m), nrow=nrow(X), byrow=TRUE)
    Y <- X %*% beta + err
    covX <- solve(crossprod(X))
    covErr <- cov(err)
    expr <- matrix(runif(m*nrow(X))<p, ncol=m, byrow=TRUE)
    list(Y=Y*expr, X=X, cov=kronecker(covErr, covX))
}

fd@keep.names <- FALSE
fd2 <- fd[, 1:20]

context("Bootstrap")
test_that("Only return coef works", {
    zzinit2 <- suppressWarnings(zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2, onlyCoef=TRUE))
    expect_that(zzinit2, is_a('array'))
    expect_equal(dim(zzinit2)[1], ncol(fd2))
})

cl <- parallel::makeForkCluster(2)
test_that("Bootstrap", {
    zf <- suppressWarnings(zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2))
    boot <- hushWarning(pbootVcov1(cl, zf, R=3),fixed('never estimible'))
    expect_is(boot, 'array')
    ## rep, genes, coef, comp
    expect_equal(dim(boot),c(3, dim(coef(zf, 'D')), 2))
})

context("Bootstrap consistency")
set.seed(1234)
N <- 200
m <- 20
middle <- floor(seq(from=m/3, to=2*m/3))
end <- floor(seq(from=2*m/3, m))
p <- 2
X <- getX(p, 40, N)
beta <- t(cbind(15, rep(3, m)))
pvec <- seq(.05, .95, length.out=m)
Y <- simYs(m, X, beta, rho=1, sigma=1, p=pvec)
cData <- data.frame(group=attr(X, 'group'))
sca <- suppressMessages(suppressWarnings(FromMatrix('SingleCellAssay', Y$Y, cData=cData)))
zfit <- suppressWarnings(zlm.SingleCellAssay(~group, sca=sca))
test_that('Expression frequencies are close to expectation', {
    expect_less_than(mean((freq(sca)-pvec)^2), 1/(sqrt(N)*m))
})

test_that('Discrete group coefficient is close to zero', {
    expect_less_than(
        abs(mean(coef(zfit, 'D')[middle,'groupB'], na.rm=TRUE)),
        10/(sqrt(N))
        )
})

test_that('Continuous group coefficient is close to expected', {
    expect_less_than(
        mean((coef(zfit, 'C')[end,'groupB']-beta[2,end])^2, na.rm=TRUE),
        3.5*Y$cov[2,2] #expected covariance of groupB
        )
})
clusterEvalQ(cl, set.seed(1234))
boot <- pbootVcov1(cl, zfit, R=50)
bootmeans <- colMeans(boot, na.rm=TRUE, dims=1)
## m2 <- 4  #top 4 expressed genes in simulation
## end4 <- (m-m2+1):m
## sca4 <- sca[,end4]
## sca4 <- subset(sca4, rowMeans(exprs(sca4)==0)==0)
## mlm <- lm(exprs(sca4) ~ group, cData(sca4))

test_that('Bootstrap is unbiased', {
    expect_less_than(mean((bootmeans[,,'C']-coef(zfit, 'C'))^2), 3*sum(abs(Y$cov[2,2])))
})

covInterceptC <- cov(boot[,,'(Intercept)','C'], use='pairwise')
## expectedCovInterceptC <- vcov(mlm)[seq(1, p*(m2), by=p), seq(1, p*(m2), by=p)]
expectedCovInterceptC <- Y$cov[seq(1, p*m, by=p), seq(1, p*m, by=p)]

test_that('Bootstrap recovers covariance', {
    sub <- covInterceptC[end,end]
    esub <- expectedCovInterceptC[end,end]
    ## approximately 40% tolerance
    expect_less_than(abs(log(mean(sub[upper.tri(sub)])/mean(esub[upper.tri(esub)]))), .4)
})

    

## M <- melt(boot[,,'groupB','C'])
## ggplot(M, aes(x=value))+geom_density() + facet_wrap(~X2)


          


