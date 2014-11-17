getX <- function(groups, minsize, N){
    gid <- LETTERS[seq_len(groups)]
    fidMin <- rep(gid, each=minsize)
    fid <- sample(gid[-1], size=N-length(fidMin), replace=TRUE)
    fid <- data.frame(group=c(fid, fidMin))
    X <- model.matrix(~group, data=fid)
    X
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

context("Bootstrap")
test_that("Only return coef works", {
    zzinit2 <- zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2, onlyCoef=TRUE)
    expect_that(zzinit2, is_a('array'))
    expect_equal(dim(zzinit2)[1], ncol(fd2))
})

test_that("Bootstrap", {
    zf <- zlm.SingleCellAssay( ~ Population*Stim.Condition, fd2)
    boot <- bootVcov1(zf, R=10)
    expect_is(boot, 'array')
    ## rep, genes, coef, comp
    expect_equal(dim(boot),c(10, dim(coef(zf, 'D')), 2))
})

context("Bootstrap consistency")
N <- 100
m <- 30
p <- 2
X <- getX(p, 30, N)
beta <- t(cbind(rnorm(m)+10, rnorm(m)))
Y <- simYs(m, X, beta, rho=1, sigma=1, p=seq(.1, 1, length.out=m))
