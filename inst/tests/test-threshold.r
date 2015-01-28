set.seed(1234)
## nobs and ngenes
N <- p <- 400
## freq expression
pi <- runif(N)
## mean is reflected normal + 3*expression
mu1 <- abs(rnorm(N, 5, 3)+3*pi)
## variance of signal
sigma1 <- 1/rgamma(N, 3, scale=.75)
## mean/variance of noise (will be truncated at zero)
mu0 <- rnorm(N, 1, .5)
sigma0 <- rnorm(N, .5, .5)
mu0 <- ifelse(sigma0>0 & mu0>0, mu0, 0)
sigma0 <- ifelse(sigma0>0 & mu0>0, sigma0, 0)

## simulate one gene
sim1SCA <- function(N, mu1, mu0, sigma1, sigma0, pi){
    y1 <- rnorm(N, mu1, sigma1)
    y1 <- ifelse(y1<0, 0, y1)
    y0 <- rnorm(N, mu0, sigma0)
    y0 <- ifelse(y0<0, 0, y0)
    as.matrix(ifelse(runif(N)<pi, y1, y0))
}

## make matrix of genes
simMat <- aaply(1:p, 1, function(i){
    sim1SCA(N, mu1[i], mu0[i], sigma1[i], sigma0[i], pi[i])
})

tt <- thresholdSCRNACountMatrix(2^(simMat)-1, nbins=10, plot=T)
