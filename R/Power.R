

simulateX <- function(n){
    rep(1, n)
}


simulate1Unit <- function(x, betaC, betaD, sigmaC){
    etaC <- x %*% betaC
    etaD <- expit(x %*% betaD)
    yC <- rnorm(length(etaC), sd=sigmaC)+etaC
    yD <- runif(length(etaD))<etaD
    list(y=yC*yD, x=x)
}
