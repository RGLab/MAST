pooledModel <- function(pi, observed, ncells, grad=FALSE, info=FALSE, loglik=TRUE){
  ## pi is parameter we are trying to estimate
  ## suppose there are n pools
  ## observed, n-vector: 1 if there was expression observed in the pool, 0 else
  ## ncells, n-vector: the number of cells in the pool
  #stopifnot(all(observed %in% c(0, 1)))
  negset <- ncells[observed==0]
  posset <- ncells[observed>0]
  gradval <- -sum(negset)/(1-pi) +         #negative pools
    sum(posset*(1-pi)^(posset-1)/ #positive pools
                                    (1- (1-pi)^posset))
  infoval <- (1-pi)^(-2)*sum(negset) +
    sum( (1-pi)^(2*posset -2) * posset^2/ (1 - (1-pi)^posset)^2 +
        (1-pi)^(posset-2)*(posset-1)*posset/(1-(1-pi)^posset))

  leqval <- sum(negset*log(1-pi)) + sum(log(1 - (1-pi)^posset))
  if(loglik) return(leqval)
  if(grad) return(gradval)
  if(info) return(infoval)
}




testSimulation <- function(){
piT <- c(.001, .005, .01, .05, .1, .5, 1)
ncells <- c(rep(1, 50), rep(10, 10), rep(100, 10))
nt <- 100
naive <- info <- res <- matrix(NA, ncol=length(piT), nrow=nt, dimnames=list(NULL, piT))
for(j in seq_along(piT)){
  for(i in 1:nt){
    obs <- rbinom(ncells, ncells, piT[j])
    tmp <- optimize(pooledModel, c(0, 1), maximum=TRUE, observed=obs, ncells=ncells)
    m <- tmp$maximum
    res[i, j] <- m
    info[i,j] <- pooledModel(m, obs, ncells, info=TRUE, loglik=FALSE)
    naive[i,j] <- if(any(obs[ncells==1]>0)){
      mean(obs[ncells==1]>0)} else if(any(obs[ncells==10]>0)){
        mean(obs[ncells==10]>0)/10
      } else{
        mean(obs[ncells==100]>0)/100
      }
  }
}
histogram(~log10(value)|X2, melt(res), scale='free')

apply(res, 2, mean)
bl <- apply(t((t(res)-piT)^2), 2, mean)
bl2 <- apply(t((t(naive)-piT)^2), 2, mean)
apply(1/sqrt(info), 2, mean)
}
