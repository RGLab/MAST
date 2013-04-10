
zlm <- function(formula, data, ...){
  init <- lm(formula, data, method='model.frame', ...) #run lm initially just to get pull response vector
  data$pos <-   model.response(init)>0
  cont <- try(glm(formula, data, subset=pos, family='gaussian', ...))
  if(inherits(cont, 'try-error')){
    warning('Some factors were not present among the positive part')
    
  }                                     
  
  init[[1L]] <- (init[[1L]]>0)*1            #apparently the response goes first in the model.frame
  disc <- glm(formula, init, family='binomial')
  out <- list(cont=cont, disc=disc)
  class(out) <- 'zlm'
  out
}

summary.zlm <- function(out){
  summary(out$cont)
  summary(out$disc)
}

## Do LR test separately on continuous and discrete portions
test.zlm <- function(model, scope){
  tt <- try({
    if(missing(scope)){
      d.c <- drop1(model$cont, scope=scope, test='LRT')
} else{
 d.c <- drop1(model$cont, test='LRT')
}
  }, TRUE)
  
  if(inherits(tt, 'try-error')){
    d.c <- 0
  }

  tt <- try({
    if(missing(scope)){
      d.d <- drop1(model$disc, scope=scope, test='LRT')
} else{
 d.d <- drop1(model$disc, test='LRT')
}
  }, TRUE)
  
  if(inherits(tt, 'try-error')){
    d.d <- 0
  }

  
  
  d.n <- d.d+d.c
  
  d.n[,'Pr(>Chi)'] <- NA
  d.n
}

zlm.SingleCellAssay <- function(formula, sca, scope, ...){
    probeid <- getMapping(sca@mapping,"primerid")[[1]]
    measure <- getMapping(sca@mapping,"measurement")[[1]]

    m <- SingleCellAssay:::melt(sca)
    ssca <- split(m, m[,probeid], drop=TRUE)
    ll <- lapply(ssca, function(x){
      raw <- test.zlm(zlm(formula, x, ...), scope)
      out <- raw[-1,][,c('Df', 'LRT')]
      out <- rename(out, c('LRT'='lrstat'))
      out$p.value <- pchisq(out$lrstat, df=out$Df, lower.tail=FALSE)
      out
    })
    do.call(rbind, ll)
}
