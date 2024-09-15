# sampling with replacement functions

bootsampling <- function(x,boot.replicates=B){
  x = as.matrix(x)
  nx = nrow(x)
  bootsamples = replicate(boot.replicates,
                x[sample(nx,replace=TRUE),])
}

bootsampling.1 <- function(x,boot.replicates=B){
  x = as.matrix(x)
  nx = nrow(x)
  bootsamples = replicate(boot.replicates,
                          x[ceiling(runif(nx,0,1)*nx),])
}

