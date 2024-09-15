#Weibell Distribution: It is characterized by its shape parameter k > 0 and its scale parameter .
#Its probability density function is given by
setwd("~/Documents/Stats 102B")
source("bootsample.R")
source("bootstats.R")
k=1.5
lambda = 1
n = c(100,250) 
#and 250,
set.seed(0)
s_size100 = rweibull(n[1],shape=k,scale=lambda)
s_size250 = rweibull(n[2],shape=k,scale=lambda)
s_100_med = median(s_size100)
s_250_med = median(s_size250)
# Your results should be based on M = 100, 500 random samples and B = 5000
# bootstrap replicates.

M = c(100,500)
B = 5000
alpha=.05
zval=qnorm(p=alpha/2, lower.tail=FALSE)

xbar = function(w){median(w)}

#1. Normal bootstrap confidence interval.
set.seed(0)
size_100_CIs = data.frame(matrix(ncol = 2, nrow = 0))
size_250_CIs = data.frame(matrix(ncol = 2, nrow = 0))


set.seed(0)
for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size100, 100, replace = TRUE)
    boot.samples = bootsampling(samp,B)
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = median(samp)-zval*boot.x.bar$se
    upper_bound = median(samp)+zval*boot.x.bar$se
    cum_length = cum_length + upper_bound-lower_bound
    if(s_100_med <= upper_bound && s_100_med >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_100_CIs = rbind(size_100_CIs, cbind(param_included,cum_length))
}

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size250, 250, replace = TRUE)
    boot.samples = bootsampling(samp,B)
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = median(samp)-zval*boot.x.bar$se
    upper_bound = median(samp)+zval*boot.x.bar$se
    cum_length = cum_length+(upper_bound-lower_bound)
    if(s_250_med <= upper_bound && s_250_med >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_250_CIs = rbind(size_250_CIs, cbind(param_included,cum_length))
}

print(size_100_CIs[1,]/100)
print(size_100_CIs[2,]/500)
print(size_250_CIs[1,]/100)
print(size_250_CIs[2,]/500)

#2. Basic bootstrap confidence interval.
size_100_CIs = data.frame(matrix(ncol = 2, nrow = 0))
size_250_CIs = data.frame(matrix(ncol = 2, nrow = 0))

set.seed(0)

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size100, 100, replace = TRUE)
    boot.samples = bootsampling(samp,B)
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = 2*median(samp)-quantile(boot.x.bar$theta,probs=(1-alpha/2))
    upper_bound = 2*median(samp)-quantile(boot.x.bar$theta,probs=(alpha/2))
    cum_length = cum_length+upper_bound-lower_bound
    if(s_100_med <= upper_bound && s_100_med >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_100_CIs = rbind(size_100_CIs, cbind(param_included,cum_length))
}

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size250, 250, replace = TRUE)
    boot.samples = bootsampling(samp,B)    
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = 2*median(samp)-quantile(boot.x.bar$theta,probs=(1-alpha/2))
    upper_bound = 2*median(samp)-quantile(boot.x.bar$theta,probs=(alpha/2))
    cum_length = cum_length+(upper_bound-lower_bound)
    if(s_250_med <= upper_bound && s_250_med >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_250_CIs = rbind(size_250_CIs, cbind(param_included,cum_length))
}

print(size_100_CIs[1,]/100)
print(size_100_CIs[2,]/500)
print(size_250_CIs[1,]/100)
print(size_250_CIs[2,]/500)

#3. Percentile bootstrap confidence interval.
size_100_CIs = data.frame(matrix(ncol = 2, nrow = 0))
size_250_CIs = data.frame(matrix(ncol = 2, nrow = 0))
set.seed(0)
for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size100, 100, replace = TRUE)
    boot.samples = bootsampling(samp,B)
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = quantile(boot.x.bar$theta,probs=(alpha/2))
    upper_bound = quantile(boot.x.bar$theta,probs=(1-alpha/2))
    cum_length = cum_length+upper_bound-lower_bound
    if(s_100_med <= upper_bound && s_100_med >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_100_CIs = rbind(size_100_CIs, cbind(param_included,cum_length))
}

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size250, 250, replace = TRUE)
    boot.samples = bootsampling(samp,B) 
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = quantile(boot.x.bar$theta,probs=(alpha/2))
    upper_bound = quantile(boot.x.bar$theta,probs=(1-alpha/2))
    cum_length = cum_length+(upper_bound-lower_bound)
    if(s_250_med <= upper_bound & s_250_med >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_250_CIs = rbind(size_250_CIs, cbind(param_included,cum_length))
}

print(size_100_CIs[1,]/100)
print(size_100_CIs[2,]/500)
print(size_250_CIs[1,]/100)
print(size_250_CIs[2,]/500)

B=1000
size_100_CIs = data.frame(matrix(ncol = 2, nrow = 0))
size_250_CIs = data.frame(matrix(ncol = 2, nrow = 0))
set.seed(0)
for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size100, 100, replace = TRUE)
    boot.samples = bootsampling(samp,B)
    boot.x.bar = boot.stats(boot.samples,xbar)
    boot.samples.SE = rep(0,ncol(boot.samples))
    for (l in 1:ncol(boot.samples)) {
      # obtain bootstrap samples from original bootstrap sample l
      iter.boot.samples = bootsampling(boot.samples[,l],B)
      # get the std dev of the bootstrap sampling distribution based
      # on the l-th bootstrap sampe
      boot.samples.SE[l] = boot.stats(iter.boot.samples,xbar)$se
    }
    
    T = (boot.x.bar$theta-median(s_size100)) / boot.samples.SE
    qval = quantile(T,probs=c(alpha/2,1-alpha/2))
    
    lower_bound = median(samp)-qval[2]*boot.x.bar$se
    upper_bound = median(samp)-qval[1]*boot.x.bar$se
    cum_length = cum_length+(upper_bound-lower_bound)
    if(s_100_med <= upper_bound && s_100_med >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_100_CIs = rbind(size_100_CIs, cbind(param_included,cum_length))
}

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size250, 250, replace = TRUE)
    boot.samples = bootsampling(samp,B) 
    boot.x.bar = boot.stats(boot.samples,xbar)
    boot.samples.SE = rep(0,ncol(boot.samples))
    for (l in 1:ncol(boot.samples)) {
      # obtain bootstrap samples from original bootstrap sample l
      iter.boot.samples = bootsampling(boot.samples[,l],B)
      # get the std dev of the bootstrap sampling distribution based
      # on the l-th bootstrap sampe
      boot.samples.SE[l] = boot.stats(iter.boot.samples,xbar)$se
    }
    T = (boot.x.bar$theta-median(s_size250)) / boot.samples.SE
    qval = quantile(T,probs=c(alpha/2,1-alpha/2))
    
    lower_bound = median(samp)-qval[2]*boot.x.bar$se
    upper_bound = median(samp)-qval[1]*boot.x.bar$se
    cum_length = cum_length+(upper_bound-lower_bound)
    if(s_250_med <= upper_bound && s_250_med >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_250_CIs = rbind(size_250_CIs, cbind(param_included,cum_length))
}

print(size_100_CIs[1,]/100)
print(size_100_CIs[2,]/500)
print(size_250_CIs[1,]/100)
print(size_250_CIs[2,]/500)

## Part b)

k=1.5
lambda = 1
n = c(100,250) 
#and 250,

s_size100 = rweibull(n[1],shape=k,scale=lambda)
s_size250 = rweibull(n[2],shape=k,scale=lambda)
sd_100 = sd(s_size100)
sd_250 = sd(s_size250)
# Your results should be based on M = 100, 500 random samples and B = 5000
# bootstrap replicates.

M = c(100,500)
B = 5000
alpha=.05
zval=qnorm(p=alpha/2, lower.tail=FALSE)

xbar = function(w){sd(w)}

#1. Normal bootstrap confidence interval.
set.seed(0)
size_100_CIs = data.frame(matrix(ncol = 2, nrow = 0))
size_250_CIs = data.frame(matrix(ncol = 2, nrow = 0))

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size100, 100, replace = TRUE)
    boot.samples = bootsampling(samp,B)
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = sd(samp)-zval*boot.x.bar$se
    upper_bound = sd(samp)+zval*boot.x.bar$se
    cum_length = cum_length + upper_bound-lower_bound
    if(sd_100 <= upper_bound && sd_100 >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_100_CIs = rbind(size_100_CIs, cbind(param_included,cum_length))
}

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size250, 250, replace = TRUE)
    boot.samples = bootsampling(samp,B)
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = sd(samp)-zval*boot.x.bar$se
    upper_bound = sd(samp)+zval*boot.x.bar$se
    cum_length = cum_length+(upper_bound-lower_bound)
    if(sd_250 <= upper_bound && sd_250 >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_250_CIs = rbind(size_250_CIs, cbind(param_included,cum_length))
}

print(size_100_CIs[1,]/100)
print(size_100_CIs[2,]/500)
print(size_250_CIs[1,]/100)
print(size_250_CIs[2,]/500)


#2. Basic bootstrap confidence interval.
size_100_CIs = data.frame(matrix(ncol = 2, nrow = 0))
size_250_CIs = data.frame(matrix(ncol = 2, nrow = 0))

set.seed(0)

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size100, 100, replace = TRUE)
    boot.samples = bootsampling(samp,B)
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = 2*sd(samp)-quantile(boot.x.bar$theta,probs=(1-alpha/2))
    upper_bound = 2*sd(samp)-quantile(boot.x.bar$theta,probs=(alpha/2))
    cum_length = cum_length+upper_bound-lower_bound
    if(sd_100 <= upper_bound && sd_100 >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_100_CIs = rbind(size_100_CIs, cbind(param_included,cum_length))
}

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size250, 250, replace = TRUE)
    boot.samples = bootsampling(samp,B)    
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = 2*sd(samp)-quantile(boot.x.bar$theta,probs=(1-alpha/2))
    upper_bound = 2*sd(samp)-quantile(boot.x.bar$theta,probs=(alpha/2))
    cum_length = cum_length+(upper_bound-lower_bound)
    if(sd_250 <= upper_bound && sd_250 >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_250_CIs = rbind(size_250_CIs, cbind(param_included,cum_length))
}

print(size_100_CIs[1,]/100)
print(size_100_CIs[2,]/500)
print(size_250_CIs[1,]/100)
print(size_250_CIs[2,]/500)

#3. Percentile bootstrap confidence interval.
size_100_CIs = data.frame(matrix(ncol = 2, nrow = 0))
size_250_CIs = data.frame(matrix(ncol = 2, nrow = 0))
set.seed(0)
for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size100, 100, replace = TRUE)
    boot.samples = bootsampling(samp,B)
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = quantile(boot.x.bar$theta,probs=(alpha/2))
    upper_bound = quantile(boot.x.bar$theta,probs=(1-alpha/2))
    cum_length = cum_length+upper_bound-lower_bound
    if(sd_100 <= upper_bound && sd_100 >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_100_CIs = rbind(size_100_CIs, cbind(param_included,cum_length))
}

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size250, 250, replace = TRUE)
    boot.samples = bootsampling(samp,B) 
    boot.x.bar = boot.stats(boot.samples,xbar)
    lower_bound = quantile(boot.x.bar$theta,probs=(alpha/2))
    upper_bound = quantile(boot.x.bar$theta,probs=(1-alpha/2))
    cum_length = cum_length+(upper_bound-lower_bound)
    if(sd_250 <= upper_bound & sd_250 >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_250_CIs = rbind(size_250_CIs, cbind(param_included,cum_length))
}

print(size_100_CIs[1,]/100)
print(size_100_CIs[2,]/500)
print(size_250_CIs[1,]/100)
print(size_250_CIs[2,]/500)

#4. Bootstrap-t (studentized) confidence interval.
# I am using a smaller B number of replicates because my computer is not powerful enough
B=57
size_100_CIs = data.frame(matrix(ncol = 2, nrow = 0))
size_250_CIs = data.frame(matrix(ncol = 2, nrow = 0))
set.seed(0)
for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size100, 100, replace = TRUE)
    boot.samples = bootsampling(samp,B)
    boot.x.bar = boot.stats(boot.samples,xbar)
    boot.samples.SE = rep(0,ncol(boot.samples))
    for (l in 1:ncol(boot.samples)) {
      # obtain bootstrap samples from original bootstrap sample l
      iter.boot.samples = bootsampling(boot.samples[,l],B)
      # get the std dev of the bootstrap sampling distribution based
      # on the l-th bootstrap sampe
      boot.samples.SE[l] = boot.stats(iter.boot.samples,xbar)$se
    }
    
    T = (boot.x.bar$theta-sd(s_size100)) / boot.samples.SE
    qval = quantile(T,probs=c(alpha/2,1-alpha/2))
    
    lower_bound = sd(samp)-qval[2]*boot.x.bar$se
    upper_bound = sd(samp)-qval[1]*boot.x.bar$se
    cum_length = cum_length+(upper_bound-lower_bound)
    if(sd_100 <= upper_bound && sd_100 >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_100_CIs = rbind(size_100_CIs, cbind(param_included,cum_length))
}

for(m in 1:2){
  param_included = 0;
  cum_length = 0;
  for(i in 1:M[m]){
    samp = sample(s_size250, 250, replace = TRUE)
    boot.samples = bootsampling(samp,B) 
    boot.x.bar = boot.stats(boot.samples,xbar)
    boot.samples.SE = rep(0,ncol(boot.samples))
    for (l in 1:ncol(boot.samples)) {
      # obtain bootstrap samples from original bootstrap sample l
      iter.boot.samples = bootsampling(boot.samples[,l],B)
      # get the std dev of the bootstrap sampling distribution based
      # on the l-th bootstrap sampe
      boot.samples.SE[l] = boot.stats(iter.boot.samples,xbar)$se
    }
    T = (boot.x.bar$theta-sd(s_size250)) / boot.samples.SE
    qval = quantile(T,probs=c(alpha/2,1-alpha/2))
    
    lower_bound = sd(samp)-qval[2]*boot.x.bar$se
    upper_bound = sd(samp)-qval[1]*boot.x.bar$se
    cum_length = cum_length+(upper_bound-lower_bound)
    if(sd_250 <= upper_bound && sd_250 >= lower_bound){
      param_included = param_included + 1
    }
  }
  size_250_CIs = rbind(size_250_CIs, cbind(param_included,cum_length))
}

print(size_100_CIs[1,]/100)
print(size_100_CIs[2,]/500)
print(size_250_CIs[1,]/100)
print(size_250_CIs[2,]/500)

## problem 2
library(MASS)


#### a
n = 100000
p = 20
R = matrix( c( rep(0.99,p/2) , rep( 0.9 , p/2) ) , p , p )
diag(R) = 1
mean.vector= c(rep(0,20))
design = mvrnorm(n,mu=mean.vector,R)
error.term = rnorm(n, 0, 1.25)
beta_true = c(rep(3,20))
response =  design %*% beta_true  + error.term


##b

betas = lm(response~design+0)
# store and print the summary of the model
# including the beta coefficients
betas.summary=summary(betas)
betas.summary
sqrt(sum((c(betas.summary$coefficients[,1])-beta_true)^2))/sqrt(sum(beta_true^2))

#Gradient Descent

# here we define the step size
mystepsize=0.000001
# here we define the tolerance of the convergence criterion
mytol = 0.000000001 
# starting point
mystartpoint = c(rep(6,10), rep(0,10)) 

gradientDesc = function(y, X, startpoint, stepsize, 
                        conv_threshold, max_iter) {
  
  old.point = startpoint
  gradient = (t(X)%*%X%*%old.point - t(X)%*%y)
  new.point = old.point - stepsize * gradient
  
  old.value.function = t(y - X%*%old.point)%*%(y-X%*%old.point)
  
  converged = F
  iterations = 0
  
  while(converged == F) {
    ## Implement the gradient descent algorithm
    old.point = new.point
    gradient = t(X)%*%X%*%old.point - t(X)%*%y
    new.point = old.point - stepsize * gradient
    
    new.value.function = t(y - X%*%new.point)%*%(y-X%*%new.point)
    if( abs(old.value.function - new.value.function) <= conv_threshold) {
      converged = T
      break
    }
    
    data.output = data.frame(iteration = iterations,
                             old.value.function = old.value.function,
                             new.value.function = new.value.function,
                             old.point=old.point, new.point=new.point
    )
    
    if(exists("iters")) {
      iters <- rbind(iters, data.output)
    } else {
      iters = data.output
    }
    
    iterations = iterations + 1
    old.value.function = new.value.function
    if(iterations >= max_iter) break
    
  }
  return(list(converged = converged, 
              num_iterations = iterations, 
              old.value.function = old.value.function,
              new.value.function = new.value.function,
              coefs = new.point,
              iters = iters))
}

results = gradientDesc(response, design, mystartpoint, mystepsize, 
                       mytol, 30000)
print(results$coefs)
print(results$num_iterations)
sqrt(sum((results$coefs-beta_true)^2))/sqrt(sum(beta_true^2))


## c
epsilon = 0.4
gradientDesc_heavyball = function(y, X, startpoint, stepsize, 
                                  conv_threshold, max_iter, epsilon) {
  
  old.point = startpoint
  gradient = (t(X)%*%X%*%old.point - t(X)%*%y)
  new.point = old.point - stepsize * gradient
  
  old.value.function = t(y - X%*%old.point)%*%(y-X%*%old.point)
  
  converged = F
  iterations = 0
  
  while(converged == F) {
    ## Implement the gradient descent algorithm
    save_old_point = old.point
    old.point = new.point
    gradient = t(X)%*%X%*%old.point - t(X)%*%y
    new.point = old.point - stepsize * gradient +epsilon * (old.point - save_old_point)
    
    new.value.function = t(y - X%*%new.point)%*%(y-X%*%new.point)
    
    if( abs(old.value.function - new.value.function) <= conv_threshold) {
      converged = T
      break
    }
    
    data.output = data.frame(iteration = iterations,
                             old.value.function = old.value.function,
                             new.value.function = new.value.function,
                             old.point=old.point, new.point=new.point
    )
    
    if(exists("iters")) {
      iters <- rbind(iters, data.output)
    } else {
      iters = data.output
    }
    
    iterations = iterations + 1
    old.value.function = new.value.function
    if(iterations >= max_iter) break
    
  }
  return(list(converged = converged, 
              num_iterations = iterations, 
              old.value.function = old.value.function,
              new.value.function = new.value.function,
              coefs = new.point,
              iters = iters))
}

results = gradientDesc_heavyball(response, design, mystartpoint, mystepsize, 
                                 mytol, 30000,epsilon)

print(results$num_iterations)
print(results$coefs)
sqrt(sum((results$coefs-beta_true)^2))/sqrt(sum(beta_true^2))



## part d

mini_batch_size = 2000
mystepsize=0.000012

gradientDesc_heavyball_stochastic = function(y, X, startpoint, stepsize, 
                                             conv_threshold, max_iter, epsilon, Q) {
  
  z=cbind(y,X)
  
  # shuffle the data and select a mini batch
  
  shuffle.sample = z[sample(Q,replace=FALSE),]  
  
  ys=shuffle.sample[,1]
  Xs=shuffle.sample[,2:ncol(shuffle.sample)]
  
  old.point = startpoint
  gradient = (t(Xs)%*%Xs%*%old.point - t(Xs)%*%ys)
  
  new.point = old.point - stepsize * gradient
  
  old.value.function = t(ys - Xs%*%old.point)%*%(ys-Xs%*%old.point)
  
  converged = F
  iterations = 0
  
  while(converged == F) {
    ## Implement the gradient descent algorithm
    
    # shuffle the data and select a mini batch
    shuffle.sample = z[sample(Q,replace=FALSE),]  
    
    ys=shuffle.sample[,1]
    Xs=shuffle.sample[,2:ncol(shuffle.sample)]
    
    ## Implement the stochastic gradient descent algorithm
    save_old_point = old.point
    old.point = new.point
    
    gradient = t(Xs)%*%Xs%*%old.point - t(Xs)%*%ys
    
    new.point = old.point - stepsize * gradient+epsilon * (old.point - save_old_point)
    
    new.value.function = t(ys - Xs%*%new.point)%*%(ys-Xs%*%new.point)
    
    
    if( abs(old.value.function - new.value.function) <= conv_threshold) {
      converged = T
      break
    }
    
    data.output = data.frame(iteration = iterations,
                             old.value.function = old.value.function,
                             new.value.function = new.value.function,
                             old.point=old.point, new.point=new.point
    )
    
    if(exists("iters")) {
      iters <- rbind(iters, data.output)
    } else {
      iters = data.output
    }
    
    iterations = iterations + 1
    old.value.function = new.value.function
    if(iterations >= max_iter) break
    
  }
  return(list(converged = converged, 
              num_iterations = iterations, 
              old.value.function = old.value.function,
              new.value.function = new.value.function,
              coefs = new.point,
              iters = iters))
}


results = gradientDesc_heavyball_stochastic(response, design, mystartpoint, mystepsize, 
                                            mytol, 30000,epsilon, mini_batch_size)

print(results$num_iterations)
print(results$coefs)
sqrt(sum((results$coefs-beta_true)^2))/sqrt(sum(beta_true^2))












