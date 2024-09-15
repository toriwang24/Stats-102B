# stochastic gradient descent for regression
# with backtracking line search for selecting the step size
library(MASS)
set.seed(200)
# simulate data
n = 1000000 # sample size
p = 4  # number of predictors
# create correlation matrix for regressors
R=matrix(c(1, 0.85, 0.55, 0.75,
           0.85, 1, 0.45, 0.65,
           0.55, 0.45, 1, 0.7,
           0.75, 0.65, 0.7, 1), 
         nrow = 4, ncol = 4, byrow = TRUE)
mean.vector = c(0,0,0,0)
# generate design matrix X
design.orig = mvrnorm(n,mu=mean.vector,R)
intercept = rep(1, n)
design = cbind(intercept,design.orig)
# generate error term
error.term = rnorm(n,0,1)
# generate beta
beta_true = c(1,1.5,2,-3,-2)
# generate response y
response = design%*%beta_true + error.term


# here we define the step size
mystepsize=5
# here we define the tolerance of the convergence criterion
mytol = 1e-15
# epsilon for backtracking line search
myepsilon = 0.5
# tau for backtracking line search
mytau = 0.5
# minibatch size
mymb = 0.001
# starting point
mystartpoint = c(1,2,3,5,6) 


SGD_BLS = function(y, X, startpoint, stepsize, conv_threshold, 
                   epsilon, tau, mb, max_iter) {

  mini_batch=ceiling(length(y)*mb)
  
  z=cbind(y,X)
  
  # shuffle the data and select a mini batch
  
  shuffle.sample = z[sample(mini_batch,replace=FALSE),]  
  
  ys=shuffle.sample[,1]
  Xs=shuffle.sample[,2:ncol(shuffle.sample)]
  
  old.point = startpoint
  gradient = (t(Xs)%*%Xs%*%old.point - t(Xs)%*%ys)
  
  # determine stepsize by backtracking line search
  
  while (t(ys - Xs%*%(old.point-stepsize*gradient))%*%(ys-Xs%*%(old.point-stepsize*gradient)) >
                t(ys - Xs%*%old.point)%*%(ys-Xs%*%old.point) - epsilon * stepsize * t(gradient) %*% gradient )            
     {
    stepsize = tau * stepsize
  }
  
  new.point = old.point - stepsize * gradient


  old.value.function = t(ys - Xs%*%old.point)%*%(ys-Xs%*%old.point)
  
  converged = F
  iterations = 0
  
  while(converged == F) {
    
    # shuffle the data and select a mini batch
    shuffle.sample = z[sample(mini_batch,replace=FALSE),]  
    
    ys=shuffle.sample[,1]
    Xs=shuffle.sample[,2:ncol(shuffle.sample)]
    
    ## Implement the stochastic gradient descent algorithm
    old.point = new.point
  
    gradient = t(Xs)%*%Xs%*%old.point - t(Xs)%*%ys
    
    # determine stepsize by backtracking line search
    
    while (t(ys - Xs%*%(old.point-stepsize*gradient))%*%(ys-Xs%*%(old.point-stepsize*gradient)) >
           t(ys - Xs%*%old.point)%*%(ys-Xs%*%old.point) - epsilon * stepsize * t(gradient) %*% gradient )            
    {
      stepsize = tau * stepsize
    }
    
    new.point = old.point - stepsize * gradient
    
    new.value.function = t(ys - Xs%*%new.point)%*%(ys-Xs%*%new.point)
    
    
    if( abs(old.value.function - new.value.function) <= conv_threshold) {
      converged = T
    }
    
    data.output = data.frame(iteration = iterations,
                       old.value.function = old.value.function,
                       new.value.function = new.value.function,
                       old.point=old.point, new.point=new.point,
                       stepsize = stepsize
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
              stepsize = stepsize,
              iters = iters))
}

start = Sys.time()
results = SGD_BLS(response, design, mystartpoint, mystepsize, mytol, 
                       myepsilon, mytau, mymb, 30000)


print(Sys.time()-start)

library(ggplot2)

quartz(width=12,height=8)
par(mfrow=c(1,1))

ggplot(data = results$iters, mapping = aes(x = iteration, y = new.value.function))+
  geom_line()

quartz(width=12,height=8)
par(mfrow=c(1,1))

ggplot(data = results$iters, mapping = aes(x = iteration, y = stepsize))+
  geom_line()

start = Sys.time()
summary(lm(response ~ design.orig))
print(Sys.time()-start)