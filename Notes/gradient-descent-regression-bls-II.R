# gradient descent for regression
# with backtracking line search for selecting the step size
library(MASS)
set.seed(200)
# simulate data
n = 100 # sample size
p = 3 # number of predictors
# create correlation matrix for regressors
R=matrix(c(1, 0.85, 0.55,
           0.85, 1, 0.45,
           0.55, 0.45, 1), 
         nrow = 3, ncol = 3, byrow = TRUE)
mean.vector = c(0,0,0)
# generate design matrix X
design.orig = mvrnorm(n,mu=mean.vector,R)
intercept = rep(1, n)
design = cbind(intercept,design.orig)
# generate error term
error.term = rnorm(n,0,1)
# generate beta
beta_true = c(1,1.5,2,-3)
# generate response y
response = design%*%beta_true + error.term


# here we define the step size
mystepsize=5
# here we define the tolerance of the convergence criterion
mytol = 0.000000000000001 
# epsilon for backtracking line search
myepsilon = 0.5
# tau for backtracking line search
mytau = 0.5
# starting point
mystartpoint = c(1,2,3,5) 

gradientDescBLS = function(y, X, startpoint, stepsize.original, 
                        conv_threshold, epsilon, tau, max_iter) {
  
  old.point = startpoint
  gradient = (t(X)%*%X%*%old.point - t(X)%*%y)

  
  # determine stepsize by backtracking line search
  # reset stepsize to the original value at every iteration
  stepsize=stepsize.original
  
  while (t(y - X%*%(old.point-stepsize*gradient))%*%(y-X%*%(old.point-stepsize*gradient)) >
         t(y - X%*%old.point)%*%(y-X%*%old.point) - epsilon * stepsize * t(gradient) %*% gradient )            
  {
    stepsize = tau * stepsize
  }

  new.point = old.point - stepsize * gradient

  old.value.function = t(y - X%*%old.point)%*%(y-X%*%old.point)
  
  converged = F
  iterations = 0
  
  while(converged == F) {
    ## Implement the gradient descent algorithm
    old.point = new.point
  
    gradient = t(X)%*%X%*%old.point - t(X)%*%y
    
    # determine stepsize by backtracking line search
    # reset stepsize to the original value at every iteration
    stepsize=stepsize.original
    
    while (t(y - X%*%(old.point-stepsize*gradient))%*%(y-X%*%(old.point-stepsize*gradient)) >
           t(y - X%*%old.point)%*%(y-X%*%old.point) - epsilon * stepsize * t(gradient) %*% gradient )            
    {
      stepsize = tau * stepsize
    }
    
    new.point = old.point - stepsize * gradient
    
    new.value.function = t(y - X%*%new.point)%*%(y-X%*%new.point)
    
    
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

results = gradientDescBLS(response, design, mystartpoint, mystepsize, mytol, 
                       myepsilon, mytau, 30000)



library(ggplot2)

quartz(width=12,height=8)
par(mfrow=c(1,1))

ggplot(data = results$iters, mapping = aes(x = iteration, y = new.value.function))+
  geom_line()

quartz(width=12,height=8)
par(mfrow=c(1,1))

ggplot(data = results$iters, mapping = aes(x = iteration, y = stepsize))+
  geom_line()