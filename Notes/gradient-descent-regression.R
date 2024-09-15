# gradient descent for regression
# with constant step size
library(MASS)
set.seed(200)
# simulate data
n = 20 # sample size
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
mystepsize=0.066
# here we define the tolerance of the convergence criterion
mytol = 0.000000000000001 
# starting point
mystartpoint = c(1,2,3,5) 

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

results = gradientDesc(response, design, mystartpoint, 
                       mystepsize, mytol, 30000)

library(ggplot2)
ggplot(data = results$iters, mapping = aes(x = iteration, y = new.value.function))+
  geom_line()