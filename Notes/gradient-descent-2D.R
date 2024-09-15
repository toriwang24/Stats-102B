# gradient descent for 2D differentiable functions

#define objective function to be optimized
objective.function = function(x) {
            3*x[1]^2 + 0.5*x[2]^2 + 2*x[1]*x[2] 
    }

# define the derivative for a 2D objective function
derivative = function(x) {
  c(6*x[1]+2*x[2], x[2]+2*x[1])
}

# here we define the step size
mystepsize=0.29
# here we define the tolerance of the convergence criterion
mytol = 0.000000001 
# starting point
mystartpoint = c(1,2) 

gradientDesc = function(obj.function, startpoint, 
                        stepsize, conv_threshold, max_iter) {
  
  old.point = startpoint
  gradient = derivative(old.point)
  new.point = c(old.point[1] - stepsize*gradient[1], 
                old.point[2] - stepsize*gradient[2])
  

  old.value.function = obj.function(new.point)
  
  converged = F
  iterations = 0
  
  while(converged == F) {
    ## Implement the gradient descent algorithm
    old.point = new.point
    print(old.point)
  
    gradient = derivative(old.point)
    new.point = c(old.point[1] - stepsize*gradient[1], 
                  old.point[2] - stepsize*gradient[2])
    
    new.value.function = obj.function(new.point)
    
    
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

results = gradientDesc(objective.function, mystartpoint, mystepsize, 
                       mytol, 30000)

library(ggplot2)
ggplot(data = results$iters, mapping = aes(x = iteration, y = new.value.function))+
  geom_line()